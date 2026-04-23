#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.peak_depth = 25
params.outdir = "results"
params.strandedness = "default" // if not default we can put "reverse" to induce behaviour
params.libtype = "10x" // or "argentag", etc. for adjusting peak filtering
params.trans_spliced = null // adjust peak filtering for trans leader (ctenophore) case

params.chr = null // default runs all chromosomes
params.bam = null 
params.gtf = null
params.homopol = null

params.sj = null // optional: path to STAR or otherwise derived SJ.out.tab format
params.refgenome = null // optional

params.ip_bin = "/users/mirimia/jchamberlin/software/intronProspector/bin/intronProspector"
params.scaf_keyword = null // e.g. 'scaf' or 'scaffold': contigs whose name contains this string are batched into one job


if( !params.homopol && !params.refgenome ) {
    error "Please provide either --homopol or --refgenome (to find homopolymers)"
} else if (params.refgenome) {
    params.homopol = "${params.outdir}/homopolymers.bed"
}

if( !params.gtf ) { error "Please provide --gtf" }
if( !params.bam ) { error "Please provide --bam" }

workflow {
    bam_ch = Channel
        .fromPath("${params.bam}{,.bai}")
        .collect()
        .map { files ->
            def bam = files.find { it.name.endsWith('.bam') }
            def bai = files.find { it.name.endsWith('.bai') }
            if (!bai) error "Index file (.bai) not found for ${bam}"
            return tuple(bam.baseName, bam, bai)
        }

    if (params.sj) {
        sj_ch = bam_ch.map { bam_id, bam, bai -> tuple(bam_id, file(params.sj)) }
    } else {
        sj_ch = extract_sj(bam_ch, params.gtf ?: [])
    }

    peaks_result = paraclu(bam_ch)

    peak_connections(
        peaks_result.join(bam_ch),
        file(params.gtf ?: []),
        params.libtype
    )

    plot_qc(
        peak_connections.out.connections
            .join(peak_connections.out.annotations)
    )

    chrs_ch = get_chromosomes(bam_ch, file(params.gtf))
        .flatMap { bam_id, chrs ->
            def all_chrs = chrs.trim().split('\n').toList()
            if (params.scaf_keyword) {
                def batched = all_chrs.findAll { it.contains(params.scaf_keyword) }
                def solo    = all_chrs.findAll { !it.contains(params.scaf_keyword) }
                def items   = solo.collect { [bam_id, it] }
                if (batched) items << [bam_id, batched.join(' ')]
                return items
            } else {
                return all_chrs.collect { [bam_id, it] }
            }
        }

    if (params.chr) {
        chrs_ch = chrs_ch.filter { bam_id, chr -> chr == params.chr }
    }

    cadena_inputs = bam_ch
        .join(peak_connections.out.gtf)
        .join(peak_connections.out.connections)
        .join(sj_ch)
        .combine(chrs_ch, by: 0)

    cadena(cadena_inputs, params.strandedness, params.libtype)

    if (!params.chr) {
        merge_cadena(
            cadena.out.cadena_out
                .groupTuple()
        )
    }
}


process merge_cadena {
    tag "${bam_id}"
    publishDir "${params.outdir}/cadena", mode: 'copy'

    input:
    tuple val(bam_id), path(cadena_files)

    output:
    tuple val(bam_id), path("${bam_id}_cadena_reads.tsv.gz"),   emit: reads
    tuple val(bam_id), path("${bam_id}_cadena_pro_tss.ids.tsv"), emit: intervals
    tuple val(bam_id), path("${bam_id}_cadena_sjs.ids.tsv"),     emit: sjs

    script:
    """
    python ${workflow.projectDir}/bin/merge_cadena.py ${bam_id}
    """
}

process cadena { // aka assign reads to TSS-SJs-TES, based on hilgers-lab/laser
    tag "${bam_id}:${chr.split(' ')[0]}${chr.contains(' ') ? '+' + (chr.split(' ').size() - 1) + ' more' : ''}"
    publishDir "${params.outdir}/cadena", mode: 'copy', enabled: params.chr as boolean

    input:
    tuple val(bam_id), path(bam), path(bai), path(custom_intervals_gtf), path(connections), path(sj_file), val(chr)
    val strandedness
    val libtype

    output:
    tuple val(bam_id), path("${bam_id}_*_cadena_*"), optional: true, emit: cadena_out

    script:
    """
    if [[ "${chr}" == *" "* ]]; then
        Rscript ${workflow.projectDir}/bin/cadena.R \
            ${bam} \
            ${custom_intervals_gtf} \
            ${sj_file} \
            ${strandedness} \
            ${chr} \
            ${bam_id}_scaf_cadena_ \
            ${libtype} \
            ${chr}
    else
        Rscript ${workflow.projectDir}/bin/cadena.R \
            ${bam} \
            ${custom_intervals_gtf} \
            ${sj_file} \
            ${strandedness} \
            ${chr} \
            ${bam_id}_${chr}_cadena_ \
            ${libtype}
    fi
    """
}


process plot_qc {
    cache false
    tag "${bam_id}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(bam_id), path(connections), path(annotations)

    output:
    path "${bam_id}_qc_plots.pdf"

    script:
    """
    python ${workflow.projectDir}/bin/plot_qc.py \
        ${connections} \
        ${annotations} \
        ${bam_id}_qc_plots.pdf
    """
}

process paraclu {
    tag "${bam_id}"
    publishDir "${params.outdir}/peaks", mode: 'copy'

    input:
    tuple val(bam_id), path(bam), path(bai)

    output:
    tuple val(bam_id), path("${bam_id}.d${params.peak_depth}.peaks.bed")

    shell:
    '''
    PARACLU="/users/mirimia/jchamberlin/software/paraclu/paraclu"

    python !{workflow.projectDir}/bin/prep_paraclu.py !{bam} ./

    pas_input=$(ls *pas_count_input_1based.tsv)
    tss_input=$(ls *tss_count_input_1based.tsv)

    $PARACLU !{params.peak_depth} "$pas_input" > pas.tmp
    $PARACLU !{params.peak_depth} "$tss_input" > tss.tmp

    ${PARACLU}-cut.sh pas.tmp > pas.tsv
    ${PARACLU}-cut.sh tss.tmp > tss.tsv

    awk -F'\t' 'BEGIN{OFS="\t"} {
        start=$3-1
        print $1, start, $4, "pas_" $6 "_" $1 ":" start, $6, $2
    }' pas.tsv > pas.bed

    awk -F'\t' 'BEGIN{OFS="\t"} {
        start=$3-1
        print $1, start, $4, "tss_" $6 "_" $1 ":" start, $6, $2
    }' tss.tsv > tss.bed

    cat pas.bed tss.bed > !{bam_id}.d!{params.peak_depth}.peaks.bed
    '''
}


process peak_connections {
    tag "${bam_id}"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(bam_id), path(peaks_bed), path(bam), path(bai)
    path gtf_file
    val libtype

    output:
    val bam_id,                                                          emit: bam_id
    tuple val(bam_id), path("${bam_id}_peak_pairwise_connections.tsv"), emit: connections
    tuple val(bam_id), path("${bam_id}_peak_annotation.tsv"),           emit: annotations
    tuple val(bam_id), path("${bam_id}_peak_intervals.gtf"),            emit: gtf


    script:
    """
    python ${workflow.projectDir}/bin/connect_peaks.py \
        ${peaks_bed} \
        ${gtf_file} \
        ${bam} \
        ${bam_id}_ \
        ${params.homopol} \
        ${params.trans_spliced ?: 'false'} \
        ${libtype}
    """
}

process get_chromosomes {
    input:
    tuple val(bam_id), path(bam), path(bai)
    path gtf

    output:
    tuple val(bam_id), stdout

    script:
    """
    samtools idxstats ${bam} | awk '\$3 > 0 {print \$1}' | grep -v '\\*' | sort > bam_chrs.txt
    grep -v '^#' ${gtf} | cut -f1 | sort -u > gtf_chrs.txt
    comm -12 bam_chrs.txt gtf_chrs.txt
    """
}

// generate SJs from GTF or BAM if not provided by user
process extract_sj {
    tag "${bam_id}"
    publishDir "${params.outdir}/prep", mode: 'copy'
    input:
    tuple val(bam_id), path(bam), path(bai)
    val gtf_file

    output:
    tuple val(bam_id), path("${bam_id}_SJ.out.tab")

    script:
    def gtf_arg = gtf_file ? "--gtf ${gtf_file}" : "--bam ${bam} --genome ${params.refgenome}"
    """
    Rscript ${workflow.projectDir}/bin/extract_SJs.R \
        ${gtf_arg} \
        --out ${bam_id}_SJ.out.tab
    """
}