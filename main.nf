#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.peak_depth = 50 
params.bam = null 
params.gtf = null
params.homopol = null
params.outdir = "results"
params.refgenome = null

if( !params.homopol && !params.refgenome ) { 
    error "Please provide either --homopol or --refgenome (to find homopolymers)"
} else if (params.refgenome) {
    params.homopol = "${params.outdir}/homopolymers.bed"
}

if( !params.bam || !params.gtf ) { error "Please provide --bam and --gtf" }

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

    peaks_result = paraclu(bam_ch)

    peak_connections(
        peaks_result.join(bam_ch),
        file(params.gtf)
    )

    plot_qc(
        peak_connections.out.bam_id
            .merge(peak_connections.out.connections)
            .merge(peak_connections.out.annotations)
    )
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

    output:
    val bam_id,                                     emit: bam_id
    path "${bam_id}_peak_pairwise_connections.tsv", emit: connections
    path "${bam_id}_peak_annotation.tsv",           emit: annotations
    path "${bam_id}_peak_intervals.gtf",            emit: gtf

    script:
    """
    python ${workflow.projectDir}/bin/connect_peaks.py \
        ${peaks_bed} \
        ${gtf_file} \
        ${bam} \
        ${bam_id}_ \
        ${params.homopol}
    """
}