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

    prep_results = preparaclu(bam_ch)
    peaks_result = paraclu_peaks(prep_results)

    peak_connections(
        peaks_result.join(bam_ch),
        file(params.gtf)
    )
}

process preparaclu {
    tag "${bam_id}"
    publishDir "${params.outdir}/prep", mode: 'copy'
    
    input:
    tuple val(bam_id), path(bam), path(bai)

    output:
    tuple val(bam_id), path("prep_out/*")

    script:
    """
    mkdir -p prep_out
    python ${workflow.projectDir}/bin/prep_paraclu.py ${bam} prep_out/
    """
}

process paraclu_peaks {
    tag "${bam_id}"
    publishDir "${params.outdir}/peaks", mode: 'copy'

    input:
    tuple val(bam_id), path(prep_files)

    output:
    tuple val(bam_id), path("paraclu.d${params.peak_depth}.peaks.bed")

    shell:
    '''
    PARACLU="/users/mirimia/jchamberlin/software/paraclu/paraclu"
    
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

    cat pas.bed tss.bed > paraclu.d!{params.peak_depth}.peaks.bed
    '''
}

process peak_connections {
    tag "${bam_id}"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(bam_id), path(peaks_bed), path(bam), path(bai)
    path gtf_file

    output:
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