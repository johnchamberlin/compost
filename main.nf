#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.peak_depth = 50 
params.bam = null 
params.gtf = null
params.outdir = "results"

// Validation
if( !params.bam || !params.gtf ) { error "Please provide --bam and --gtf" }

workflow {
    // 1. Channel setup (BAM + Index)
    bam_ch = Channel
        .fromPath("${params.bam}{,.bai}")
        .collect() 
        .map { files -> 
            def bam = files.find { it.name.endsWith('.bam') }
            def bai = files.find { it.name.endsWith('.bai') }
            if (!bai) error "Index file (.bai) not found for ${bam}"
            return tuple(bam.baseName, bam, bai)
        }

    // 2. Execution Flow
    prep_results = preparaclu(bam_ch)
    peaks_result = paraclu_peaks(prep_results)
    
    // 3. Final Step: Peak Connections
    // We pass the gtf file using the file() helper to ensure it's treated as a path
    peak_connections(peaks_result, file(params.gtf))
}

process preparaclu {
    tag "${bam_id}"
    publishDir "${params.outdir}/prep", mode: 'copy'
    
    input:
    tuple val(bam_id), path(bam), path(bai)

    output:
    // Capture the files directly rather than just the directory for better staging
    tuple val(bam_id), path(bam), path(bai), path("prep_out/*")

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
    // 'prep_files' now receives the list of files from the previous output
    tuple val(bam_id), path(bam), path(bai), path(prep_files)

    output:
    tuple val(bam_id), path(bam), path(bai), path("paraclu.d${params.peak_depth}.peaks.bed")

shell:
    '''
    PARACLU="/users/mirimia/jchamberlin/software/paraclu/paraclu"
    
    pas_input=$(ls *pas_count_input_1based.tsv)
    tss_input=$(ls *tss_count_input_1based.tsv)

    # 1. Run Paraclu
    $PARACLU !{params.peak_depth} "$pas_input" > pas.tmp
    $PARACLU !{params.peak_depth} "$tss_input" > tss.tmp

    # 2. Run Paraclu-cut
    ${PARACLU}-cut.sh pas.tmp > pas.tsv
    ${PARACLU}-cut.sh tss.tmp > tss.tsv

    # 3. Format as BED6 
    # Paraclu columns: 1:Chrom, 2:Strand, 3:Start, 4:End, 5:Score
    # BED6 target:    1:Chrom, 2:Start, 3:End, 4:Name, 5:Score, 6:Strand
    
    # We also subtract 1 from the Start ($3) because BED is 0-based
    awk 'BEGIN{OFS="\t"} {print $1, $3-1, $4, "pas_" $5 "_" NR, $5, $2}' pas.tsv > pas.bed
    awk 'BEGIN{OFS="\t"} {print $1, $3-1, $4, "tss_" $5 "_" NR, $5, $2}' tss.tsv > tss.bed

    # 4. Concatenate the properly formatted BED files
    cat pas.bed tss.bed > paraclu.d!{params.peak_depth}.peaks.bed
    '''
}


process peak_connections {
    tag "${bam_id}"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    // This tuple matches the output of paraclu_peaks
    tuple val(bam_id), path(bam), path(bai), path(peaks_bed)
    path gtf_file

    output:
    path "${bam_id}_peak_pairwise_connections.tsv", emit: connections
    path "${bam_id}_peak_annotation.tsv",           emit: annotations

    script:
    // Using the bash variable ${bam_id} for the prefix argument
    """
    python ${workflow.projectDir}/bin/connect_peaks.py \
        ${peaks_bed} \
        ${gtf_file} \
        ${bam} \
        ${bam_id}_
    """
}