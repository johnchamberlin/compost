#' Title
#'
#' @param pathSJ.out: object of junctions from STAR
#' @param customAnnotPath : peak-to-peak transcript intervals
#'
#' @return
#' @export
#'
#' @examples

custom_reference_junctions <- function(pathSJ.out, customAnnotPath){
    # peak to peak transcript intervals
    custom_annotation = rtracklayer::import.gff(customAnnotPath)

    gene_ranges <- custom_annotation %>%
        as.data.frame() %>%
        group_by(seqnames, strand, gene_id) %>%
        summarise(
            start = min(start),
            end = max(end),
            .groups = "drop"
        ) %>%
        mutate(
            source = "custom",
            type = "gene"
        )

    # Convert back to GRanges
    custom_genes <- GenomicRanges::makeGRangesFromDataFrame(
        gene_ranges,
        keep.extra.columns = TRUE
    )

    star_sjs = data.table::fread(sj_path,col.names = c("chromosome",
                                                        "start",
                                                        "end",
                                                        "strand",
                                                        "intron_motif",
                                                        "annotated",
                                                        "unique_reads",
                                                        "multi_mapping_reads",
                                                        "maximum_overhang")) %>% 
            mutate(strand = recode(strand, `0` = "*", `1` = "+",`2` = "-")) %>%
            GenomicRanges::makeGRangesFromDataFrame()

    # intersect gene intervals with sjs
    hits <- GenomicRanges::findOverlaps(star_sjs, custom_genes)
    star_sjs$gene_id <- NA_character_
    star_sjs$gene_id[GenomicRanges::queryHits(hits)] <- custom_genes$gene_id[GenomicRanges::subjectHits(hits)]
    star_sjs = star_sjs %>%  mutate(juncID = paste0(seqnames,":",start,"-",end))
    return(star_sjs)
}
