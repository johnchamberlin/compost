#' Count full length reads
#'
#' @param alignments: alignments file from minimap2
#' @param linksDatabase: created with prepare_links_database
#'
#' @return
#' @export
#'
#' @examples
countLinks <- function(alignments, linksDatabase) {

  allReads <- data.frame(name = alignments$name)

  # make single nt starts
  startsAlignments <- prepareForCountStarts(alignments, 1)
  startAlignments <-
    readTSSassignment(startsAlignments, linksDatabase$TSSCoordinate.base)
  startAlignments <-
    startAlignments %>% as.data.frame(.) %>% dplyr::select(name, promoter_id)
  # make single nt ends
  endsAlignments <- prepareForCountEnds(alignments, 1)
  endsAlignments <-
    readTESassignment(endsAlignments, linksDatabase$TESCoordinate.base)
  endsAlignments <-
    endsAlignments %>% as.data.frame(.) %>% dplyr::select(name, tes_id)

  # join from all reads so unassigned ends are kept as non_peak
  pairsTested <- allReads %>%
    left_join(endsAlignments, by = "name") %>%
    left_join(startAlignments, by = "name") %>%
    dplyr::mutate(
      promoter_id = as.character(ifelse(is.na(promoter_id), "non_peak", promoter_id)),
      tes_id      = as.character(ifelse(is.na(tes_id),      "non_peak", tes_id))
    ) %>%
    dplyr::filter(
      tes_id == "non_peak" | promoter_id == "non_peak" |
      gsub("\\:.*", "", tes_id) == gsub("\\:.*", "", promoter_id)
    )

  # get gene_id from tes_id where available, fall back to promoter_id
  gene_from_tes <- linksDatabase$pairDataBase %>%
    dplyr::distinct(tes_id, .keep_all = TRUE) %>%
    dplyr::select(gene_id, tes_id)
  gene_from_promoter <- linksDatabase$pairDataBase %>%
    dplyr::distinct(promoter_id, .keep_all = TRUE) %>%
    dplyr::select(gene_id, promoter_id)

  pairsTested <- pairsTested %>%
    left_join(gene_from_tes, by = "tes_id") %>%
    left_join(gene_from_promoter, by = "promoter_id", suffix = c("", ".from_promoter")) %>%
    dplyr::mutate(gene_id = ifelse(is.na(gene_id), gene_id.from_promoter, gene_id)) %>%
    dplyr::select(-gene_id.from_promoter) %>%
    mutate(pairs_id = paste0(
      gene_id, ":",
      gsub(".*:", "", promoter_id), ":",
      gsub(".*:", "", tes_id)
    ))

  countedPairs <- pairsTested %>% group_by(pairs_id) %>% tally()
  countedPairsFinal <-
    left_join(
      countedPairs,
      pairsTested %>% dplyr::distinct(pairs_id, .keep_all = TRUE) %>% dplyr::select(gene_id, pairs_id),
      by = "pairs_id"
    )
  # normalize expression in CPMs
  if (nrow(countedPairsFinal) == 0) {
    countedPairsFinal$pairs_cpm <- numeric(0)
    return(list(countedPairsFinal = countedPairsFinal, pairsTested = pairsTested))
  }
  dd <- edgeR::DGEList(counts = countedPairsFinal$n)
  dge <- edgeR::calcNormFactors(dd)
  countedPairsFinal$pairs_cpm <- as.numeric(edgeR::cpm(dge))
  res <- list()
  res$countedPairsFinal <- countedPairsFinal
  res$pairsTested <- pairsTested
  return(res)
}
#' prepare counts trimming the reads to their most 5' into a 1 nt window for counting
#'
#' @param x alignments: alignments file from minimap2
#' @param window window for trimming the reads to their most 5'/3'
#'
#' @return
#' @export
#'
#' @examples
prepareForCountStarts <- function(x, window) {
  alignments <- x
  pos <- alignments[alignments@strand == "+",]
  neg <- alignments[alignments@strand == "-",]
  GenomicRanges::end(pos) <- GenomicRanges::start(pos) + window
  GenomicRanges::start(neg) <- GenomicRanges::end(neg) - window
  shortstarts <- c(pos, neg)
  return(shortstarts)
}
#' prepare counts trimming the reads to their most 3' into a 1 nt window for counting
#'
#' @param x  alignments: alignments file from minimap2
#' @param window window for trimming the reads to their most 5'/3'
#' @return
#' @export
#'
#' @examples
prepareForCountEnds <- function(x, window) {
  alignments <- x
  pos <- alignments[alignments@strand == "+",]
  neg <- alignments[alignments@strand == "-",]
  GenomicRanges::start(pos) <- GenomicRanges::end(pos) - window
  GenomicRanges::end(neg) <- GenomicRanges::start(neg) + window
  shortstarts <- c(pos, neg)
  return(shortstarts)
}
#' Title
#'
#' @param startsAlignments  trimmed reads  output from prepareForCountStarts
#' @param TSSCoordinate.base reference database for TSS
#'
#' @return
#' @export
#'
#' @examples
readTSSassignment <- function(startsAlignments, TSSCoordinate.base) {
  tssDb <- TSSCoordinate.base
  ovlps <- GenomicRanges::findOverlaps(startsAlignments , tssDb , maxgap = 50)
  promoterIds <- tssDb[subjectHits(ovlps),]$count
  promoterStarts <- GenomicRanges::start(tssDb[subjectHits(ovlps),])
  startsAlignments2 <- startsAlignments[queryHits(ovlps),]
  startsAlignments2$promoter_id <- promoterIds
  startsAlignments2$promoterStarts <- promoterStarts
  startsAlignments2$dist2Assignment <-
    abs(GenomicRanges::start(startsAlignments2) - startsAlignments2$promoterStarts)
  startsAlignments2 <-
    startsAlignments2[order(startsAlignments2$dist2Assignment, decreasing = FALSE),]
  startsAlignments2 <-
    startsAlignments2[!duplicated(startsAlignments2$name),]
  return(startsAlignments2)
}
#' Title
#'
#' @param endsAlignements trimmed reads  output from prepareForCountEnds
#' @param TESCoordinate.base reference database for TES
#'
#' @return
#' @export
#'
#' @examples
readTESassignment <- function(endsAlignements, TESCoordinate.base) {
  tesDb <- TESCoordinate.base
  ovlps <- findOverlaps(endsAlignements , tesDb , maxgap = 150)
  tesIds <- tesDb[subjectHits(ovlps),]$count
  endStarts <- GenomicRanges::start(tesDb[subjectHits(ovlps),])
  endsAlignments2 <- endsAlignements[queryHits(ovlps),]
  endsAlignments2$tes_id <- tesIds
  endsAlignments2$endStarts <- endStarts
  endsAlignments2$dist2Assignment <-
    abs(GenomicRanges::start(endsAlignments2) - endsAlignments2$endStarts)
  endsAlignments2 <-
    endsAlignments2[order(endsAlignments2$dist2Assignment, decreasing = FALSE),]
  endsAlignments2 <-
    endsAlignments2[!duplicated(endsAlignments2$name),]
  return(endsAlignments2)
}

#' Identify full length reads for downstream analysis
#'
#' @param linksDatabase: created with prepare_links_database
#' @param alignmentsFile : bam file alignments object
#' @param tss.ntwindow : window for the database assignment for full length at thes TSS
#' @param tes.ntwindow : window for the database assignment for full length at thes TES
#' @param reverse_strandedness : strandedness reversal for opposite strand alignmen
#'
#' @return
#' @export
#'
#' @examples

get_sailor_custom <- function(linksDatabase, alignmentsFile, tss.ntwindow, tes.ntwindow, reverse_strandedness=TRUE) {
  alignmentsFile <- GenomicRanges::GRanges(alignmentsFile)
  if (reverse_strandedness){
    message("Reversing strand information")
    current_strand <- strand(alignmentsFile)
    flipped_strand <- ifelse(current_strand == "+", "-", 
                            ifelse(current_strand == "-", "+", "*"))
    strand(alignmentsFile) <- flipped_strand  
    }
  alignmentsFile$name <- names(alignmentsFile)
  names(alignmentsFile) <- NULL
  countsLinks <- countLinks(alignmentsFile, linksDatabase)
  return(countsLinks)
  }


