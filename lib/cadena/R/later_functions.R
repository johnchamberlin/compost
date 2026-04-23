# currently not integrated in the cadena package

#' estimatePromoterDominance
#'
#' @param LATER
#' @param IsoformDatabase
#' @param method
#'
#' @return
#' @export estimatePromoterDominance
#' @import GenomicFeatures GenomicAlignments S4Vectors dplyr
estimatePromoterDominance <- function(LATER, IsoformDatabase, method) {
  dominance <- calculatePromoterDominance(LATER, IsoformDatabase)
  transcriptional_bias <- estimateTranscriptionalBias(dominance, method)
  dominance(LATER) <- dominance
  #result(LATER) <- transcriptional_bias %>% tibble::as_tibble()
  stats(LATER) <- transcriptional_bias %>% tibble::as_tibble()
  return(LATER)
}


#' Caculate promoter dominance
#'
#' @param countData counts data from 5'-3' links
#' @param pairsDataBase database of 5'-3' link isoforms
#'
#' @return
#' @export
#' @import GenomicFeatures GenomicAlignments S4Vectors dplyr
#' @examples
calculatePromoterDominance <- function(LATER, IsoformDatabase) {
  countData <- isoformCounts(LATER)
  pairsDataBase <- showLinks(IsoformDatabase)
  # Assign count data to their respective 3'end and 5'ends.
  annotPairsExp <-
    left_join(
      countData %>% dplyr::select(!gene_id),
      pairsDataBase %>% dplyr::distinct(pairs_id,
                                        .keep_all = TRUE),
      by = "pairs_id"
    )
  # Calculate promoter dominance and other estimates.
  # promoterDominance: fraction of reads within a TSS-3'end site divided by the total number of the associated 3'end in the gene.
  # endDominance: Fraction of reads within a given TSS that expressed a given a 3'end.
  # endFraction: Fraction of 3'end expression over the total number of reads of the gene
  # startFraction: Fraction of TSS expression compare to total gene expression.
  annotPairsExp <-
    annotPairsExp %>%
    dplyr::group_by(tes_id) %>%
    dplyr::mutate(tes_cpm = sum(cpm)) %>% # note this
    group_by(promoter_id) %>%
    dplyr::mutate(tss_cpm = sum(cpm)) %>%
    dplyr::group_by(pairs_id) %>%
    dplyr::mutate(pairs_cpm = sum(cpm)) %>%
    group_by(gene_id) %>%
    mutate(gene_cpm = sum(cpm)) %>%
    dplyr::mutate(
      promoterDominance = pairs_cpm / tes_cpm,
      endDominance = pairs_cpm / tss_cpm,
      endFraction = tes_cpm / gene_cpm,
      startFraction = tss_cpm / gene_cpm
    )
  annotPairsExp <- annotPairsExp %>%
    dplyr::select(-c(cpm)) %>% dplyr::rename(pairs_read_counts = read_counts) %>%
    mutate(pairType=ifelse(is.na(gene_id), "novelPair", "known"))
  return(annotPairsExp)
}


#' estimateTranscriptionalBias
#' transcriptional bias estimates as in (Alfonso-Gonzalez, C, et al., 2022)
#' @param promoterDominanceTable Table of promoter dominance summaring transcript pairs expression
#' @param method option "chisq" or "fisher" test
#'
#' @return table of significantly biased genes
#' @export
#' @import GenomicFeatures GenomicAlignments S4Vectors dplyr
#' @examples
estimateTranscriptionalBias  <-
  function(promoterDominanceTable, method) {
    apa.genes <-
      promoterDominanceTable %>%
      distinct(tes_id, .keep_all = TRUE) %>%
      group_by(gene_id) %>%
      filter(n() > 1) %>%
      pull(gene_id)
    atss.genes <-
      promoterDominanceTable %>%
      distinct(promoter_id, .keep_all = TRUE) %>%
      group_by(gene_id) %>%
      filter(n() > 1) %>%
      pull(gene_id)
    atss.apa.genes <- intersect(apa.genes, atss.genes)
    headData <-
      promoterDominanceTable %>%
      filter(gene_id %in% atss.apa.genes) %>%
      dplyr::select(gene_id, tes_id, pairs_cpm, promoter_id)
    perGeneList <- split(headData, f = headData$gene_id)
    couplingsmatrix <- lapply(perGeneList, function(x) {
     # Produce combination matrix of every TES and Promoter of the gene
      x1 <- x %>%
        maditr::dcast(tes_id ~ promoter_id, value.var = "pairs_cpm") %>%
        as.data.frame(.)
      rownames(x1) <- x1$tes_id
      x1$tes_id <- NULL
      return(x1)
    })
    couplingsmatrix <-
      lapply(couplingsmatrix, function(x) {
        x[is.na(x)] <- 0
        x
      })
    # Add pseudocount
    couplingsmatrix <- lapply(couplingsmatrix, function(x) {
      x + 0.7
    })
    # For every per gene matrix perform chisq test and monte carlo simulation to obtain better p-value estimates.
    if(method=="chisq"){
      res <- list()
      # make chisq test
      perGeneChisqTest <- lapply(couplingsmatrix, function(x) {
        chisqRes <- list()
        x2 <- chisq.test(x, simulate.p.value = TRUE)
        result <- data.frame(p.value = x2$p.value, x_squared = x2$statistic)
        return(result)
      })
      # summarize genes
      res$affectedGenes <-
        as.data.frame(do.call(rbind, perGeneChisqTest)) %>%
        dplyr::rename(p.value.chisq = 1) %>%
        dplyr::mutate(gene_id = rownames(.))
      res$affectedGenes$p.adj.chisq <-
        p.adjust(res$affectedGenes$p.value.chisq , method = "BH")
      res$affectedGenes <- res$affectedGenes %>%
        dplyr::select(c(gene_id,
                        x_squared,
                        p.value.chisq,
                        p.adj.chisq)) %>%
        group_by(gene_id) %>%
        arrange(gene_id)
      #sqrt <-
      #  as.data.frame(do.call(rbind,
      #                        perGeneChisqTest$sum_sqrt)) %>%
       # dplyr::rename(x_squared = 1) %>%
       # dplyr::mutate(gene_id = rownames(.))
      #res$affectedGenes <- left_join(res$affectedGenes,
       #                              sqrt,
        #                             by = "gene_id")
      # get residuals per pair
      residuals <- lapply(couplingsmatrix, function(x) {
        td <- chisq.test(x, simulate.p.value = TRUE)
        #x1 <- x2$residuals
        tdo <- as.data.frame(td$observed) %>%
          dplyr::mutate(new_junID=rownames(.))
        tdo <- reshape2::melt( tdo) %>%
          mutate(pairs_id = paste0(new_junID,":",variable)) %>%
          dplyr::select(pairs_id, value) %>%
          dplyr::rename(observed=value)
        tde <- as.data.frame(td$expected) %>%
          dplyr::mutate(new_junID=rownames(.))
        tde <- reshape2::melt( tde) %>%
          mutate(pairs_id = paste0(new_junID,":",variable)) %>%
          dplyr::select(pairs_id, value) %>%
          dplyr::rename(expected=value)
        tdr <- as.data.frame(td$residuals) %>%
          dplyr::mutate(new_junID=rownames(.))
        tdr <- reshape2::melt( tdr) %>%
          dplyr::mutate(pairs_id = paste0(new_junID,":",variable),
                        gene_id=gsub("\\:.*","",new_junID) ) %>%
          dplyr::select(pairs_id, value, gene_id) %>%
          dplyr::rename(residuals=value)
        d<- left_join(tdo, tde, by="pairs_id")
        d<- left_join(d, tdr, by="pairs_id")
        return(d)
      }
      )
      res$stats <- residuals %>% do.call(rbind,.) %>%
        group_by(gene_id) %>%
        arrange(gene_id)
      return(res)
    }else{
      # For every per gene matrix perform fisher test and monte carlo simulation to obtain better p-value estimates.
      perGeneFisher <- lapply(couplingsmatrix, function(x) {
        res.fish <- fisher.test(x, simulate.p.value = TRUE)
        res.fish <- res.fish$p.value
        return(res.fish)
      })
      resFisher <-
        as.data.frame(do.call(rbind, perGeneFisher)) %>%
        dplyr::rename(p.value.fisher =1) %>%
        dplyr::mutate(gene_id = rownames(.))
      resFisher$p.adj.fisher <-
        p.adjust(resFisher$p.value.fisher , method = "BH")
      return(resFisher)
    }
  }