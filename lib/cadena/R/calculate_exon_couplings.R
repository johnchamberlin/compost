#' Calculate exon couplings
#'
#' @param junctions: junctions obtained from Long read data assigned to read_ids
#' @param reference_junctions: set of junctions of interest classified using SaiLoR
#'
#' @return couplings result list file with couplings for TSS and TES
#' @export
#'
#' @examples
calculate_exon_couplings <- function(junctions,reference_junctions){
  couplingsResult <- list()
  res_couplings_tss <- compute_exon_tss_couplings(junctions,reference_junctions)
  res_couplings_tes <- compute_exon_3end_couplings(junctions,reference_junctions)
  couplingsResult$TSS <- res_couplings_tss
  couplingsResult$TES <- res_couplings_tes
  couplingsResult <- unlist(couplingsResult, recursive = FALSE)
  return(couplingsResult)
}


#' Calculate exon couplings
#' # I changed this function uses monte carlo simulation based chisq test to compute couplings
#' @param junctions: junctions obtained from Long read data assigned to read_ids
#' @param reference_junctions: set of junctions of interest classified using SaiLoR
#' @param B: number of replicates for monte carlo simulation in chisq test
#' @return couplings result list file with couplings for TSS and TES
#' @export
#'
#' @examples
calculate_exon_couplings_custom <- function(junctions,reference_junctions,B=2000){
  couplingsResult <- list()
  res_couplings_tss <- compute_exon_tss_couplings_custom(junctions,reference_junctions,B=B)[[1]]
  couplingsmatrix_tss <- compute_exon_tss_couplings_custom(junctions,reference_junctions,B=B)[[2]]
  res_couplings_tes <- compute_exon_3end_couplings_custom(junctions,reference_junctions,B=B)[[1]]
  couplingsmatrix_tes <- compute_exon_3end_couplings_custom(junctions,reference_junctions,B=B)[[2]]
  couplingsResult$TSS <- res_couplings_tss
  couplingsResult$TES <- res_couplings_tes
  couplingsResult <- unlist(couplingsResult, recursive = FALSE)
  return(list(couplingsResult,list(TSS_matrix=couplingsmatrix_tss,TES_matrix=couplingsmatrix_tes)))
}
