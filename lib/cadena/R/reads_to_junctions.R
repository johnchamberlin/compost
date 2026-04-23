#' Title
#'
#' @param bamJunctions: bam junctions obtained from read_to_junctions functions
#' @param shortReadJunctions: read junctions from STAR
#' @param dist.5ss : distance to 5ss splice site to compare to reference to be corrected in the long reads
#' @param dist.3ss : distance to 3ss splice site to compare to reference to be corrected in the long reads
#' @importFrom rlang .data
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @return
#' @export
#'
#' @examples
read_refjunction_correction <- function(bamJunctions,
                                        shortReadJunctions,
                                        dist.5ss,
                                        dist.3ss) {
  # Take bam file and summarize all junctions
  tmpRef <- unique(unlist(bamJunctions))
  if (length(tmpRef) == 0) {
    message("No junctions found in reads — returning empty result")
    return(data.frame(raw_juncID = character(0), new_junID = character(0)))
  }
  tmpRef$raw_juncID <-
    paste0(seqnames(tmpRef), ":",
           start(tmpRef), "-",
           end(tmpRef))
  refJunctions <- tmpRef
  # Assign found junctions in bam to reference junctions
  hitsovlps <- findOverlaps(refJunctions, shortReadJunctions)
  if (length(hitsovlps) == 0) {
    message("No junction overlaps found — returning empty result")
    return(data.frame(raw_juncID = character(0), new_junID = character(0)))
  }
  junctionsInShort <- refJunctions[queryHits(hitsovlps), ]
  message("Assigned to junctions reads: " , length(unique(junctionsInShort)))
  junctionsInShort$start_short <-
    start(shortReadJunctions[subjectHits(hitsovlps),])
  junctionsInShort$end_short <-
    end(shortReadJunctions[subjectHits(hitsovlps),])
  junctionsInShort$strand_short <-
    strand(shortReadJunctions[subjectHits(hitsovlps),])
  junctionsInShort$dist_start <-
    abs(start(shortReadJunctions[subjectHits(hitsovlps),]) -  start(junctionsInShort))
  junctionsInShort$dist_end <-
    abs(end(shortReadJunctions[subjectHits(hitsovlps),]) -  end(junctionsInShort))
  junctionsShortAssignments <-
    junctionsInShort %>% as.data.frame(., row.names = NULL) %>%
    arrange(dist_end + dist_start) %>%
    dplyr::mutate (totalDist = dist_end + dist_start) %>%
    dplyr::distinct(raw_juncID, .keep_all = TRUE)
  # filter by distance and assign new id
  message("Junctions before filtering " ,
      nrow(junctionsShortAssignments))
  junctionsShortAssignments <-
    junctionsShortAssignments %>%
    dplyr::filter(dist_start <= dist.5ss &
                    dist_end <= dist.3ss) %>%
    dplyr::mutate(new_junID = paste0(seqnames, ":",
                                     start_short, "-",
                                     end_short))
  message("Junctions after filtering " ,
      nrow(junctionsShortAssignments))
  return(junctionsShortAssignments)
}

#' merge junctions to reads
#'
#' @param bamJunctions
#'
#' @return
#' @export
#'
#' @examples
make_junction_database <- function(bamJunctions) {
  junBam <- unlist(bamJunctions)
  junBam$raw_juncID <-
    paste0(seqnames(junBam),":",
           start(junBam),"-",
           end(junBam))
  junBam$read_id <- names(unlist(bamJunctions))
  # add later
  #junBam[junBam$read_id %in% fullLength,]
  junBam <-
    junBam %>% as.data.frame(., row.names = NULL) %>% dplyr::select(raw_juncID, read_id)
  return(junBam)
}


#' custom tss-pas junctions to reads assignments
#'
#' @param bamPath: path to minimap2 bam file
#' @param refJunAnnot: reference_junctions: set of junctions of interest classified using SaiLoR
#' @param linksDatabase: path to tss-pas links from prepareLinksDatabase
#'
#' @return
#' @export
#'
#' @examples
read_to_junctions_custom <- function(bamPath, refJunAnnot, linksDatabase, fixstrand = "default", chr = "all") {

  if(chr == "all") {
    # Extract all tags from BAM (this keeps CL, CB, NM, AS, etc.)
    message("Extracting per-read tags from all chrs (including CL and CB)")
    bam_tags <- Rsamtools::scanBam(
      bamPath,
      param = Rsamtools::ScanBamParam(tag = c("CL", "CB"), what = "qname")
    )[[1]]
    
    # NEW: Safely extract both CL and CB tags (some reads may not have one of them)
    tag_df <- tibble(
      read_id = bam_tags$qname,
      CL = if (!is.null(bam_tags$tag$CL)) bam_tags$tag$CL else NA,
      CB = if (!is.null(bam_tags$tag$CB)) bam_tags$tag$CB else NA
    )
  }
  if(chr != "all") {
    message(paste("Extracting per-read tags from chr:", chr, "(including CL and CB)"))

    # Get chromosome lengths (to define range)
    hdr <- Rsamtools::scanBamHeader(bamPath)
    chr_len <- hdr[[1]]$targets[chr]

    # Define which region to read (entire chromosome)
    param <- Rsamtools::ScanBamParam(
      which = GRanges(chr, IRanges(1, chr_len)),
      what = "qname",
      tag = c("CL", "CB")
    )

    # Read BAM for that chromosome only
    bam_tags <- Rsamtools::scanBam(bamPath, param = param)[[1]]

    # Build a tidy tibble of read IDs and tags
    tag_df <- tibble(
      read_id = bam_tags$qname,
      CL = if (!is.null(bam_tags$tag$CL)) bam_tags$tag$CL else NA,
      CB = if (!is.null(bam_tags$tag$CB)) bam_tags$tag$CB else NA
    )
  }

  q = 10  # MAPQ filter
 message("Loading minimap2 alignments (primary only, MAPQ filtered)")
  flag_filter <- Rsamtools::scanBamFlag(
    isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = FALSE,
    isDuplicate = FALSE
  )

  if(chr == "all") {
    param <- Rsamtools::ScanBamParam(
      flag = flag_filter,
      what = c("qname", "flag", "mapq"),
      mapqFilter = q
    )
  } else {
    hdr <- Rsamtools::scanBamHeader(bamPath)
    chr_len <- hdr[[1]]$targets[chr]
    param <- Rsamtools::ScanBamParam(
      which = GRanges(chr, IRanges(1, chr_len)),
      flag = flag_filter,
      what = c("qname", "flag", "mapq"),
      mapqFilter = q
    )
  }

  bamAlignments <- GenomicAlignments::readGAlignments(bamPath, param = param, use.names = TRUE)

  
  if (fixstrand == "default") {
    message("Reversing strand information")
    current_strand <- strand(bamAlignments)
    flipped_strand <- ifelse(current_strand == "+", "-",
                             ifelse(current_strand == "-", "+", "*"))
    strand(bamAlignments) <- flipped_strand
  }

  # Filter reads by full length
  message("Filtering only full length reads")
  fullLengths <- get_sailor_custom(linksDatabase, bamAlignments, tss.ntwindow = 50, tes.ntwindow = 150)
  bamAlignments <- bamAlignments[names(bamAlignments) %in% fullLengths$pairsTested$name, ]

  message("Extract junctions per reads")
  read_junctions <- GenomicAlignments::junctions(bamAlignments, use.mcols = TRUE)

  message("Assign junctions to reference junctions")
  reference_junction_data <- read_refjunction_correction(read_junctions, refJunAnnot, 10, 10)

  message("Create database per read")
  read.features <- fullLengths$pairsTested %>%
    group_by(name) %>%
    distinct(name, .keep_all = TRUE) %>%
    dplyr::rename(read_id = name)

  if (nrow(reference_junction_data) > 0) {
    bam_recovered_junctions <- make_junction_database(read_junctions)
    junction_assignments <- inner_join(bam_recovered_junctions, reference_junction_data, by = "raw_juncID")
    reads_to_junctions <- right_join(junction_assignments, read.features, by = "read_id")
  } else {
    reads_to_junctions <- read.features %>% dplyr::mutate(new_junID = NA_character_, strand = NA_character_)
  }

  # --- NEW SECTION: add start and end of the read alignment ---
  message("Adding alignment start and end coordinates")
  read_positions <- tibble(
    read_id = names(bamAlignments),
    read_start = start(bamAlignments),
    read_end = end(bamAlignments)
  )
  reads_to_junctions <- left_join(reads_to_junctions, read_positions, by = "read_id")

  # Join CL and CB tags
  reads_to_junctions <- left_join(reads_to_junctions, tag_df, by = "read_id")

  message("Per read database data (with CL and CB tags)")
  return(reads_to_junctions)
}
