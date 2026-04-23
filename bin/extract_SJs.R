# This script extracts splice junctions from a GTF annotation (preferred)
# If no GTF, it extracts them from the BAM (required) using intronProspector

# Usage:
# From GTF: Rscript extract_SJs.R --gtf annotation.gtf --out annotation_SJ.out.tab
# From BAM: Rscript extract_SJs.R --bam reads.bam --genome genome.fa --out reads_SJ.out.tab
#           

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  i <- which(args == flag)
  if (length(i) == 0) return(NULL)
  args[i + 1]
}

gtf_file    <- get_arg("--gtf")
bam_file    <- get_arg("--bam")
genome_file <- get_arg("--genome")
sj_out_file <- get_arg("--out")

ip_bin <- "intronProspector"
if (is.null(sj_out_file)) stop("--out is required")
if (is.null(gtf_file) && is.null(bam_file)) stop("Provide --gtf or --bam")
if (!is.null(bam_file) && is.null(genome_file)) stop("--genome is required when using --bam")

library(dplyr)
library(tidyr)
library(readr)

strand_map <- c("+" = 1, "-" = 2, "." = 0)

if (!is.null(gtf_file)) {

  # --- Extract SJs from GTF ---
  message("Extracting splice junctions from GTF: ", gtf_file)

  gtf <- rtracklayer::import(gtf_file) %>% data.frame()

  exons <- gtf %>%
    filter(type == "exon") %>%
    group_by(transcript_id) %>%
    arrange(start, .by_group = TRUE)

  introns <- exons %>%
    summarise(
      seqnames = first(seqnames),
      strand = first(strand),
      intron_starts = list(end[-length(end)] + 1),
      intron_ends = list(start[-1] - 1),
      .groups = "drop"
    ) %>%
    unnest(cols = c(intron_starts, intron_ends))

  sj_out <- introns %>%
    mutate(
      strand_code = strand_map[strand],
      motif = 0,
      annotated = 1,
      unique_reads = 0,
      multi_reads = 0,
      max_overhang = 0
    ) %>%
    select(seqnames, intron_starts, intron_ends, strand_code,
           motif, annotated, unique_reads, multi_reads, max_overhang)

} else {

  # --- Extract SJs from BAM using intronProspector ---
  # https://github.com/diekhans/intron-prospector
  # BED6 output columns: chrom, start (0-based), end (exclusive), name, score, strand
  message("Extracting splice junctions from BAM: ", bam_file)

  tmp_bed <- tempfile(fileext = ".bed")
  ret <- system2(ip_bin,
                 args = c(
                   "--genome-fasta", genome_file,
                   paste0("--intron-bed6=", tmp_bed),
                   "-C", "0.0","-r","10",
                   bam_file
                 ),
                 stdout = TRUE, stderr = TRUE)

  if (!file.exists(tmp_bed) || file.info(tmp_bed)$size == 0) {
    stop("intronProspector failed or produced no output:\n", paste(ret, collapse = "\n"))
  }

  ip <- read_tsv(tmp_bed,
                 col_names = c("chrom", "start", "end", "name", "score", "strand"),
                 col_types = "ciiici",
                 comment = "#")

  sj_out <- ip %>%
    mutate(
      intron_starts = start + 1L,  # BED 0-based -> STAR 1-based
      intron_ends   = end,          # BED exclusive end = STAR 1-based inclusive
      strand_code   = strand_map[strand],
      motif         = 0L,
      annotated     = 0L,
      unique_reads  = score,
      multi_reads   = 0L,
      max_overhang  = 0L
    ) %>%
    select(chrom, intron_starts, intron_ends, strand_code,
           motif, annotated, unique_reads, multi_reads, max_overhang)

}

write_tsv(sj_out, sj_out_file, col_names = FALSE)
cat("SJ.out.tab created:", sj_out_file, "\n")
