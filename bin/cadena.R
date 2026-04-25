script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE))))
devtools::load_all(file.path(script_dir, "..", "lib", "cadena"))
library(GenomicRanges)
library(Rsamtools)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)

# 
# bampath
# custom gene intervals from peak-to-peak map, after merging with reference to rescue intergenic peaks
# star splice junctions
# custom links
# strandedness(default no change, alternative reverse the bam)
# chr(optional), do one chr at a time


bampath = args[1]
custom_intervals_path = args[2]
sj_path = args[3]
strandedness = args[4] # "default" or "reverse"
chr = args[5]
outpath = args[6]
libtype = if (length(args) >= 7) args[7] else "10x"
chr_list = if (length(args) >= 8) strsplit(args[8], " ")[[1]] else chr

# Step 1: import custom gene intervals
message("importing peak to peak intervals...")

custom_tx = rtracklayer::import.gff(custom_intervals_path)
tss_pas_links = prepareLinksDatabase(custom_tx,1,1) # default tss.window and tes.window are now set to 0


# infer the gene interval based on earliest start and latest end of the txs
gene_ranges <- custom_tx %>%
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
custom_genes <- makeGRangesFromDataFrame(
  gene_ranges,
  keep.extra.columns = TRUE
)

# Step 2: import splice junctions and add gene id and juncID

message("importing splice junctions and adding gene ids...")
sjs = data.table::fread(sj_path,col.names = c(
    "chromosome",
    "start",
    "end",
    "strand",
    "intron_motif",
    "annotated",
    "unique_reads",
    "multi_mapping_reads",
    "maximum_overhang"
    )) %>% mutate(strand = recode(strand, `0` = "*", `1` = "+",`2` = "-")) %>%
        makeGRangesFromDataFrame()

# add gene id to splice junctions
hits <- findOverlaps(sjs, custom_genes)
sjs$gene_id = NA_character_
sjs$gene_id[queryHits(hits)] = custom_genes$gene_id[subjectHits(hits)]
# add juncID following LATER naming syntax
sjs = sjs %>%  mutate(juncID = paste0(seqnames,":",start,"-",end))

# write custom interval transcripts, peaks, and labels:
custom_links_to_ids = tss_pas_links$pairDataBase %>% dplyr::select(transcript_id, pairs_id) %>% 
    left_join(custom_tx %>% data.frame() %>% select(transcript_id, tss_peak_name, tss_valid, pas_peak_name, pas_valid))
data.table::fwrite(custom_links_to_ids, paste0(outpath,"pro_tss.ids.tsv"),sep="\t")

for (current_chr in chr_list) {
  current_outpath = if (length(chr_list) > 1) paste0(outpath, current_chr, "_") else outpath

  # Step 3: annotate reads with junctions and tss-pas
  #################### Main counting function ####################

  message("annotating reads to junctions and links for chr ", current_chr)

  ###### main action:
  reads_annotated = read_to_junctions_custom(bampath, sjs, tss_pas_links, fixstrand=strandedness, chr=current_chr)

  if (nrow(reads_annotated) == 0) {
    message("No reads annotated for chr ", current_chr, " — skipping")
    next
  }

  reads_annotated <- reads_annotated %>%
    arrange(gene_id, new_junID) %>%
    group_by(gene_id) %>%
    mutate(J_id = if_else(
      is.na(new_junID),
      NA_character_,
      paste0("J", as.integer(if_else(strand == "+", dense_rank(new_junID), dense_rank(desc(new_junID)))))
    )) %>%
    ungroup()
    
  reads_annotated = reads_annotated %>% rename("junction" = "new_junID")
  # new_juncID is the corrected reference junction
  sjs2 = reads_annotated %>% filter(!is.na(junction)) %>% dplyr::select(gene_id, strand, junction, J_id) %>% distinct()
  data.table::fwrite(sjs2 %>% data.frame(), paste0(current_outpath,"sjs.ids.tsv"),sep="\t")

  # condense reads_annotated into the junction chain:
  reads_annotated_condensed = reads_annotated %>%
    group_by(read_id, gene_id, tes_id, promoter_id) %>%
    summarise(
      ejc = (function(j) {
        j <- j[!is.na(j)]
        if (length(j) == 0) return(NA_character_)
        paste(j[order(as.integer(sub("^.", "", j)))], collapse = ",")
      })(unique(J_id)),
      .groups = "drop"
    ) %>% mutate(tes_id = gsub(".*:","",tes_id),
                promoter_id = gsub(".*:","",promoter_id)) %>% distinct()

  outfile = paste0(current_outpath, current_chr, "_reads.tsv.gz")

  message("writing output to ", outfile)
  data.table::fwrite(reads_annotated_condensed, outfile, sep="\t", compress = "gzip", na = "NA")
}


