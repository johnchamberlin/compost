#!/usr/bin/env python

import pysam
import pandas as pd
import pyranges as pr
import numpy as np
from intervaltree import IntervalTree
from scipy.sparse import dok_matrix
import sys
from collections import defaultdict
from statistics import median
from Bio.Seq import Seq 

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from cadena_py.utils import make_peaks_gtf
from cadena_py.utils import annotate_peaks

bedpath=sys.argv[1] # TSS and PAS peaks bed
genepath=sys.argv[2] # original gene annotation
bampath=sys.argv[3] # BAM file of long reads
prefix=sys.argv[4] # output path prefix
homopol=sys.argv[5] # homopolymer annotation file path
trans_spliced=sys.argv[6].lower() == "true" if len(sys.argv) > 6 else False
libtype=sys.argv[7] if len(sys.argv) > 7 else "10x"

def get_softclip_lengths(read):
    if not read.cigartuples:
        return 0, 0
    left_len = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
    right_len = read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0
    return left_len, right_len


def get_align_and_softclip(read): # get the 8bp at softclip:read boundary at the "TSS"
    '''
    For analyzing the 5' cap signature, we want to look at the softclip:aligned boundary 
    sequence at the TSS end of the read. If it contains 4G, evidence of 5p cap
    If only 3G, lack of 5p cap. Other sequence, other weirdness.
    Assumes reads align opposite to gene strand.
    '''

    readseq = read.query_sequence
    if read.is_reverse:
        # we want the left clip and the first four bp of the readseq
        tss_clip_len = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
        boundary_seq = readseq[tss_clip_len-4:tss_clip_len+4]
    else:
        tss_clip_len = read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0
        boundary_seq = readseq[-(tss_clip_len+4):-(tss_clip_len-4)]
        # if read is forward strand, we need to reverse complement
        boundary_seq = str(Seq(boundary_seq).reverse_complement())
    return(boundary_seq)

def count_motif(boundary_seq): # analyze the 5p cap signature
    BSEQ = boundary_seq.upper()
    if 'GGGG' in BSEQ:
        return 'GGGG'
    elif 'GGG' in BSEQ:
        return 'GGG'
    else:
        return 'other'

# import the peaks and genes, assign bam file
peaks_pr = pr.read_bed(bedpath)
# keep genes only
genes = pr.read_gtf(genepath)
genes = genes[genes.Feature == "gene"].copy()


bam = pysam.AlignmentFile(bampath, "rb" )

# filter out the antisense artifact peaks
# i.e. remove peaks with antisense partner >100x stronger
peaks = peaks_pr[~peaks_pr.Name.isin(peaks_pr.join(peaks_pr, strandedness="opposite",suffix="_AS").
                               df.query("Score_AS >= 100 * Score").Name)]

peaks = peaks_pr.df

# split tss and pas peaks, reset index so index = row number
tss_peaks = peaks[peaks.Name.str.contains("tss")].reset_index(drop=True)
pas_peaks = peaks[peaks.Name.str.contains("pas")].reset_index(drop=True)

tss_softclip_dists = [defaultdict(int) for _ in range(len(tss_peaks))]
pas_softclip_dists = [defaultdict(int) for _ in range(len(pas_peaks))]

# generate the tree dictionary thingies
tss_trees = {}
pas_trees = {}

for peak in tss_peaks.itertuples():
    key = (peak.Chromosome, peak.Strand)
    tss_trees.setdefault(key, IntervalTree()).addi(
        int(peak.Start),
        int(peak.End),
        peak.Index
    )
    
for peak in pas_peaks.itertuples():
    key = (peak.Chromosome, peak.Strand)
    pas_trees.setdefault(key, IntervalTree()).addi(
        int(peak.Start),
        int(peak.End),
        peak.Index
    )

contact_matrix = dok_matrix((len(tss_peaks),len(pas_peaks)), dtype=np.int64)
# empty vector to hold tss counts:
tss_counts = defaultdict(int)
pas_counts = defaultdict(int)

# keep track of the alignment lengths between peak pairs
contact_length_counts = defaultdict(lambda: defaultdict(int))
# keep track of the 5p cap signatures at TSS peaks
tss_motif_counts = defaultdict(lambda: defaultdict(int))

n = 0
for read in bam.fetch():
    if not(read.is_secondary or read.is_supplementary or read.is_unmapped or read.mapping_quality == 0):
        n+=1
        left_len, right_len = get_softclip_lengths(read)

        # assign read boundaries to pas/tss points
        if read.is_reverse:
            read_strand = "-"
            peak_strand = "+"
            read_pas_pos = read.reference_end - 1
            read_tss_pos = read.reference_start
            read_soft_at_tss = left_len
            read_soft_at_pas = right_len
        else: 
            read_strand = "+"
            peak_strand = "-"
            read_pas_pos = read.reference_start
            read_tss_pos = read.reference_end - 1
            read_soft_at_tss = right_len
            read_soft_at_pas = left_len
        chrom = read.reference_name
        
        # for the peak indices:
        start_pid = None
        end_pid = None

        # if the read is on a chromosome with peaks
        if (chrom, peak_strand) in tss_trees:
            # Find start peak at 5' position
            tss_hits = tss_trees[(chrom, peak_strand)].overlap(read_tss_pos, read_tss_pos + 1)
            if tss_hits:
                start_pid = list(tss_hits)[0].data
                tss_counts[start_pid] += 1
                tss_softclip_dists[start_pid][read_soft_at_tss] += 1
                cap_seq = get_align_and_softclip(read)
                motif = count_motif(cap_seq)
                tss_motif_counts[tss_peaks.Name[start_pid]][motif] += 1
        # Find end peak at 3' position
        if (chrom, peak_strand) in pas_trees:
            pas_hits = pas_trees[(chrom, peak_strand)].overlap(read_pas_pos, read_pas_pos + 1)
            if pas_hits:
                end_pid = list(pas_hits)[0].data
                pas_counts[end_pid] += 1
                pas_softclip_dists[end_pid][read_soft_at_pas] += 1

        if start_pid is not None and end_pid is not None:
            contact_matrix[start_pid, end_pid] += 1
            contact_length_counts[(start_pid, end_pid)][read.query_alignment_length] += 1
            

# convert the sparse matrix to a DataFrame
rows, cols = zip(*contact_matrix.keys())
values = list(contact_matrix.values())

links_df = pd.DataFrame({
    "row_idx": rows,
    "col_idx": cols,
    "count": values,
})

# Add row and column names by mapping indices
links_df["tss_peak_name"] = links_df["row_idx"].map(lambda i: tss_peaks.Name[i])
links_df["pas_peak_name"] = links_df["col_idx"].map(lambda i: pas_peaks.Name[i])

#####
medians = {k: median([val for val, cnt in v.items() for _ in range(cnt)])
           for k, v in contact_length_counts.items()}
links_df["median_aln_length"] = medians.values()
links_df["idx_check"] = medians.keys()

links_df["tss_read_count_fromname"] = links_df["tss_peak_name"].str.split("_").str[1].astype(int)
links_df["pas_read_count_fromname"] = links_df["pas_peak_name"].str.split("_").str[1].astype(int)

# Map read counts
links_df["tss_read_count"] = links_df["row_idx"].map(tss_counts)
links_df["pas_read_count"] = links_df["col_idx"].map(pas_counts)

# links_df is now pairwise connections counts, plus align length info, etc
links_df.to_csv(prefix+"peak_pairwise_connections.tsv", index=False, sep="\t", header=True)
print(f"Total reads checked: {n}")
print(f"Total reads in peaks matrix: {contact_matrix.sum()}")

# tss 5p cap dataframe:

rows = []
for peak_name, counts in tss_motif_counts.items():
    row = {"peak_name": peak_name}
    row.update(counts)
    rows.append(row)

df_motifs = pd.DataFrame(rows).fillna(0)
df_motifs['total'] = df_motifs[['GGG','GGGG','other']].sum(axis=1)
df_motifs['frac_GGG'] = df_motifs['GGG'] / df_motifs['total']
df_motifs['frac_GGGG'] = df_motifs['GGGG'] / df_motifs['total']
df_motifs['frac_other'] = df_motifs['other'] / df_motifs['total']


# tss and pas soft clip lengths
def weighted_median(lengths, counts):
    data = np.repeat(lengths, counts)
    return np.median(data)

tss_medians = {}
for i, dist in enumerate(tss_softclip_dists):
    if dist:  # skip empty peaks
        lengths = np.array(list(dist.keys()))
        counts = np.array(list(dist.values()))
        tss_medians[tss_peaks.Name[i]] = weighted_median(lengths, counts)

# PAS medians
pas_medians = {}
for i, dist in enumerate(pas_softclip_dists):
    if dist:
        lengths = np.array(list(dist.keys()))
        counts = np.array(list(dist.values()))
        pas_medians[pas_peaks.Name[i]] = weighted_median(lengths, counts)

pas_clip_df = pd.DataFrame(list(pas_medians.items()), 
    columns=['peak_name','median_clip']).sort_values("median_clip")
pas_clip_df["peak_type"]="pas"

tss_clip_df = pd.DataFrame(list(tss_medians.items()),
    columns=['peak_name','median_clip']).sort_values("median_clip")
tss_clip_df["peak_type"]="tss"

# combine the median soft clip length data into one df
clip_df = pd.concat([tss_clip_df,pas_clip_df],axis=0,ignore_index=True)

# combine soft clip df wtih the 5p cap signature df, pas peaks will be NaN
annot_df = clip_df.merge(df_motifs,on="peak_name",how="left")




## annotate_peaks(peakpath,peakconnectionspath,peakannotationpath,

peaks_with_genes = annotate_peaks(peaks_pr,
    links_df,
    annot_df,
    genepath,
    homopol,
    libtype=libtype,
    trans_spliced=trans_spliced)

peaks_with_genes.to_csv(prefix+"peak_annotation.tsv", index=False, sep="\t", header=True)
print(f"Output written to {prefix}peak_annotation.tsv")

# create the intervals GTF

make_peaks_gtf(
    peaks_with_genes,
    f"{prefix}peak_pairwise_connections.tsv",
    f"{prefix}peak_intervals.gtf"
)