import pysam as ps
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import sys


# python prep_paraclu.py ${bam_file} results/${bam_file.baseName}


print(sys.argv[1], sys.argv[2])

bam = ps.AlignmentFile(sys.argv[1], "rb")
outpath = sys.argv[2]

tss_count = Counter()
pas_count = Counter()
# TSS is read.ref_start for - strand alignments (+ strand gene)
read_counts = 0
for read in bam.fetch():
    read_counts += 1
    if not (read.is_secondary or read.is_supplementary or read.is_unmapped or read.mapping_quality == 0):
        chrom = read.reference_name
        tss = read.reference_start + 1  if read.is_reverse else read.reference_end
        pas = read.reference_end if read.is_reverse else read.reference_start + 1

        strand = "+" if read.is_reverse else "-" # gene strand
        
        key_tss = (chrom, strand, tss)
        key_pas = (chrom, strand, pas)
        
        tss_count[key_tss] += 1
        pas_count[key_pas] += 1

print("total reads: " + str(read_counts))
print("total tss passing: " + str(len(tss_count)))
print("total pas passing: " + str(len(pas_count)))

df_tss = pd.DataFrame([
    {"chrom": k[0], "strand": k[1], "start": k[2], "count": v}
    for k, v in tss_count.items()
])

df_pas = pd.DataFrame([
    {"chrom": k[0], "strand": k[1], "end": k[2], "count": v}
    for k, v in pas_count.items()
])


# write output

df_tss.sort_values(by=['chrom','strand','start']).to_csv(
   outpath + "tss_count_input_1based.tsv",
   sep="\t",index=False,header=False)

df_pas.sort_values(by=['chrom','strand','end']).to_csv(
   outpath + "pas_count_input_1based.tsv",
   sep="\t",index=False,header=False)