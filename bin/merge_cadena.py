import sys
import glob
import pandas as pd

bam_id = sys.argv[1]

reads_files    = sorted(glob.glob("*_reads.tsv.gz"))
interval_files = sorted(glob.glob("*_pro_tss.ids.tsv"))
sj_files       = sorted(glob.glob("*_sjs.ids.tsv"))

reads = pd.concat([pd.read_csv(f, compression="gzip", sep="\t") for f in reads_files], ignore_index=True)
reads["pairs_id"] = reads[["gene_id", "promoter_id", "tes_id"]].astype(str).agg(":".join, axis=1)
reads.to_csv(f"{bam_id}_cadena_reads.tsv.gz", index=False, sep="\t", compression="gzip")

intervals = pd.concat([pd.read_csv(f, sep="\t") for f in interval_files], ignore_index=True)
intervals.to_csv(f"{bam_id}_cadena_pro_tss.ids.tsv", index=False, sep="\t")

sjs = pd.concat([pd.read_csv(f, sep="\t") for f in sj_files], ignore_index=True)
sjs.to_csv(f"{bam_id}_cadena_sjs.ids.tsv", index=False, sep="\t")
