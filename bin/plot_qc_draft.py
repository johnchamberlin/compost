import pysam as ps
import pandas as pd
import pyranges as pr
import numpy as np
from intervaltree import IntervalTree
from scipy.sparse import dok_matrix
import sys
import seaborn as sns

import matplotlib.pyplot as plt


import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)



from cadena.utils import read_connections


# this plots the spanning read fraction for each sample
# i.e., for a given pair of peaks, we have the total reads (tss + pas - spanning)
# and we have the spanning reads, so we can compute the fraction of spanning reads for each sample and plot it
# we also ahve the average length of the spanning alignments

def add_gene_to_peaks(contacts, pann):
    df = contacts.merge(
    pann[["peak_name","gene_id"]],
    left_on="tss_peak_name",
    right_on="peak_name",
    how="left"
    ).rename(columns={"gene_id":"tss_gene_id"}).merge(
        pann[["peak_name","gene_id"]],
        left_on="pas_peak_name",
        right_on="peak_name",
        how="left"
    ).rename(columns={"gene_id":"pas_gene_id"})[["count","tss_peak_name","pas_peak_name","median_aln_length",
                                                "tss_read_count","pas_read_count","span_frac","tss_gene_id","pas_gene_id"]]
    return(df)


connections = sys.argv[1] # path to the links tsv
annotations = sys.argv[2] # path to the annotation tsv
outpdf = sys.argv[3] # path to output pdf

# tss:pas ratio omit if >100 either direction, or span count < 50
spans = read_connections(connections, min_reads=50, max_ratio=100)
peakann = pd.read_csv(annotations, sep="\t")


# plot 1: soft clip lengths
d = peakann
# Plot into the correct subplot

d["above500"] = d["Score"]>=500

ylim=100 # this needs to adjust for the species somehow
sns.violinplot(
    x="peak_type",
    y="median_clip",
    hue="above500",
    data=d[d.median_clip < ylim],
    palette="Set2",
    split=True,
    inner="quart"
)


# plot 2: spanning fraction
spans2 = add_gene_to_peaks(spans, peakann)

def logistic_decreasing(x, L, k, x0, c):
    return c + L / (1 + np.exp(k * (x - x0)))

def x_at_y(y_target, L, k, x0, c):
    if not (c < y_target < c + L):
        return np.nan
    return x0 + (1/k) * np.log(L / (y_target - c) - 1)

num_bins = 20
quantile = 0.95
# Filter
tdf = spans2[spans2["tss_gene_id"] == spans2["pas_gene_id"]]

# filter by genes with only one tss and one pas?

#df = tdf[tdf.groupby('tss_gene_id')['tss_gene_id'].transform('count') == 1].copy()
#df = df[df.groupby('pas_gene_id')['pas_gene_id'].transform('count') == 1].copy()
df = tdf
x = df["median_aln_length"].to_numpy()
y = df["span_frac"].to_numpy()

# Bin and compute quantile frontier
bins = np.linspace(x.min(), x.max(), num_bins + 1)
df["bin"] = np.digitize(x, bins)
frontier = df.groupby("bin").agg({
    "median_aln_length": "mean",
    "span_frac": lambda g: np.quantile(g, quantile)
}).dropna()
frontier_x = frontier["median_aln_length"].to_numpy()
frontier_y = frontier["span_frac"].to_numpy()

# Fit
p0 = [1.0, 0.01, np.median(frontier_x), 0.0]
params, _ = curve_fit(
    logistic_decreasing, frontier_x, frontier_y,
    p0=p0,
    bounds=([0, 0, frontier_x.min(), 0], [1, np.inf, frontier_x.max(), 1])
)
L, k, x0, c = params
x50 = x_at_y(0.5, L, k, x0, c)

# Plot
fig, ax = plt.subplots()
ax.scatter(x, y, alpha=0.85, s=2, label="Peaks > 50 reads")
ax.scatter(frontier_x, frontier_y, color="red", s=4, label=f"{int(quantile*100)}th percentile")

x_fit = np.linspace(x.min(), x.max(), 500)
ax.plot(x_fit, logistic_decreasing(x_fit, *params), color="black", linewidth=2, label="Sigmoid fit")

ax.axvline(x=x50, ymin=0, ymax=0.5, color='gray', linewidth=2, alpha=0.8)
ax.axhline(y=0.5, xmin=0, xmax=(x50 - x.min()) / (x.max() - x.min()), color='gray', linewidth=2, alpha=0.8)
ax.text(x50 + 10, 0.51, f"{int(round(x50))} bp", fontsize=12, ha='left', va='bottom')

ax.set_xlabel("Transcript length")
ax.set_ylabel("Spanning fraction")
ax.legend(fontsize=8)
plt.tight_layout()
plt.show()

# plot 3: median alignment length distribution 
