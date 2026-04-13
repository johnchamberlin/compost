import pandas as pd
import numpy as np
import sys
import warnings

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

warnings.simplefilter(action='ignore', category=FutureWarning)

from cadena.utils import read_connections


def add_gene_to_peaks(contacts, pann):
    df = contacts.merge(
        pann[["peak_name", "gene_id"]],
        left_on="tss_peak_name",
        right_on="peak_name",
        how="left"
    ).rename(columns={"gene_id": "tss_gene_id"}).merge(
        pann[["peak_name", "gene_id"]],
        left_on="pas_peak_name",
        right_on="peak_name",
        how="left"
    ).rename(columns={"gene_id": "pas_gene_id"})[
        ["count", "tss_peak_name", "pas_peak_name", "median_aln_length",
         "tss_read_count", "pas_read_count", "span_frac", "tss_gene_id", "pas_gene_id"]
    ]
    return df


def logistic_decreasing(x, L, k, x0, c):
    return c + L / (1 + np.exp(k * (x - x0)))


def x_at_y(y_target, L, k, x0, c):
    if not (c < y_target < c + L):
        return np.nan
    return x0 + (1 / k) * np.log(L / (y_target - c) - 1)


connections = sys.argv[1]
annotations = sys.argv[2]
outpdf      = sys.argv[3]

spans   = read_connections(connections, min_reads=50, max_ratio=100)
peakann = pd.read_csv(annotations, sep="\t")

with PdfPages(outpdf) as pdf:

    # --- Plot 1: Soft clip lengths by peak type ---
    fig, ax = plt.subplots(figsize=(8, 5))
    d = peakann.copy()
    d["above500"] = d["Score"] >= 500
    ylim = 100
    sns.violinplot(
        x="peak_type",
        y="median_clip",
        hue="above500",
        data=d[d["median_clip"] < ylim],
        palette="Set2",
        split=True,
        inner="quart",
        ax=ax
    )
    ax.set_title("Soft clip lengths by peak type")
    ax.set_xlabel("Peak type")
    ax.set_ylabel("Median clip length")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # --- Plot 2: Spanning fraction vs transcript length ---
    spans2 = add_gene_to_peaks(spans, peakann)
    tdf = spans2[spans2["tss_gene_id"] == spans2["pas_gene_id"]].copy()
    df = tdf

    x = df["median_aln_length"].to_numpy()
    y = df["span_frac"].to_numpy()

    num_bins = 20
    quantile = 0.95

    bins = np.linspace(x.min(), x.max(), num_bins + 1)
    df["bin"] = np.digitize(x, bins)
    frontier = df.groupby("bin").agg({
        "median_aln_length": "mean",
        "span_frac": lambda g: np.quantile(g, quantile)
    }).dropna()
    frontier_x = frontier["median_aln_length"].to_numpy()
    frontier_y = frontier["span_frac"].to_numpy()

    p0 = [1.0, 0.01, np.median(frontier_x), 0.0]
    try:
        params, _ = curve_fit(
            logistic_decreasing, frontier_x, frontier_y,
            p0=p0,
            bounds=([0, 0, frontier_x.min(), 0], [1, np.inf, frontier_x.max(), 1])
        )
        L, k, x0, c = params
        x50 = x_at_y(0.5, L, k, x0, c)
        fit_succeeded = True
    except RuntimeError:
        fit_succeeded = False

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(x, y, alpha=0.15, s=2, label="Peak pairs > 50 reads")
    ax.scatter(frontier_x, frontier_y, color="red", s=4, label=f"{int(quantile*100)}th percentile frontier")

    if fit_succeeded:
        x_fit = np.linspace(x.min(), x.max(), 500)
        ax.plot(x_fit, logistic_decreasing(x_fit, *params), color="black", linewidth=2, label="Sigmoid fit")
        if not np.isnan(x50):
            ax.axvline(x=x50, ymin=0, ymax=0.5, color='gray', linewidth=2, alpha=0.8)
            ax.axhline(y=0.5, xmin=0, xmax=(x50 - x.min()) / (x.max() - x.min()), color='gray', linewidth=2, alpha=0.8)
            ax.text(x50 + 10, 0.51, f"{int(round(x50))} bp", fontsize=12, ha='left', va='bottom')

    ax.set_title("Spanning fraction vs transcript length")
    ax.set_xlabel("Transcript length (bp)")
    ax.set_ylabel("Spanning fraction")
    ax.legend(fontsize=8)
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # --- Plot 3: Median alignment length distribution ---
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(spans["median_aln_length"], bins=50, color="steelblue", edgecolor="white")
    ax.set_title("Median alignment length distribution")
    ax.set_xlabel("Median alignment length (bp)")
    ax.set_ylabel("Count")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # Plot 4: Median alignment length per read
    # --- Plot 3: Median alignment length distribution (weighted by read count) ---
    fig, ax = plt.subplots(figsize=(8, 5))
    lengths = np.repeat(spans["median_aln_length"].to_numpy(), spans["count"].to_numpy().astype(int))
    ax.hist(lengths, bins=50, color="skyblue", edgecolor="white")
    ax.set_title("Median alignment length distribution (weighted by read count)")
    ax.set_xlabel("Median alignment length (bp)")
    ax.set_ylabel("Read count")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)