# this will be adapted to output the full length read coverage rate plot
# ultimately we want to make a panel of plots for all bam files run in one nf

plot_titles = ["Islets_12d_AT","Islets_33_AT","Ret_AT","Ret_ACME_AT","Nvec","Mlei"]

dfs = [islets_H1_12d_AT_annf,islets_33_AT_annf,ret_AT_annf,ret_fx_AT_annf,nvec_annf,mlei_annf]
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

quantile = 0.95
num_bins = 50

def logistic_decreasing(x, L, k, x0, c):
    return c + L / (1 + np.exp(k * (x - x0)))

def x_at_y(y_target, L, k, x0, c):
    # Solve analytically for x where logistic_decreasing(x) = y_target
    if not (c < y_target < c + L):
        return np.nan  # y_target outside curve range
    return x0 + (1/k) * np.log(L / (y_target - c) - 1)

# --- Set up plot grid ---
fig, axes = plt.subplots(2, 3, figsize=(9,6))
axes = axes.flatten()

x50_vals = []

for i, df in enumerate(dfs):

    tdf = df[df["tss_gene_id"] == df["pas_gene_id"]]
    df = tdf[tdf.groupby('tss_gene_id')['tss_gene_id'].transform('count') == 1].copy()
    df2 = df[df.groupby('pas_gene_id')['pas_gene_id'].transform('count') == 1].copy()
    df = df2.copy()
    ax = axes[i]
 
    x = df["median_aln_length"].to_numpy()
    y = df["span_frac"].to_numpy()
    
    # Bin x
    bins = np.linspace(x.min(), x.max(), num_bins + 1)
    df["bin"] = np.digitize(x, bins)
    
    # Compute quantile in each bin
    frontier = df.groupby("bin").agg({
        "median_aln_length": "mean",
        "span_frac": lambda g: np.quantile(g, quantile)
    }).dropna()
    
    frontier_x = frontier["median_aln_length"].to_numpy()
    frontier_y = frontier["span_frac"].to_numpy()
    
    # Fit decreasing logistic
    p0 = [1.0, 0.01, np.median(frontier_x), 0.0]
    params, _ = curve_fit(
        logistic_decreasing,
        frontier_x, frontier_y,
        p0=p0,
        bounds=([0, 0, frontier_x.min(), 0], [1, np.inf, frontier_x.max(), 1])
    )
    
    L, k, x0, c = params
    x50 = x_at_y(0.5, L, k, x0, c)
    x50_vals.append(x50)
    
    # Plot all points and frontier
    ax.scatter(x, y, alpha=0.15, label="Peaks > 50 reads",s=2)
    ax.scatter(frontier_x, frontier_y, color="red", label=f"{int(quantile*100)}th percentile",s=4)
    
    # Sigmoid fit
    x_fit = np.linspace(x.min(), x.max(), 500)
    y_fit = logistic_decreasing(x_fit, *params)
    ax.plot(x_fit, y_fit, color="black", linewidth=2, label="Sigmoid fit")
    
    # --- Dot at y=0.5 ---
    y_dot = 0.5
    ax.text(
    x50+10, y_dot+.01, f"{int(round(x50))} bp",
    fontsize=12,
    ha='left', va='bottom',  # position slightly offset
    )
    # --- Dotted lines from axes ---
    ax.axvline(x=x50, ymin=0, ymax=y_dot, color='gray',linewidth=2,alpha=.8)
    ax.axhline(y=y_dot, xmin=0, xmax=(x50 - x.min())/(x.max()-x.min()), color='gray', linewidth=2,alpha=.8)
    # Labels and legend
    ax.set_title(plot_titles[i])
    ax.set_xlabel("Transcript length")
    ax.set_ylabel("Spanning fraction")
    ax.legend(fontsize=8)
    
plt.subplots_adjust(hspace=0.5)  # increase the vertical space between rows
plt.show()

