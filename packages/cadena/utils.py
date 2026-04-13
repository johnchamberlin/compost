"""
utils.py
---------
Plotting and data manipulation functions for long-read scRNA-seq
2025 Sept 8
John Chamberlin
"""

import pysam as ps
from collections import defaultdict
import polars as pl
import matplotlib.pyplot as plt
from pybedtools import BedTool
import pyranges as pr
import pandas as pd
import numpy as np

# plotting functions

def plot_connections(bampath,gtfpath,peakspath,celltype,geneid=None,coords=None):

    """
    Plot the destination of long reads according to TSS and the origin
    according to PAS based on pre-defined peaks, a gene annotation, and
    optional cell type labels

    Pass gene id or coords in format coords=[chrom,start,end,strand]
    """
    
    fig, ax = plt.subplots(1,2,figsize=(12,5)) 
    
    bamf = ps.AlignmentFile(bampath,"rb")
    gtf =  pd.read_csv(gtfpath, sep="\t", comment="#", header=None,
                  names=["Chromosome", "Source", "Feature", "Start", "End",
                         "Score", "Strand", "Frame", "Attributes"])
    
    gtf['gene_id'] = gtf['Attributes'].str.extract(r'gene_id "([^"]+)"')
    genes_df = gtf[gtf.Feature == "gene"].copy()
    genes = pr.PyRanges(genes_df[["Chromosome", "Start", "End", "Strand", "Attributes"]])
    peaks = pr.read_bed(peakspath)
    overlaps = peaks.join(genes, strandedness="same").df
    if geneid is not None:
        gene_loc = genes_df[genes_df['gene_id']==geneid]
        
        gene_chrom = gene_loc['Chromosome'].iloc[0]
        gene_start = gene_loc['Start'].iloc[0]-1
        gene_end = gene_loc['End'].iloc[0]+1
        gene_strand = gene_loc['Strand'].iloc[0]
    elif coords is not None:
        gene_chrom=coords[0]
        gene_start=coords[1]
        gene_end=coords[2]
        gene_strand = coords[3]
    else:
        return(None)

    tss_coords = []
    pas_coords = []
    aln_lengths = []
    
    for read in bamf.fetch(gene_chrom,gene_start, gene_end):
        # skip multimappred reads
        if celltype != "all":
            try:
                if read.get_tag("CL") != celltype:
                    continue
            except KeyError:
                continue
            
        if(gene_strand == "+"):
            if read.is_reverse and not (read.is_secondary or read.is_supplementary or read.mapping_quality == 0):
                tss_coords.append(read.reference_start)
                pas_coords.append(read.reference_end)
                aln_lengths.append(read.query_alignment_length)
        else: # gene is minus strand
            # only keep positive strand reads
            if not (read.is_reverse or read.is_secondary or read.is_supplementary or read.mapping_quality == 0):
                tss_coords.append(read.reference_end)
                pas_coords.append(read.reference_start)
                aln_lengths.append(read.query_alignment_length)
    # overlap 'TSS' (alignment boundary) with peaks, to associate corresponding'PAS' with TSS peak    
    tss_df = pd.DataFrame({
    "Chromosome": gene_chrom,
    "Start": tss_coords,
    "End":None,
    "Strand": gene_strand,
    "PAS": pas_coords,
    "aln_length": aln_lengths
    })
    tss_df['End']=tss_df['Start']+1
    tss_df = pl.from_pandas(pr.PyRanges(tss_df).join(peaks,strandedness="same").df)
    # overlap 'PAS' (alignment boundary) with peaks, to associate corresponding 'TSS' with PAS peak    

    pas_df = pd.DataFrame({
    "Chromosome": gene_chrom,
    "Start": pas_coords,
    "End":None,
    "Strand": gene_strand,
    "TSS": tss_coords,
    "aln_length": aln_lengths
    })
    
    pas_df['End']=pas_df['Start']+1
    pas_df = pl.from_pandas(pr.PyRanges(pas_df).join(peaks,strandedness="same").df)

    if(gene_strand == "+"):
        tss_df = tss_df.filter(pl.col("PAS") <= gene_end+1000)
        tss_df = tss_df.filter(pl.col("Start") >= gene_start-1000)
        pas_df = pas_df.filter(pl.col("TSS") >= gene_start-1000)
        pas_df = pas_df.filter(pl.col("Start") <= gene_end+1000)
    else:
        tss_df = tss_df.filter(pl.col("PAS") >= gene_start-1000)
        tss_df = tss_df.filter(pl.col("Start") <= gene_end+1000)
        pas_df = pas_df.filter(pl.col("TSS") <= gene_end+1000)
        pas_df = pas_df.filter(pl.col("Start") >= gene_start-1000)


    # now filter to keep the more abundant PAS and TSS peaks:
    major_tss = tss_df.select(['Name', 'Score']).unique()
    major_tss = major_tss.filter((pl.col('Name').str.starts_with('tss_')) & (pl.col('Score') > 500))
    
    # keep tss with 5% of total
    tot_tss = major_tss.select(pl.col("Score").sum()).item()
    major_tss = major_tss.filter(
        pl.col('Score') > 0.05 * tot_tss
    )
    major_tss = major_tss.sort("Score",descending=True).get_column('Name').to_list()
    

    major_pas = pas_df.select(['Name', 'Score']).unique()
    major_pas = major_pas.filter((pl.col('Name').str.starts_with('pas_')) & (pl.col('Score') > 500))
    
    # keep tss with 5% of total
    tot_pas = major_pas.select(pl.col("Score").sum()).item()
    major_pas = major_pas.filter(
        pl.col('Score') > 0.05 * tot_pas
    )
    major_pas = major_pas.sort("Score",descending=True).get_column('Name').to_list()
    # p lotting part
    # need to filter by the gene intervals to excluve PAS that are in timbuktu
    for peak in major_tss:
        subset = tss_df.filter(pl.col("Name") == peak)
        peak_start = subset.get_column('Start_b').unique().item()
        peak_end = subset.get_column('End_b').unique().item()
        values = subset["PAS"].to_list()
        bins = np.linspace(gene_start-1000, gene_end+1000, 100 + 1)
        n, bins, patches = ax[0].hist(values, bins=bins, alpha=0.5, label=peak,edgecolor='gray')
        color = patches[0].get_facecolor()
        ax[0].axvline( (peak_start+peak_end) / 2, linewidth = 1, linestyle='--', color=color) 

    ax[0].set_xlabel("Coordinate")
    ax[0].set_ylabel("Frequency")
    ax[0].set_xlim(gene_start-1000,gene_end+1000)
    ax[0].legend()
    bamf.close()

    for peak in major_pas:
        subset = pas_df.filter(pl.col("Name") == peak)
        peak_start = subset.get_column('Start_b').unique().item()
        peak_end = subset.get_column('End_b').unique().item()
        values = subset["TSS"].to_list()
        bins = np.linspace(gene_start-1000, gene_end + 1000, 100 + 1)
        n, bins, patches = ax[1].hist(values, bins=bins, alpha=0.5, label=peak,edgecolor='gray')
        color = patches[0].get_facecolor()
        ax[1].axvline( (peak_start+peak_end) / 2, linewidth = 1, linestyle='--', color=color) 

    ax[1].set_xlabel("Coordinate")
    ax[1].set_ylabel("Frequency")
    ax[1].set_xlim(gene_start-1000,gene_end+1000)
    ax[1].legend(loc="best")

    fig.suptitle(f"{geneid} PAS and TSS connections in {celltype}",fontsize=14)
    plt.tight_layout()
    tss_out = tss_df.to_pandas()
    pas_out = pas_df.to_pandas()

    return fig, ax, tss_out, pas_out



# data manipulation functions
# tsv path,
# min reads to keep a connection, # max ratio between tss and pas reads to keep a connection

def read_connections(tsv_path, min_reads = 50, max_ratio = 100):
    """
    Read peak to peak read count data produced by my custom pipeline
    
    """
    df = pd.read_csv(tsv_path, sep="\t")
    df["total_reads"] = df["tss_read_count"] + df["pas_read_count"] - df["count"]
    df["span_frac"] = df["count"] / df["total_reads"]
    dff = df[df["count"] > min_reads].copy()
    #dff["length_bin"] = pd.qcut(dff["median_aln_length"], q=10)
    dfff = dff[(dff["tss_read_count"] / dff["pas_read_count"] < max_ratio) & (dff["pas_read_count"] / dff["tss_read_count"] < max_ratio)].copy()
    dfff["length_bin"] = pd.qcut(dfff["median_aln_length"], q=10)

    return(dfff)


def annotate_peaks(peaks,
    peak_links,
    peak_annotation,
    gtfpath,
    homopolymerpath,libtype="ont",min_reads=50):
    """
    Annotate the peak to peak read count data with homopolymer presence (A7),
    soft clipping stats (5' and 3'), and estimate peak validity (as TSS or PAS vs artifact)
    """
    homopol = pr.read_bed(homopolymerpath)
    peaks = peaks.count_overlaps(homopol,strandedness="same")
    links = peak_links.copy()
    peakann = peak_annotation.copy()
    
    gtf = pr.read_gtf(gtfpath)
    genes = gtf[gtf.Feature == "gene"].copy()
    exons = gtf[gtf.Feature == "exon"].copy()
    
    peaks_to_exons = peaks.nearest(exons, strandedness="same")
    peaks_to_exons = peaks_to_exons[peaks_to_exons.Distance == 0]
    peaks_to_genes = peaks[~peaks.Name.isin(peaks_to_exons.Name)].nearest(genes,strandedness="same")
    peaks_to_genes = peaks_to_genes[peaks_to_genes.Distance == 0]
    
    peaks_to_features = pr.concat([peaks_to_exons,peaks_to_genes])
    peaks_to_intergenic = peaks[~peaks.Name.isin(peaks_to_features.Name)].df
    peaks_to_intergenic["Feature"] = "intergenic"
    peaks_to_intergenic["gene_id"] = None
    
    peaksGenes = pd.concat([peaks_to_features[["Chromosome","Start","End","Name","Score","Strand","Feature","gene_id","NumberOverlaps"]].df,peaks_to_intergenic])
    
    peakann = peakann.merge(peaksGenes,left_on="peak_name",right_on="Name").drop(["Name","total"],axis=1).copy()

    # filtering changes depending on the data type
    if(libtype == "ont"):
        tss_minclip = 12
        tss_maxclip = 45
        pas_minclip = 15
        pas_maxclip = 30
    elif(libtype == "argentag"):
        tss_minclip = 12
        tss_maxclip = 45
        pas_minclip = 105
        pas_maxclip = 120
    else:
        tss_minclip = 12
        tss_maxclip = 45
        pas_minclip = 23
        pas_maxclip = 37
    
    peakann["valid"] = np.where(
    (
        (peakann.peak_type == "tss")
        & (peakann.median_clip > tss_minclip)
        & (peakann.median_clip < tss_maxclip)
        & (peakann.frac_GGGG > 0.5)
        & (peakann.Score > min_reads)
    )
    | (
        (peakann.peak_type == "pas")
        & (peakann.median_clip > pas_minclip)
        & (peakann.median_clip < pas_maxclip)
        & (peakann.Score > min_reads)
        & (peakann.NumberOverlaps == 0)
    ),
    True,
    False
    )
    peakann["width"] = peakann["End"]-peakann["Start"]
    return(peakann)


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
    

def make_peaks_gtf(peak_ann, links_path, outpath):
    links = read_connections(links_path, min_reads = 10)
    linkGene = add_gene_to_peaks(links,peak_ann)
    linkGtf = linkGene[["tss_peak_name","pas_peak_name","tss_gene_id","pas_gene_id"]].copy()
    linkGtf["tss_gene_id"] = linkGtf["tss_gene_id"].combine_first(linkGtf["pas_gene_id"])
    linkGtf["pas_gene_id"] = linkGtf["pas_gene_id"].combine_first(linkGtf["tss_gene_id"])
    linkGtf = linkGtf[linkGtf["tss_gene_id"] == linkGtf["pas_gene_id"]].copy()

    tss = peak_ann[peak_ann["peak_type"] == "tss"].copy()[["Chromosome","Start","End","Strand","peak_name","valid"]]
    tss.columns = ["chromosome","tss_start","tss_end","strand","tss_peak_name","tss_valid"]
    pas = peak_ann[peak_ann["peak_type"] == "pas"].copy()[["Chromosome","Start","End","Strand","peak_name","valid"]]
    pas.columns = ["chromosome","pas_start","pas_end","strand","pas_peak_name","pas_valid"]

    linkGtf = linkGtf.merge(tss).merge(pas)
    # fix logic to keep the midpoint as the coordinate
    linkGtf["pas_mid"] = (linkGtf['tss_end'] + linkGtf['tss_start']) // 2
    linkGtf["tss_mid"] = (linkGtf['pas_end'] + linkGtf['pas_start']) // 2
    linkGtf["start"] = linkGtf[['tss_mid','pas_mid']].min(axis=1)
    linkGtf["end"] = linkGtf[['tss_mid','pas_mid']].max(axis=1)

    # Recover intergenic (NA gene) peaks:
    # if TSS is genic and PAS is NA, assign gene id from TSS, and vice versa
    linkGtf["tss_gene_id"] = linkGtf["tss_gene_id"].fillna(linkGtf["pas_gene_id"])
    linkGtf["pas_gene_id"] = linkGtf["pas_gene_id"].fillna(linkGtf["tss_gene_id"])

    # For intergenic peaks where both TSS and PAS are NA
    mask = linkGtf["tss_gene_id"].isna() & linkGtf["pas_gene_id"].isna()
    if mask.any():
        intergenic_ids = ["intergenic_" + str(i+1) for i in range(mask.sum())]
        linkGtf.loc[mask, "tss_gene_id"] = intergenic_ids
        linkGtf.loc[mask, "pas_gene_id"] = intergenic_ids

    linkGtf["gene_id"] = linkGtf["tss_gene_id"]

    # for discordant peak intervals, make a concatenated id
    linkGtf.loc[linkGtf["tss_gene_id"] != linkGtf["pas_gene_id"], "gene_id"] = linkGtf["tss_gene_id"] + "_" + linkGtf["pas_gene_id"]
 

    linksGtfFmt = linkGtf[["chromosome","start","end","strand","tss_peak_name","pas_peak_name","tss_valid","pas_valid","tss_gene_id"]].copy()
    linksGtfFmt["feature"] = "transcript"
    linksGtfFmt = linksGtfFmt.rename(columns={'tss_gene_id':'gene_id'})
    df = linksGtfFmt.copy()
    df['n'] = df.groupby('gene_id').cumcount() + 1
    
    # Create transcript_id
    df['transcript_id'] = df['gene_id'] + "_custom_" + df['n'].astype(str)
    df['seqname'] = df['chromosome']
    df['source'] = 'custom'
    df['score'] = '.'
    df['frame'] = '.'
    
    # Define column order for GTF (standard 9)
    core_cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']
    
    # All other columns become GTF attributes
    attr_cols = [c for c in df.columns if c not in core_cols + ['chromosome']]
    
    def make_attributes(row):
        # Format: key "value"; key "value";
        return '; '.join([f'{col} "{row[col]}"' for col in attr_cols]) + ';'
    
    df['attribute'] = df.apply(make_attributes, axis=1)
    
    # Keep only the 9 standard columns
    gtf = df[core_cols + ['attribute']]
    gtf.to_csv(outpath, sep='\t', header=False, index=False, quoting=3)



### drafts

def annotate_peaks_ori(peakpath,peakconnectionspath,peakannotationpath,gtfpath,homopolymerpath,libtype="ont",min_reads=50):
    """
    Annotate the peak to peak read count data with homopolymer presence (A7),
    soft clipping stats (5' and 3'), and estimate peak validity (as TSS or PAS vs artifact)
    """
    homopol = pr.read_bed(homopolymerpath)
    peaks = pr.read_bed(peakpath)
    peaks = peaks.count_overlaps(homopol,strandedness="same")
    links = read_connections(peakconnectionspath,min_reads = min_reads)
    peakann = pd.read_csv(peakannotationpath,sep="\t")
    
    gtf = pr.read_gtf(gtfpath)
    genes = gtf[gtf.Feature == "gene"].copy()
    exons = gtf[gtf.Feature == "exon"].copy()
    
    peaks_to_exons = peaks.nearest(exons, strandedness="same")
    peaks_to_exons = peaks_to_exons[peaks_to_exons.Distance == 0]
    peaks_to_genes = peaks[~peaks.Name.isin(peaks_to_exons.Name)].nearest(genes,strandedness="same")
    peaks_to_genes = peaks_to_genes[peaks_to_genes.Distance == 0]
    
    peaks_to_features = pr.concat([peaks_to_exons,peaks_to_genes])
    peaks_to_intergenic = peaks[~peaks.Name.isin(peaks_to_features.Name)].df
    peaks_to_intergenic["Feature"] = "intergenic"
    peaks_to_intergenic["gene_id"] = None
    
    peaksGenes = pd.concat([peaks_to_features[["Chromosome","Start","End","Name","Score","Strand","Feature","gene_id","NumberOverlaps"]].df,peaks_to_intergenic])
    
    peakann = peakann.merge(peaksGenes,left_on="peak_name",right_on="Name").drop(["Name","total"],axis=1).copy()

    # filtering changes depending on the data type
    if(libtype == "ont"):
        tss_minclip = 12
        tss_maxclip = 45
        pas_minclip = 15
        pas_maxclip = 30
    elif(libtype == "argentag"):
        tss_minclip = 12
        tss_maxclip = 45
        pas_minclip = 105
        pas_maxclip = 120
    else:
        tss_minclip = 12
        tss_maxclip = 45
        pas_minclip = 23
        pas_maxclip = 37
    
    peakann["valid"] = np.where(
    (
        (peakann.peak_type == "tss")
        & (peakann.median_clip > tss_minclip)
        & (peakann.median_clip < tss_maxclip)
        & (peakann.frac_GGGG > 0.5)
        & (peakann.Score > min_reads)
    )
    | (
        (peakann.peak_type == "pas")
        & (peakann.median_clip > pas_minclip)
        & (peakann.median_clip < pas_maxclip)
        & (peakann.Score > min_reads)
        & (peakann.NumberOverlaps == 0)
    ),
    True,
    False
    )
    peakann["width"] = peakann["End"]-peakann["Start"]
    return(peakann)