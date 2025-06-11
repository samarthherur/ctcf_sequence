# === Import required libraries ===
# Numerical arrays, dataframes, plotting, file operations, and bioinformatics tools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import subprocess 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import scipy.stats as stats
import motif_strength_script
from itertools import combinations

# === Suppress non-critical warnings for cleaner output ===
### Ignore warnings
import warnings 
warnings.filterwarnings("ignore")

# === Import custom pipeline functions from functions_script ===
from functions_script import *

import argparse

# === Parse command-line arguments for input/output paths and parameters ===
def parse_args():
    parser = argparse.ArgumentParser(description="Process FIMO neighbourhood analysis.")
    parser.add_argument("--neighbourhood_size", type=int, default=50, help="Size of the neighbourhood window (bp).")
    parser.add_argument("--fimo_filepath", type=str, required=True, help="Path to the FIMO TSV file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory for all outputs.")
    parser.add_argument("--genome_fasta", type=str, required=True, help="Path to the genome FASTA file.")
    parser.add_argument("--ignore_repeats", action="store_true", help="Skip repeat-masked regions.")
    parser.add_argument("--chip_filepath", type=str, required=True, help="Path to the ChIP-seq peaks BED file.")
    parser.add_argument("--boundary", type=str, required=True, help="Path to the boundary BED file.")
    return parser.parse_args()

# === Main pipeline execution ===
def main():
    args = parse_args()
    neighbourhood_size = args.neighbourhood_size
    fimo_filepath = args.fimo_filepath
    output_dir = args.output_dir
    genome_fasta = args.genome_fasta
    ignore_repeats = args.ignore_repeats
    chip_filepath = args.chip_filepath
    boundary = args.boundary

    # -- Create output directory if it does not exist --
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # -- Open summary report file for logging results --
    report_path = os.path.join(output_dir, 'summary.txt')
    report = open(report_path, 'w')

    # -- Generate neighbourhood FASTA and BED from FIMO hits --
    fimo_to_neighbourhood(fimo_filepath, output_dir, genome_fasta, neighbourhood_size=neighbourhood_size)

    # -- Plot heatmap and sequence logo for raw neighbourhood sequences --
    plot_fasta_heatmap(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}.fasta"))

    plot_logo_from_fasta(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}.fasta"))

    # -- Compute and log median score of all neighbourhood windows --
    neigh_bed = pd.read_csv(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}.bed"), sep='\t', header=None)

    median_all = neigh_bed[4].median()
    report.write(f"Median score (all neighbourhoods): {median_all}\n")

    # -- Filter neighbourhood windows by ChIP-seq peaks and re-plot --
    chip_filtered_neighbourhood(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}.bed"), chip_filepath, genome_fasta)
    plot_fasta_heatmap(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}_filtered.fasta"))
    plot_logo_from_fasta(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}_filtered.fasta"))

    # -- Split filtered windows by boundary intersection and plot subsets --
    intersect_filtered_neighbourhood_boundary(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}_filtered.bed"), boundary, genome_fasta)

    #plot the heatmap
    plot_fasta_heatmap(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}_filtered_int_boundary.fasta"))
    plot_logo_from_fasta(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}_filtered_int_boundary.fasta"))

    plot_fasta_heatmap(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}_filtered_no_int_boundary.fasta"))
    plot_logo_from_fasta(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}_filtered_no_int_boundary.fasta"))

    #read in int bed file
    int_bed = pd.read_csv(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}_filtered_int_boundary.bed"), sep='\t', header=None)

    median_int = int_bed[4].median()
    report.write(f"Median score (intersecting boundary): {median_int}\n")

    #read in no int bed file
    no_int_bed = pd.read_csv(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}_filtered_no_int_boundary.bed"), sep='\t', header=None)

    median_no_int = no_int_bed[4].median()
    report.write(f"Median score (non-intersecting boundary): {median_no_int}\n")

    # determine common range and bin edges
    min_score = min(int_bed[4].min(), no_int_bed[4].min())
    max_score = max(int_bed[4].max(), no_int_bed[4].max())
    bins = np.linspace(min_score, max_score, 30)   # you can adjust 30 to any number of bins you like

    plt.figure(figsize=(8, 6))

    # histogram outlines
    sns.histplot(int_bed[4], bins=bins, stat="density", element="step", fill=False,
                label="Intersect", linewidth=1.5)
    sns.histplot(no_int_bed[4], bins=bins, stat="density", element="step", fill=False,
                label="No Intersect", linewidth=1.5)

    # KDE overlays
    sns.kdeplot(int_bed[4], bw_adjust=1, linewidth=2, label="Intersect KDE")
    sns.kdeplot(no_int_bed[4], bw_adjust=1, linewidth=2, label="No Intersect KDE")

    plt.xlabel("Scores")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "score_distribution.png"))
    plt.close()

    #plot using histplot
    plt.figure()
    sns.histplot(int_bed[4], label='Intersect', kde=True, stat='density')
    sns.histplot(no_int_bed[4], label='No Intersect', kde=True, stat='density')
    plt.xlabel('Scores')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'score_histogram.png'))
    plt.close()

    # -- Perform Kolmogorov–Smirnov test between intersecting and non-intersecting scores --
    statistic, pvalue = stats.ks_2samp(int_bed[4], no_int_bed[4])
    report.write(f"KS test statistic: {statistic}\n")
    report.write(f"KS test p-value: {pvalue}\n")

    # -- Run motif_strength_script to classify motifs by boundary and strand relationships --
    filtered_fimo = os.path.join(output_dir, 'fimo_filtered.tsv')
    motif_strength_script.filter_fimo_by_ctcf_peaks(
        fimo_filepath,
        chip_filepath,
        filtered_fimo
    )
    classified_motifs = os.path.join(output_dir, 'classified_motifs.tsv')
    motif_strength_script.main(
        boundary,
        filtered_fimo,
        classified_motifs
    )
    report.write(f"Filtered FIMO saved to: {filtered_fimo}\n")
    report.write(f"Classified motifs saved to: {classified_motifs}\n")

    # -- Load classified motifs and log counts, medians, and KS tests per class --
    classified_motifs_filepath = classified_motifs
    classified_df = pd.read_csv(classified_motifs_filepath, sep='\t')
    total = len(classified_df)
    report.write(f"Total classified motifs: {total}\n")

    # Counts per class
    classes = np.sort(classified_df['classification'].unique())
    for cls in classes:
        count = (classified_df['classification'] == cls).sum()
        report.write(f"Count {cls}: {count}\n")

    # Median score per class
    for cls in classes:
        median_cls = classified_df.loc[classified_df['classification']==cls, 'score'].median()
        report.write(f"Median score {cls}: {median_cls}\n")

    # KS tests between classes
    for c1, c2 in combinations(classes, 2):
        data1 = classified_df.loc[classified_df['classification']==c1, 'score']
        data2 = classified_df.loc[classified_df['classification']==c2, 'score']
        stat, pval = stats.ks_2samp(data1, data2)
        report.write(f"KS {c1} vs {c2}: stat={stat}, p={pval}\n")

    # -- Plot and save density distributions for each motif class --
    plt.figure()
    for cls in classes:
        sns.kdeplot(classified_df.loc[classified_df['classification']==cls, 'score'], shade=True, label=cls)
    plt.xlabel('Scores')
    plt.ylabel('Density')
    plt.legend()
    density_fig = os.path.join(output_dir, 'classification_density.png')
    plt.savefig(density_fig)
    plt.close()
    report.write(f"Density plot saved to: {density_fig}\n")

    # -- Filter out sequences containing Ns or lowercase bases and log base counts --
    filtered_fasta = os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}_filtered.fasta")
    seqs_wr, n_wr = fasta_without_repeats(filtered_fasta)
    report.write(f"Sequences without repeats: {n_wr}\n")

    # Count nucleotides in repeat-free sequences
    count_dict = {}
    for rec in seqs_wr:
        for nt in rec.seq:
            count_dict[nt] = count_dict.get(nt, 0) + 1
    report.write("Base counts without repeats:\n")
    for nt in sorted(count_dict):
        report.write(f"{nt}: {count_dict[nt]}\n")

    # -- Write cleaned FASTA of repeat-free sequences --
    wr_fasta = filtered_fasta.replace('.fasta', '_wr.fasta')
    from Bio import SeqIO
    SeqIO.write(seqs_wr, wr_fasta, 'fasta')
    report.write(f"Saved repeat-free FASTA to: {wr_fasta}\n")

    # -- Visualize repeat-free neighbourhood sequences --
    plot_fasta_heatmap(wr_fasta)
    plot_logo_from_fasta(wr_fasta)
    report.write("Saved heatmap and logo for repeat-free sequences\n")

    # -- Generate binary overlap matrix for ChIP and boundary on each window --
    matrix_neighbourhood(os.path.join(output_dir, f"fimo_neighbourhood_{neighbourhood_size}.bed"), chip_filepath, boundary, output_dir)

    #count number of 1s in each column
    matrix = pd.read_csv(os.path.join(output_dir, f"neighbourhood_matrix.tsv"), sep='\t') 

    # -- Close the summary report file --
    report.close()

if __name__ == "__main__":
    main()