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
from functions_script import plot_intersect_distribution, count_neighbourhood_windows
from motif_strength_script import summarize_classified_motifs

import argparse

# === Parse command-line arguments for input/output paths and parameters ===
def parse_args():
    parser = argparse.ArgumentParser(description="Process FIMO neighbourhood analysis.")
    parser.add_argument("--n_a", type=int, default=50, help="Upstream extension size (bp).")
    parser.add_argument("--n_b", type=int, default=50, help="Downstream extension size (bp).")
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
    n_a = args.n_a
    n_b = args.n_b
    fimo_filepath = args.fimo_filepath
    output_dir = args.output_dir
    genome_fasta = args.genome_fasta
    ignore_repeats = args.ignore_repeats
    chip_filepath = args.chip_filepath
    boundary = args.boundary

    # -- Create output directory if it does not exist --
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # Redirect all stdout and stderr to the summary report file
    report_path = os.path.join(output_dir, 'summary.txt')
    report_file = open(report_path, 'w')
    sys.stdout = report_file
    sys.stderr = report_file

    # -- Generate neighbourhood FASTA and BED from FIMO hits --
    fimo_to_neighbourhood(fimo_filepath, output_dir, genome_fasta, n_a, n_b)

    # -- Plot heatmap and sequence logo for raw neighbourhood sequences --
    plot_fasta_heatmap(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}.fasta"))

    plot_logo_from_fasta(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}.fasta"))

    # -- Compute and log median score of all neighbourhood windows --
    neigh_bed = pd.read_csv(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}.bed"), sep='\t', header=None)

    median_all = neigh_bed[4].median()
    print(f"Median score (all neighbourhoods): {median_all}")

    # -- Filter neighbourhood windows by ChIP-seq peaks and re-plot --
    chip_filtered_neighbourhood(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}.bed"), chip_filepath, genome_fasta)
    plot_fasta_heatmap(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered.fasta"))
    plot_logo_from_fasta(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered.fasta"))

    # -- Split filtered windows by boundary intersection and plot subsets --
    intersect_filtered_neighbourhood_boundary(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered.bed"), boundary, genome_fasta)

    #plot the heatmap
    plot_fasta_heatmap(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered_int_boundary.fasta"))
    plot_logo_from_fasta(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered_int_boundary.fasta"))

    plot_fasta_heatmap(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered_no_int_boundary.fasta"))
    plot_logo_from_fasta(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered_no_int_boundary.fasta"))

    #read in int bed file
    int_bed = pd.read_csv(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered_int_boundary.bed"), sep='\t', header=None)

    median_int = int_bed[4].median()
    print(f"Median score (intersecting boundary): {median_int}")

    #read in no int bed file
    no_int_bed = pd.read_csv(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered_no_int_boundary.bed"), sep='\t', header=None)

    median_no_int = no_int_bed[4].median()
    print(f"Median score (non-intersecting boundary): {median_no_int}")

    plot_intersect_distribution(int_bed[4], no_int_bed[4], output_dir, n_a, n_b)

    # -- Perform Kolmogorov–Smirnov test between intersecting and non-intersecting scores --
    statistic, pvalue = stats.ks_2samp(int_bed[4], no_int_bed[4])
    print(f"KS test statistic: {statistic}")
    print(f"KS test p-value: {pvalue}")

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
    print(f"Filtered FIMO saved to: {filtered_fimo}")
    print(f"Classified motifs saved to: {classified_motifs}")

    summarize_classified_motifs(classified_motifs, output_dir)

    # -- Filter out sequences containing Ns or lowercase bases and log base counts --
    filtered_fasta = os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered.fasta")
    seqs_wr, n_wr = fasta_without_repeats(filtered_fasta)
    print(f"Sequences without repeats: {n_wr}")

    # Count nucleotides in repeat-free sequences
    count_dict = {}
    for rec in seqs_wr:
        for nt in rec.seq:
            count_dict[nt] = count_dict.get(nt, 0) + 1
    print("Base counts without repeats:")
    for nt in sorted(count_dict):
        print(f"{nt}: {count_dict[nt]}")

    # -- Write cleaned FASTA of repeat-free sequences --
    wr_fasta = filtered_fasta.replace('.fasta', '_wr.fasta')
    SeqIO.write(seqs_wr, wr_fasta, 'fasta')
    print(f"Saved repeat-free FASTA to: {wr_fasta}")

    # -- Visualize repeat-free neighbourhood sequences --
    plot_fasta_heatmap(wr_fasta)
    plot_logo_from_fasta(wr_fasta)
    print("Saved heatmap and logo for repeat-free sequences")

    count_neighbourhood_windows(output_dir, n_a, n_b)

    # -- Generate binary overlap matrix for ChIP and boundary on each window --
    matrix_neighbourhood(os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}.bed"), chip_filepath, boundary, output_dir)

    #count number of 1s in each column
    matrix = pd.read_csv(os.path.join(output_dir, f"neighbourhood_matrix.tsv"), sep='\t') 

    # Close the summary report file
    report_file.close()

if __name__ == "__main__":
    main()