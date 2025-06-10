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

### Ignore warnings
import warnings 
warnings.filterwarnings("ignore")

from functions_script import *

import argparse

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

def main():
    args = parse_args()
    neighbourhood_size = args.neighbourhood_size
    fimo_filepath = args.fimo_filepath
    output_dir = args.output_dir
    genome_fasta = args.genome_fasta
    ignore_repeats = args.ignore_repeats
    chip_filepath = args.chip_filepath
    boundary = args.boundary

    #create output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    report_path = os.path.join(output_dir, 'summary.txt')
    report = open(report_path, 'w')


    #print the fimo file head
    fimo = pd.read_csv(fimo_filepath, sep='\t')
    fimo.head()

    #get FASTA file
    fimo_to_neighbourhood(fimo_filepath, output_dir, genome_fasta, neighbourhood_size=neighbourhood_size)

    #plot heatmap
    plot_fasta_heatmap(output_dir + 'fimo_neighbourhood_50.fasta')

    plot_logo_from_fasta(output_dir + 'fimo_neighbourhood_50.fasta')

    #median value of score in fimo_neighbourhood_50.fasta
    neigh_bed = pd.read_csv(output_dir + 'fimo_neighbourhood_50.bed', sep='\t', header=None)

    median_all = neigh_bed[4].median()
    report.write(f"Median score (all neighbourhoods): {median_all}\n")

    chip_filtered_neighbourhood(output_dir + 'fimo_neighbourhood_50.bed', chip_filepath, genome_fasta)
    plot_fasta_heatmap(output_dir + 'fimo_neighbourhood_50_filtered.fasta')
    plot_logo_from_fasta(output_dir + 'fimo_neighbourhood_50_filtered.fasta')

    #filter the neighbourhoods that intersect the boundary
    intersect_filtered_neighbourhood_boundary(output_dir + 'fimo_neighbourhood_50_filtered.bed', boundary, genome_fasta)

    #plot the heatmap
    print('Plotting heatmap for filtered neighbourhoods that intersect the boundary')
    plot_fasta_heatmap(output_dir + 'fimo_neighbourhood_50_filtered_int_boundary.fasta')
    plot_logo_from_fasta(output_dir + 'fimo_neighbourhood_50_filtered_int_boundary.fasta')

    print('Plotting heatmap for filtered neighbourhoods that do not intersect the boundary')
    plot_fasta_heatmap(output_dir + 'fimo_neighbourhood_50_filtered_no_int_boundary.fasta')
    plot_logo_from_fasta(output_dir + 'fimo_neighbourhood_50_filtered_no_int_boundary.fasta')

    #read in int bed file
    int_bed = pd.read_csv(output_dir + 'fimo_neighbourhood_50_filtered_int_boundary.bed', sep='\t', header=None)

    median_int = int_bed[4].median()
    report.write(f"Median score (intersecting boundary): {median_int}\n")

    #read in no int bed file
    no_int_bed = pd.read_csv(output_dir + 'fimo_neighbourhood_50_filtered_no_int_boundary.bed', sep='\t', header=None)

    median_no_int = no_int_bed[4].median()
    report.write(f"Median score (non-intersecting boundary): {median_no_int}\n")

    #plot the distribution of the scores
    plt.figure()
    sns.distplot(int_bed[4], label='Intersect')
    sns.distplot(no_int_bed[4], label='No Intersect')
    plt.xlabel('Scores')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'score_distribution.png'))
    plt.close()

    #plot using histplot
    plt.figure()
    sns.histplot(int_bed[4], label='Intersect', kde=True, stat='density')
    sns.histplot(no_int_bed[4], label='No Intersect', kde=True, stat='density')
    plt.xlabel('Scores')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'score_histogram.png'))
    plt.close()

    statistic, pvalue = stats.ks_2samp(int_bed[4], no_int_bed[4])
    report.write(f"KS test statistic: {statistic}\n")
    report.write(f"KS test p-value: {pvalue}\n")

    #inputs for matrix generation - fimo_neighbourhood_50.bed, chip_filepath, boundary, output_dir
    matrix_neighbourhood(output_dir + 'fimo_neighbourhood_50.bed', chip_filepath, boundary, output_dir)

    #count number of 1s in each column
    matrix = pd.read_csv(output_dir + 'neighbourhood_matrix.tsv', sep='\t') 

    matrix.head()

    #unique values in the chip column
    print('Unique values in the chip column')
    print(matrix['chip'].unique())

    #unique values in the boundary column
    print('Unique values in the boundary column')
    print(matrix['boundary'].unique())


    #count number of 1s in each column
    #number of 1s in matrix['chip'}
    print('Number of 1s in the chip column')
    print(matrix['chip'].sum())

    #number of 1s in matrix['boundary']
    print('Number of 1s in the boundary column')
    print(matrix['boundary'].sum())

    report.close()

if __name__ == "__main__":
    main()