import os
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import argparse
import json
import subprocess

from nplb_helper_functions import *

def parse_args():
    parser = argparse.ArgumentParser(description="Compute avg gene distances per cluster")
    parser.add_argument("--tss_bed", "-t", required=True, help="Path to TSS BED file")
    parser.add_argument("--nplb_bed", "-n", required=True, help="Path to nplb_clustered.bed")
    parser.add_argument("--phastcons_bw", "-p", required=True,
                        help="Path to the phastCons bigWig file")
    return parser.parse_args()

def main():
    args = parse_args()
    tss_bed = args.tss_bed
    nplb_bed = args.nplb_bed
    phastcons_bw = args.phastcons_bw

    base_dir = os.path.dirname(nplb_bed)
    work_dir = os.path.join(base_dir, "gene_distances")
    os.makedirs(work_dir, exist_ok=True)
    p_thresh  = 1e-6

    # pipeline
    tss_sorted     = preprocess_tss(tss_bed)
    closest_bed_path = os.path.join(work_dir, "closest_gene.bed")
    closest_bed    = compute_closest(nplb_bed, tss_sorted, closest_bed_path)
    df_bed_path = os.path.join(work_dir, "closest_gene_df.bed")
    df_bed         = filter_columns_bed(closest_bed, [3,4,10,11], df_bed_path)
    hypergeom_df, clusters, p_vals = run_hypergeom_tests(df_bed, p_thresh)
    plot_pvalue_heatmaps(clusters, p_vals, work_dir, p_thresh)
    hypergeom_df.to_csv(os.path.join(work_dir, 'hypergeom_results.tsv'), sep='\t', index=False)

    # Compute average gene distance per cluster
    # Read the dataframe again with Distance column
    data = pd.read_csv(df_bed, sep='\t', header=None)
    data.columns = ['Cluster','motif_direction','gene_direction','Distance']
    avg_dist = data.groupby('Cluster')['Distance'].mean().to_dict()

    # -- Proximal binding count --
    cluster_window = 20000
    proximal_dir = os.path.join(base_dir, 'proximal_binding_count')
    extended_bed = extend_cluster_windows(nplb_bed, proximal_dir, cluster_window)
    count_bed = count_proximal_binding(extended_bed, nplb_bed, proximal_dir)
    avg_binding = compute_average_proximal(count_bed)

    # -- Compute average phastCons per cluster in its own subdir --
    phast_dir = os.path.join(base_dir, 'phastcons')
    # Run bigWigAverageOverBed and compute averages
    phast_tab = run_phastcons_average(phastcons_bw, nplb_bed, phast_dir)
    avg_phast = compute_average_phastcons(phast_tab, nplb_bed)
    
    
if __name__ == "__main__":
    main()
