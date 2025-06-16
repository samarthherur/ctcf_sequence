import os
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import argparse
import json
import subprocess
import csv

from nplb_helper_functions import *

def parse_args():
    parser = argparse.ArgumentParser(description="Compute avg gene distances per cluster")
    parser.add_argument("--tss_bed", "-t", required=False, default=None,
                        help="(Optional) Path to TSS BED file")
    parser.add_argument("--nplb_bed", "-n", required=False, default=None,
                        help="(Optional) Path to nplb_clustered.bed")
    parser.add_argument("--phastcons_bw", "-p", required=False, default=None,
                        help="(Optional) Path to the phastCons bigWig file")
    parser.add_argument(
        "--metric", "-m",
        choices=["avg_dist","avg_binding","avg_phastcons","hybrid","none"],
        default="none",
        help="Metric to order clusters by"
    )
    parser.add_argument(
        "--delete_clusters", "-d",
        type=str,
        default="",
        help="Comma-separated list of cluster IDs to delete (mapped to None)"
    )
    parser.add_argument(
        "--flip_clusters", "-f",
        type=str,
        default="",
        help="Comma-separated list of cluster IDs whose strand should be flipped"
    )
    parser.add_argument(
        "--cluster_map_tsv",
        required=False,
        default=None,
        help="Optional path to cluster mapping TSV"
    )
    return parser.parse_args()

def main():
    args = parse_args()
    tss_bed = args.tss_bed
    nplb_bed = args.nplb_bed
    phastcons_bw = args.phastcons_bw

    # Parse new CLI flags
    metric = args.metric
    delete_clusters = [s.strip() for s in args.delete_clusters.split(",") if s.strip()]
    flip_clusters = [s.strip() for s in args.flip_clusters.split(",") if s.strip()]

    # Validate metric requirements
    import sys
    if metric in ("avg_dist","hybrid") and not tss_bed:
        sys.exit("ERROR: --metric %s requires --tss_bed" % metric)
    if metric in ("avg_binding","hybrid") and not nplb_bed:
        sys.exit("ERROR: --metric %s requires --nplb_bed" % metric)
    if metric == "avg_phastcons" and not phastcons_bw:
        sys.exit("ERROR: --metric %s requires --phastcons_bw" % metric)
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

    # ---- Metric mapping and ordering ----
    # clusters is a list of cluster IDs (as str)
    metrics_df = pd.DataFrame({
        'Cluster': clusters,
        'avg_dist': [avg_dist.get(c) for c in clusters],
        'avg_binding': [avg_binding.get(int(c)) for c in clusters],
        'avg_phastcons': [avg_phast.get(int(c)) for c in clusters]
    })
    # Hybrid is average of normalized avg_dist and avg_binding
    if metric == "hybrid":
        md = metrics_df['avg_dist']
        mb = metrics_df['avg_binding']
        metrics_df['norm_dist'] = 1 - ((md - md.min())/(md.max()-md.min()))
        metrics_df['norm_binding'] = (mb - mb.min())/(mb.max()-mb.min())
        metrics_df['hybrid'] = (metrics_df['norm_dist'] + metrics_df['norm_binding'])/2
    # Determine ordering values
    if metric in ["avg_dist","avg_binding","avg_phastcons","hybrid"]:
        metrics_df['order_val'] = metrics_df[metric]
    else:
        metrics_df['order_val'] = None
    # Rank clusters (higher order_val = higher rank; for avg_dist lower is better so invert)
    if metric == "avg_dist":
        metrics_df = metrics_df.sort_values('order_val', ascending=True)
    elif metric in ["avg_binding","avg_phastcons","hybrid"]:
        metrics_df = metrics_df.sort_values('order_val', ascending=False)

    # Build mapping table
    map_df = pd.DataFrame({
        'original_cluster_id': metrics_df['Cluster'].astype(int).astype(str),
        'mapped_cluster_id': metrics_df['Cluster'].astype(int).astype(str),
        'flip': False,
        'rank': list(range(1, len(metrics_df)+1))
    })
    # Apply deletions
    map_df.loc[map_df['original_cluster_id'].isin(delete_clusters), ['mapped_cluster_id','rank']] = (None, None)
    # Apply flips
    map_df.loc[map_df['original_cluster_id'].isin(flip_clusters), 'flip'] = True

    # Save mapping to tsv in base_dir
    map_out = os.path.join(base_dir, "cluster_mapping.tsv")
    map_df.to_csv(map_out, sep='\t', index=False)
    print(f"Saved cluster mapping to {map_out}")

    # Update architecture details file if mapping TSV provided
    update_architecture_details(base_dir, args.cluster_map_tsv)
    
if __name__ == "__main__":
    main()
