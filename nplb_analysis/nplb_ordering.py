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

def update_nplb_bed(base_dir: str, cluster_map_tsv: str):
    """
    Reads nplb_clustered.bed in base_dir, applies mapping TSV (drops None mappings,
    flips strand if flagged), and writes nplb_clustered_updated.bed.
    """
    import os, csv, pandas as pd

    bed_in = os.path.join(base_dir, 'nplb_clustered.bed')
    if not (os.path.exists(bed_in) and cluster_map_tsv and os.path.exists(cluster_map_tsv)):
        return

    # Load mapping with flip
    cluster_map = {}
    with open(cluster_map_tsv, newline='') as mf:
        reader = csv.DictReader(mf, delimiter='\t')
        for row in reader:
            orig = str(row['original_cluster_id'])
            mapped = row['mapped_cluster_id']
            flip_flag = row['flip'].lower() == 'true'
            if mapped in (None, '', 'None'):
                mapped_val = None
            else:
                mapped_val = str(mapped)
            cluster_map[orig] = (mapped_val, flip_flag)

    # Read and transform BED
    df = pd.read_csv(bed_in, sep='\t', header=None)
    new_rows = []
    for _, row in df.iterrows():
        orig = str(row[3])
        if orig not in cluster_map:
            continue
        mapped_val, flip_flag = cluster_map[orig]
        if mapped_val is None:
            continue
        strand = row[4]
        if flip_flag:
            strand = '+' if strand == '-' else '-'
        new_rows.append([row[0], row[1], row[2], mapped_val, strand])

    if not new_rows:
        return

    df_out = pd.DataFrame(new_rows, columns=[0,1,2,3,4])
    out_bed = os.path.join(base_dir, 'nplb_clustered_updated.bed')
    df_out.to_csv(out_bed, sep='\t', header=False, index=False)
    print(f"Saved updated NPLB BED to {out_bed}")


def main():
    args = parse_args()
    tss_bed = args.tss_bed
    nplb_bed = args.nplb_bed
    phastcons_bw = args.phastcons_bw

    metric = args.metric
    # Enforce metric none when using an external cluster_map_tsv
    if args.cluster_map_tsv and metric != "none":
        import sys
        sys.exit("ERROR: --cluster_map_tsv requires --metric none")
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

    # pipeline: conditional steps based on provided files
    clusters = []
    avg_dist = {}
    avg_binding = {}
    avg_phast = {}

    # 1) Gene distance & hypergeom (requires both tss_bed and nplb_bed)
    if tss_bed and nplb_bed:
        tss_sorted = preprocess_tss(tss_bed)
        closest_bed_path = os.path.join(work_dir, "closest_gene.bed")
        closest_bed = compute_closest(nplb_bed, tss_sorted, closest_bed_path)
        df_bed_path = os.path.join(work_dir, "closest_gene_df.bed")
        df_bed = filter_columns_bed(closest_bed, [3,4,10,11], df_bed_path)
        hypergeom_df, clusters, p_vals = run_hypergeom_tests(df_bed, p_thresh)
        plot_pvalue_heatmaps(clusters, p_vals, work_dir, p_thresh)
        hypergeom_df.to_csv(os.path.join(work_dir, 'hypergeom_results.tsv'), sep='\t', index=False)
        data = pd.read_csv(df_bed, sep='\t', header=None)
        data.columns = ['Cluster','motif_direction','gene_direction','Distance']
        avg_dist = data.groupby('Cluster')['Distance'].mean().to_dict()
    elif nplb_bed:
        # derive cluster list if no TSS analysis
        df_nb = pd.read_csv(nplb_bed, sep='\t', header=None)
        clusters = sorted(df_nb[3].astype(str).unique())

    # 2) Proximal binding count (requires nplb_bed)
    if nplb_bed:
        cluster_window = 20000
        proximal_dir = os.path.join(base_dir, 'proximal_binding_count')
        os.makedirs(proximal_dir, exist_ok=True)
        extended_bed = extend_cluster_windows(nplb_bed, proximal_dir, cluster_window)
        count_bed = count_proximal_binding(extended_bed, nplb_bed, proximal_dir)
        avg_binding = compute_average_proximal(count_bed)
    else:
        avg_binding = {c: None for c in clusters}

    # 3) phastCons average (requires phastcons_bw)
    if phastcons_bw:
        phast_dir = os.path.join(base_dir, 'phastcons')
        os.makedirs(phast_dir, exist_ok=True)
        phast_tab = run_phastcons_average(phastcons_bw, nplb_bed, phast_dir)
        avg_phast = compute_average_phastcons(phast_tab, nplb_bed)
    else:
        avg_phast = {c: None for c in clusters}

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
    elif metric in ["avg_binding", "avg_phastcons", "hybrid"]:
        metrics_df = metrics_df.sort_values('order_val', ascending=False)

    # Build mapping table: original clusters, apply deletions and flips, then assign continuous mapped IDs
    map_df = pd.DataFrame({
        'original_cluster_id': metrics_df['Cluster'].astype(int).astype(str)
    })

    # Apply deletions: drop any original IDs in delete_clusters
    if delete_clusters:
        map_df = map_df[~map_df['original_cluster_id'].isin(delete_clusters)].copy()

    # Initialize flip column, then apply flips based on original clusters
    map_df['flip'] = False
    if flip_clusters:
        mask_flip = map_df['original_cluster_id'].isin(flip_clusters)
        map_df.loc[mask_flip, 'flip'] = True

    # Assign continuous mapped_cluster_id starting from 1 in the sorted order
    map_df = map_df.reset_index(drop=True)
    map_df['mapped_cluster_id'] = map_df.index + 1

    # Reorder columns to original_cluster_id, mapped_cluster_id, flip
    map_df = map_df[['original_cluster_id', 'mapped_cluster_id', 'flip']]

    # Save mapping to TSV
    map_out = os.path.join(base_dir, "cluster_mapping.tsv")
    map_df.to_csv(map_out, sep='\t', index=False)
    print(f"Saved cluster mapping to {map_out}")

    # Choose mapping TSV: use provided one if available, else use computed map_out
    cluster_map_path = args.cluster_map_tsv if args.cluster_map_tsv else map_out
    update_architecture_details(base_dir, cluster_map_path)
    update_nplb_bed(base_dir, cluster_map_path)

if __name__ == "__main__":
    main()