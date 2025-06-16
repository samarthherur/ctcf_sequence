import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import json


def parse_nplb_train_args():
    """
    Parse arguments for training: only input FASTA and output prefix.
    """
    parser = argparse.ArgumentParser(
        description="Run promoterLearn training on neighbourhood FASTA"
    )
    parser.add_argument(
        "--fasta", "-f",
        required=True,
        help="Input FASTA file (neighbourhood sequences) for training"
    )
    parser.add_argument(
        "--output_prefix", "-o",
        required=True,
        help="Prefix (directory+basename) for all training outputs"
    )
    return parser.parse_args()

def parse_nplb_classify_args():
    """
    Parse arguments for classification: input FASTA, model, output prefix,
    and a TSV mapping file with columns 'original_cluster_id', 'mapped_cluster_id', 'flip'.
    """
    parser = argparse.ArgumentParser(
        description="Run promoterClassify on neighbourhood FASTA with cluster mapping"
    )
    parser.add_argument(
        "--fasta", "-f",
        required=True,
        help="Input FASTA file (neighbourhood sequences) for classification"
    )
    parser.add_argument(
        "--model", "-m",
        required=True,
        help="Path to the trained promoterClassify model (.p file)"
    )
    parser.add_argument(
        "--output_prefix", "-o",
        required=True,
        help="Prefix (directory+basename) for all classification outputs"
    )
    parser.add_argument(
        "--cluster_map_tsv",
        required=False,
        default=None,
        help="Optional TSV file mapping original_cluster_id to mapped_cluster_id, with a 'flip' column"
    )
    return parser.parse_args()

def run_promoter_learn(fasta_path: str, output_prefix: str):
    cmd = f"promoterLearn -f {fasta_path} -o {output_prefix}"
    print(f"Executing: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"promoterLearn failed with exit code {ret}")

def run_promoter_classify(fasta_path: str, model_path: str, output_prefix: str):
    """
    Run promoterClassify on the given FASTA using a pre-trained model.
    """
    cmd = f"promoterClassify -f {fasta_path} -m {model_path} -o {output_prefix}"
    print(f"Executing: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"promoterClassify failed with exit code {ret}")
    
    
###ordering functions
def preprocess_tss(in_bed: str) -> str:
    """Drop duplicates and sort TSS BED; return sorted path."""
    out = in_bed.replace('.bed', '_no_duplicates.bed')
    df = pd.read_csv(in_bed, sep='\t', header=None)
    # drop duplicates on columns 1 and 2
    df = df.drop_duplicates(subset=1).drop_duplicates(subset=2)
    df = df.sort_values(by=5)
    df.to_csv(out, sep='\t', header=False, index=False)
    sorted_out = out.replace('.bed', '_sorted.bed')
    os.system(f'bedtools sort -i {out} > {sorted_out}')
    return sorted_out

def compute_closest(nplb_bed: str, tss_sorted: str, out_bed: str) -> str:
    """Sort nplb_bed and compute closest to TSS; return path."""
    cmd = f'bedtools sort -i {nplb_bed} | bedtools closest -a - -b {tss_sorted} -D ref -t first > {out_bed}'
    os.system(cmd)
    return out_bed

def filter_columns_bed(in_bed: str, cols: list, out_bed: str) -> str:
    """Select only columns by zero-based index in cols."""
    df = pd.read_csv(in_bed, sep='\t', header=None)
    df = df.iloc[:, cols]
    df.to_csv(out_bed, sep='\t', header=False, index=False)
    return out_bed

def run_hypergeom_tests(df_bed: str, p_thresh: float):
    """Perform hypergeometric test per cluster; return DataFrame, clusters array, p-values."""
    data = pd.read_csv(df_bed, sep='\t', header=None)
    data.columns = ['Cluster','motif_direction','gene_direction','Distance']
    total = len(data)
    clusters = np.sort(data['Cluster'].unique())
    zero = data[data['Distance']==0]
    total_zero = len(zero)

    results = []
    for c in clusters:
        n_cluster = len(data[data['Cluster']==c])
        n_zero_cluster = len(zero[zero['Cluster']==c])
        pval = hypergeom.sf(n_zero_cluster-1, total, n_cluster, total_zero)
        results.append({
            'Cluster': c,
            'Total Entries in Cluster': n_cluster,
            'Total Zero Distance Entries': total_zero,
            'Zero Distance Entries in Cluster': n_zero_cluster,
            'Hypergeometric P-Value': pval
        })
    df_res = pd.DataFrame(results)
    p_vals = df_res['Hypergeometric P-Value'].values
    return df_res, clusters, p_vals

def plot_pvalue_heatmaps(clusters, p_vals, work_dir: str, p_thresh: float):
    """Plot raw, log, and color-coded heatmaps of p-values."""
    # raw
    fig, ax = plt.subplots()
    cax = ax.matshow(p_vals.reshape(1,-1), cmap='coolwarm', aspect='auto')
    fig.colorbar(cax)
    ax.set_xticks(range(len(clusters))); ax.set_xticklabels(clusters)
    ax.set_yticks([])
    plt.title('Hypergeometric P-Values')
    fig.savefig(os.path.join(work_dir, 'p_values_heatmap.png'))
    plt.close(fig)

    # log scale
    fig, ax = plt.subplots()
    cax = ax.matshow(-np.log10(p_vals).reshape(1,-1), cmap='coolwarm', aspect='auto')
    fig.colorbar(cax)
    ax.set_xticks(range(len(clusters))); ax.set_xticklabels(clusters)
    ax.set_yticks([])
    plt.title('-log10(P-Values)')
    fig.savefig(os.path.join(work_dir, 'log_p_values_heatmap.png'))
    plt.close(fig)

    # color-coded
    sig = clusters[p_vals < p_thresh]
    opp = clusters[p_vals > 1 - p_thresh]
    colors = ['red' if c in sig else 'blue' if c in opp else 'grey' for c in clusters]
    fig, ax = plt.subplots()
    cax = ax.matshow(p_vals.reshape(1,-1), cmap='coolwarm', aspect='auto')
    fig.colorbar(cax)
    ax.set_xticks(range(len(clusters))); ax.set_xticklabels(clusters)
    for i, lbl in enumerate(ax.get_xticklabels()):
        lbl.set_color(colors[i])
    ax.set_yticks([])
    plt.title('Color-coded P-Values')
    fig.savefig(os.path.join(work_dir, 'color_p_values_heatmap.png'))
    plt.close(fig)
    
# --- Helper functions for proximal binding ---
def extend_cluster_windows(bed_file: str, output_dir: str, window: int) -> str:
    """
    Extend each interval in the BED by +/- window bp, floor at zero.
    Returns path to the extended BED.
    """
    import pandas as pd, os
    df = pd.read_csv(bed_file, sep='\t', header=None)
    df[1] = df[1] - window
    df[2] = df[2] + window
    df.loc[df[1] < 0, 1] = 0
    df.loc[df[2] < 0, 2] = 0
    out_bed = os.path.join(output_dir, 'nplb_clustered_extended.bed')
    df.to_csv(out_bed, sep='\t', header=False, index=False)
    return out_bed

def count_proximal_binding(extended_bed: str, original_bed: str, output_dir: str) -> str:
    """
    Count overlaps of extended regions with original BED.
    Returns path to the cluster-count BED.
    """
    import os, subprocess
    os.makedirs(output_dir, exist_ok=True)
    out_bed = os.path.join(output_dir, 'cluster_count.bed')
    cmd = f'bedtools intersect -a {extended_bed} -b {original_bed} -c > {out_bed}'
    subprocess.call(cmd, shell=True)
    return out_bed

def compute_average_proximal(count_bed: str) -> dict:
    """
    Read the cluster-count BED and compute average binding count per cluster.
    Returns a dict {cluster_id: avg_count}.
    """
    import pandas as pd
    df = pd.read_csv(count_bed, sep='\t', header=None)
    df.columns = ['chr','start','end','cluster_id','strand','cluster_count']
    return df.groupby('cluster_id')['cluster_count'].mean().to_dict()


##phastcons
def run_phastcons_average(phastcons_bw: str, bed_file: str, output_dir: str) -> str:
    """
    Run bigWigAverageOverBed to compute per-region phastCons scores.
    Returns the path to the output table.
    """
    import os
    import subprocess

    os.makedirs(output_dir, exist_ok=True)
    out_tab = os.path.join(output_dir, 'phastcons_over_clusters.tab')
    cmd = f"bigWigAverageOverBed {phastcons_bw} {bed_file} {out_tab}"
    print(f"Executing: {cmd}")
    ret = subprocess.call(cmd, shell=True)
    if ret != 0:
        raise RuntimeError(f"bigWigAverageOverBed failed with exit code {ret}")
    return out_tab


def compute_average_phastcons(phast_tab: str, bed_file: str) -> dict:
    """
    Read bigWigAverageOverBed output and compute average phastCons per cluster.
    Returns a dict mapping cluster_id → average phastCons score.
    """
    import pandas as pd

    # Load the phastCons summary table
    df = pd.read_csv(phast_tab, sep='\t', header=None)

    # Reload the original BED to get cluster IDs in the same order
    df_regions = pd.read_csv(
        bed_file,
        sep='\t',
        header=None,
        names=['chr', 'start', 'end', 'cluster_id', 'strand']
    )

    # Append the cluster_id column
    df['cluster_id'] = df_regions['cluster_id']

    # Column index 4 (0-based) contains the "mean0" value from bigWigAverageOverBed
    avg_phast = df.groupby('cluster_id')[4].mean().to_dict()

    return avg_phast

# --- Function to update architectureDetails.txt based on mapping TSV ---
def update_architecture_details(base_dir: str, cluster_map_tsv: str):
    """
    Reads architectureDetails.txt in base_dir, applies the cluster_map_tsv mapping,
    and writes architectureDetails_updated.txt.
    """
    import os, csv, pandas as pd

    if not cluster_map_tsv:
        return

    # Load mapping from TSV
    cluster_mapping = {}
    with open(cluster_map_tsv, newline='') as mf:
        reader = csv.DictReader(mf, delimiter='\t')
        for row in reader:
            orig = str(row['original_cluster_id'])
            mapped = row['mapped_cluster_id']
            if mapped in (None, '', 'None'):
                mapped_val = None
            else:
                mapped_val = int(mapped)
            cluster_mapping[orig] = mapped_val

    arch_path = os.path.join(base_dir, 'architectureDetails.txt')
    if os.path.exists(arch_path):
        arch = pd.read_csv(arch_path, sep='\t', header=None)
        arch.iloc[:,0] = arch.iloc[:,0].astype(str).map(cluster_mapping)
        arch = arch.dropna(subset=[0])
        arch.iloc[:,0] = arch.iloc[:,0].astype(int)
        updated_path = os.path.join(base_dir, 'architectureDetails_updated.txt')
        arch.to_csv(updated_path, sep='\t', header=False, index=False)
        print(f"Saved updated architecture details to {updated_path}")
        
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