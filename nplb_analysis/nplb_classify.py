# nplb_classify.py

import os
import json
import csv
import pandas as pd
import argparse
from nplb_helper_functions import *

def parse_nplb_classify_args():
    """
    Parse arguments for classification: input FASTA, model, output prefix,
    and a TSV mapping file with columns 'original_cluster_id', 'mapped_cluster_id', 'flip'.
    """
    parser = argparse.ArgumentParser(
        description="Run promoterClassify on neighbourhood FASTA with cluster mapping"
    )
    parser.add_argument(
        "--nplb_classify_dir", "-o",
        dest="nplb_classify_dir",
        required=True,
        help="Directory where promoterClassify output resides"
    )
    parser.add_argument(
        "--cluster_map_tsv",
        required=False,
        default=None,
        help="Optional TSV file mapping original_cluster_id to mapped_cluster_id, with a 'flip' column"
    )
    return parser.parse_args()

def main():
    args = parse_nplb_classify_args()
    # Update architectureDetails.txt before classification
    base_dir = args.nplb_classify_dir
    
    # -- Build nplb_clustered.bed from architectureDetails.txt --
    arch_file = os.path.join(base_dir, 'architectureDetails.txt')
    assert os.path.exists(arch_file), \
        f"Architecture details file not found: {arch_file}. Please run promoterLearn first."
    # Parse architecture details into a DataFrame
    nplb_clustered = pd.read_csv(arch_file, sep='\t', header=None)
    df_nplb = pd.DataFrame(columns=['chr', 'start', 'end', 'cluster', 'strand'])
    for i in range(len(nplb_clustered)):
        entry = nplb_clustered.iloc[i, 1]
        chrom, coords = entry.split(':', 1)
        parts = coords.split('-')
        start = parts[0]
        if len(parts) < 3:
            end = parts[1][:-3]
            strand = '+'
        else:
            end = parts[1][:-1]
            strand = '-'
        cluster = nplb_clustered.iloc[i, 0]
        df_nplb.loc[i] = [chrom, start, end, cluster, strand]
    bed_path = os.path.join(base_dir, 'nplb_clustered.bed')
    df_nplb.to_csv(bed_path, sep='\t', index=False, header=False)
    print(f"Saved clustered BED to: {bed_path}")
    
    update_architecture_details(base_dir, args.cluster_map_tsv)

    # Remap and flip BED according to the provided mapping TSV
    update_nplb_bed(base_dir, args.cluster_map_tsv)

if __name__ == "__main__":
    main()