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
    update_architecture_details(base_dir, args.cluster_map_tsv)

    # Load cluster mapping TSV into a dict: original -> (mapped, flip)
    if args.cluster_map_tsv:
        cluster_map = {}
        with open(args.cluster_map_tsv, newline='') as mapf:
            reader = csv.DictReader(mapf, delimiter='\t')
            for row in reader:
                orig = str(row['original_cluster_id'])
                mapped = row['mapped_cluster_id']
                flip_flag = row['flip'].lower() == 'true'
                if mapped in (None, '', 'None'):
                    mapped_val = None
                else:
                    mapped_val = str(mapped)
                cluster_map[orig] = (mapped_val, flip_flag)
    else:
        # default mapping: identity, no flip
        cluster_map = None


    # -- Parse architectural details if present --
    arch_file = os.path.join(base_dir, 'architectureDetails.txt')
    out_dir = base_dir
    if os.path.exists(arch_file):
        nplb_clustered = pd.read_csv(arch_file, sep='\t', header=None)

        # Prepare DataFrame with duplicate strand column
        df_nplb = pd.DataFrame(columns=[
            'chr', 'start', 'end', 'cluster', 'strand', 'strand_dup'
        ])

        for i in range(nplb_clustered.shape[0]):
            entry = nplb_clustered.iloc[i, 1]
            chrom, coords = entry.split(':', 1)
            parts = coords.split('-')
            start = parts[0]

            # Determine end and initial strand
            if len(parts) < 3:
                end = parts[1][:-3]
                strand = '+'
            else:
                end = parts[1][:-1]
                strand = '-'

            # Remap or skip based on cluster_map
            raw_cluster = str(nplb_clustered.iloc[i, 0])
            if cluster_map is not None:
                # use provided mapping
                if raw_cluster not in cluster_map:
                    continue
                mapped_val, flip_flag = cluster_map[raw_cluster]
                if mapped_val is None:
                    continue
                cluster = mapped_val
            else:
                # identity mapping, no flip
                cluster = raw_cluster
                flip_flag = False

            # Apply flip if requested
            if flip_flag:
                if strand == '+':
                    strand = '-'
                else:
                    strand = '+'

            # Duplicate strand
            strand_dup = strand

            df_nplb.loc[len(df_nplb)] = [
                chrom, start, end, cluster, strand, strand_dup
            ]

        # Save as BED
        bed_out = os.path.join(out_dir, 'nplb_clustered_classified.bed')
        df_nplb.to_csv(bed_out, sep='\t', index=False, header=False)
        print(f"Saved clustered BED to: {bed_out}")

if __name__ == "__main__":
    main()