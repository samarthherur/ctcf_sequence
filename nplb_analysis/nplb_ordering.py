# nplb_classify.py

import os
import json
import csv
import pandas as pd
from nplb_helper_functions import *

def update_architecture_details():
    pass

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
    args = parse_nplb_classify_args()
    fasta         = args.fasta
    model         = args.model
    output_prefix = args.output_prefix

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

    # 1) Run promoterClassify
    run_promoter_classify(fasta, model, output_prefix)

    # -- Parse architectural details if present --
    arch_file = os.path.join(os.path.dirname(output_prefix), 'architectureDetails.txt')
    
    if os.path.exists(arch_file):
        nplb_clustered = pd.read_csv(arch_file, sep='\t', header=None)
        out_dir = os.path.dirname(arch_file)

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
        bed_out = os.path.join(out_dir, 'nplb_clustered_updated.bed')
        df_nplb.to_csv(bed_out, sep='\t', index=False, header=False)
        print(f"Saved clustered BED to: {bed_out}")

        # Update the NPLB clustered bed using the same mapping TSV
        update_nplb_bed(out_dir, args.cluster_map_tsv)

if __name__ == "__main__":
    main()