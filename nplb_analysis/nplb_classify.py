# nplb_classify.py

import os
import json
import csv
import pandas as pd
from nplb_helper_functions import parse_nplb_classify_args, run_promoter_classify

def main():
    args = parse_nplb_classify_args()
    fasta         = args.fasta
    model         = args.model
    output_prefix = args.output_prefix

    # Load cluster mapping TSV into a dict: original -> (mapped, flip)
    cluster_map = {}
    with open(args.cluster_map_tsv, newline='') as mapf:
        reader = csv.DictReader(mapf, delimiter='\t')
        for row in reader:
            orig = str(row['original_cluster_id'])
            mapped = row['mapped_cluster_id']
            flip_flag = row['flip'].lower() == 'true'
            # Normalize mapped values
            if mapped in (None, '', 'None'):
                mapped_val = None
            else:
                mapped_val = str(mapped)
            cluster_map[orig] = (mapped_val, flip_flag)

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
            if raw_cluster not in cluster_map:
                continue
            mapped_val, flip_flag = cluster_map[raw_cluster]
            if mapped_val is None:
                continue
            cluster = mapped_val

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