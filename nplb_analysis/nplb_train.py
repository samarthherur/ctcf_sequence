import os
import pandas as pd
from nplb_helper_functions import *

def main():
    args = parse_nplb_train_args()
    fasta = args.fasta
    output_prefix = args.output_prefix

    # # 1) Run promoterLearn
    # run_promoter_learn(fasta, output_prefix)
    
    os.system("module unload python-2.7")

    # -- Parse architectural details if file exists in output directory --
    arch_file = os.path.join(os.path.dirname(output_prefix),
                             'architectureDetails.txt')
    if os.path.exists(arch_file):
        # Load raw architecture details
        nplb_clustered = pd.read_csv(arch_file, sep='\t', header=None)
        out_dir = os.path.dirname(arch_file)

        # Build DataFrame of chr, start, end, cluster, strand
        df_nplb = pd.DataFrame(columns=['chr','start','end','cluster','strand'])
        for i in range(nplb_clustered.shape[0]):
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

        # Save as BED
        bed_out = os.path.join(out_dir, 'nplb_clustered.bed')
        df_nplb.to_csv(bed_out, sep='\t', index=False, header=False)
        print(f"Saved clustered BED to: {bed_out}")

if __name__ == "__main__":
    main()