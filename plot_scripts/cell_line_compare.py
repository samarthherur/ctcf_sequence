import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

from plot_helper_functions import parse_architectural_details

def main():
    parser = argparse.ArgumentParser(description="CTCF neighbourhood vs ChIP-seq analysis")
    parser.add_argument("--ctcf_neighbourhood_file", "-n", required=True,
                        help="Path to CTCF neighbourhood BED")
    parser.add_argument("--architectural_details", "-a", required=True,
                        help="Path to architectureDetails.txt")
    parser.add_argument("--chip_dir", "-c", required=True,
                        help="Directory containing ChIP-seq peak BED files")
    parser.add_argument("--output_directory", "-o", required=True,
                        help="Directory for outputs")
    args = parser.parse_args()

    input_ctcf_neighbourhood_file = args.ctcf_neighbourhood_file
    architectural_details = args.architectural_details
    chip_dir = args.chip_dir
    output_directory = args.output_directory

    # ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # find all .bed files in chip_dir
    input_chip_files = sorted([
        os.path.join(chip_dir, f)
        for f in os.listdir(chip_dir)
        if f.endswith(".bed")
    ])
    print("Found ChIP BED files:", input_chip_files)

    # extract base filenames
    input_just_filenames = [os.path.basename(f) for f in input_chip_files]

    # intersect each ChIP with the CTCF neighbourhoods
    for file, filename in zip(input_chip_files, input_just_filenames):
        out_path = os.path.join(output_directory, f"ctcf_neighbourhood_h1_intersect_{filename}")
        cmd = f"bedtools intersect -a {input_ctcf_neighbourhood_file} -b {file} -u > {out_path}"
        os.system(cmd)

    # collect the intersected files
    input_intersect_files = sorted([
        os.path.join(output_directory, f)
        for f in os.listdir(output_directory)
        if f.startswith("ctcf_neighbourhood_h1_intersect_")
    ])
    print("Intersected files:", input_intersect_files)

    # load CTCF neighbourhoods
    ctcf_neighbourhood = pd.read_csv(input_ctcf_neighbourhood_file, sep="\t", header=None)
    N = ctcf_neighbourhood.shape[0]

    # build overlap matrix
    neighbourhood_matrix = np.zeros((N, len(input_intersect_files)), dtype=int)
    intersect_dfs = []
    for idx, file in enumerate(input_intersect_files):
        df = pd.read_csv(file, sep="\t", header=None)
        print(df.shape, "rows in", file)
        intersect_dfs.append(df)
        for _, row in df.iterrows():
            matches = ctcf_neighbourhood[
                (ctcf_neighbourhood[0]==row[0]) &
                (ctcf_neighbourhood[1]==row[1]) &
                (ctcf_neighbourhood[2]==row[2])
            ].index
            neighbourhood_matrix[matches, idx] = 1

    # make DataFrame and add sum
    neighbourhood_df = pd.DataFrame(neighbourhood_matrix, columns=input_just_filenames)
    neighbourhood_df['sum'] = neighbourhood_df.sum(axis=1)

    # plot overlap distribution
    sns.barplot(
        x=neighbourhood_df['sum'].value_counts().index,
        y=neighbourhood_df['sum'].value_counts().values
    )
    plt.xlabel('Number of overlaps')
    plt.ylabel('Frequency')
    plt.title('Number of overlaps between CTCF neighbourhood and ChIP-seq peaks')
    plt.savefig(os.path.join(output_directory, 'overlap_distribution.png'))
    plt.close()

    # plot overlaps per sample (first four)
    sum2 = neighbourhood_df.sum(axis=0)
    plt.figure(figsize=(12,4))
    plt.xlabel('Cell line')
    plt.ylabel('Number of overlaps')
    sns.barplot(x=sum2.index[:4], y=sum2.values[:4])
    plt.savefig(os.path.join(output_directory, 'overlap_by_cell_line.png'))
    plt.close()

    # save neighbourhood matrix
    neighbourhood_df.to_csv(os.path.join(output_directory, 'neighbourhood_df.csv'), index=False)

    # load and parse architectural details
    architectural = pd.read_csv(architectural_details, sep="\t", header=None)
    parsed = parse_architectural_details(architectural)

    # prepare intersect_dfs for clustering
    for df in intersect_dfs:
        df.sort_index(inplace=True)
        df.rename(columns={0:'chr',1:'start',2:'end'}, inplace=True)
        df.drop(columns=['cluster'], errors='ignore', inplace=True)

    def ensure_dtype_consistency(df):
        df['chr']   = df['chr'].astype(str)
        df['start'] = df['start'].astype(int)
        df['end']   = df['end'].astype(int)
        return df

    def assign_clusters(cluster_df, intersect_dfs):
        updated = []
        for df in intersect_dfs:
            merged = pd.merge(df, cluster_df, on=['chr','start','end'], how='left')
            merged['cluster'] = pd.to_numeric(merged['cluster'], errors='coerce').astype('Int64')
            updated.append(merged)
        return updated

    parsed = ensure_dtype_consistency(parsed)
    new_intersect_dfs = assign_clusters(parsed, intersect_dfs)
    new_intersect_dfs = [df.dropna(subset=['cluster']) for df in new_intersect_dfs]

    # enrichment statistics
    pop_size = len(parsed)
    unique = parsed['cluster'].nunique()
    p_values = np.zeros((unique, len(new_intersect_dfs)))
    cdf_p_values = np.zeros((unique, len(new_intersect_dfs)))
    sig_matrix = np.zeros((unique, len(new_intersect_dfs)))
    threshold = 0.001

    for i in range(1, unique+1):
        pop_white = (parsed['cluster']==i).sum()
        for j, df in enumerate(new_intersect_dfs):
            sample_size = len(df)
            sample_white = (df['cluster']==i).sum()
            p = stats.hypergeom.sf(sample_white-1, pop_size, pop_white, sample_size)
            cdf = stats.hypergeom.cdf(sample_white-1, pop_size, pop_white, sample_size)
            p_values[i-1,j] = p
            cdf_p_values[i-1,j] = cdf
            if p < cdf and p <= (threshold/(unique*len(new_intersect_dfs))):
                sig_matrix[i-1,j] = 1
            elif cdf <= (threshold/(unique*len(new_intersect_dfs))):
                sig_matrix[i-1,j] = -1
            else:
                sig_matrix[i-1,j] = 0

    # save matrices
    pd.DataFrame(p_values).to_csv(os.path.join(output_directory, 'p_values_df.csv'), index=False)
    pd.DataFrame(cdf_p_values).to_csv(os.path.join(output_directory, 'cdf_p_values_df.csv'), index=False)
    pd.DataFrame(sig_matrix).to_csv(os.path.join(output_directory, 'sig_matrix_df.csv'), index=False)

    # plot significance heatmap
    plt.figure(figsize=(12,4))
    sns.heatmap(sig_matrix, cmap='coolwarm', cbar=False, linewidths=1, linecolor='black')
    plt.xlabel('Cell line')
    plt.ylabel('Cluster')
    plt.title('Significance matrix')
    plt.yticks(ticks=np.arange(unique)+0.5, labels=np.arange(1, unique+1))
    plt.savefig(os.path.join(output_directory, 'sig_matrix.png'))
    plt.close()

    print("Analysis complete. Outputs written to:", output_directory)

if __name__ == "__main__":
    main()