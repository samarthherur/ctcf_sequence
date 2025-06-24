import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import hypergeom


def main():
    parser = argparse.ArgumentParser(description="CTCF neighbourhood vs atac-seq analysis with ATAC pre-filter")
    parser.add_argument("--nplb_bed", "-n", required=True,
                        help="Path to nplb_clustered.bed")
    parser.add_argument("--atac_bed", "-a", required=True,
                        help="Path to ATAC-seq narrowPeak BED file for the same cell line as nplb clusters")
    parser.add_argument("--atac_dir", "-b", required=True,
                        help="Directory containing atac-seq peak BED files")
    args = parser.parse_args()

    # inputs
    nplb_bed = args.nplb_bed
    atac_bed = args.atac_bed
    atac_dir = args.atac_dir

    # determine output location
    nplb_dir = os.path.dirname(nplb_bed)
    subdir_name = "cell_line_compare_atac"
    output_directory = os.path.join(nplb_dir, subdir_name)
    os.makedirs(output_directory, exist_ok=True)

    # step 0: pre-filter NPLB clusters by ATAC peaks
    atac_filtered_bed = os.path.join(output_directory, f"nplb_atac_filtered.bed")
    os.system(f"bedtools intersect -a {nplb_bed} -b {atac_bed} -u > {atac_filtered_bed}")

    # step 1: intersect nplb (ATAC-filtered) clusters with each atac file
    intersect_files = []
    for peak_file in sorted(os.listdir(atac_dir)):
        if not peak_file.endswith('.bed'):
            continue
        atac_path = os.path.join(atac_dir, peak_file)
        out_path = os.path.join(output_directory, f"nplb_intersect_{peak_file}")
        os.system(f"bedtools intersect -a {atac_filtered_bed} -b {atac_path} -u > {out_path}")
        intersect_files.append(out_path)

    # step 2: load master NPLB (pre-filtered) for stats
    nplb_df = pd.read_csv(atac_filtered_bed, sep="\t", header=None,
                          names=["chr","start","end","cluster","strand"])
    pop_size = len(nplb_df)
    cluster_counts = nplb_df["cluster"].value_counts().sort_index()
    unique_clusters = len(cluster_counts)

    # step 3: prepare result matrices
    p_values = np.zeros((unique_clusters, len(intersect_files)))
    cdf_p_values = np.zeros((unique_clusters, len(intersect_files)))
    sig_matrix = np.zeros((unique_clusters, len(intersect_files)))
    threshold = 0.001

    # step 4: hypergeometric test per cluster per file
    for j, f in enumerate(intersect_files):
        df = pd.read_csv(f, sep="\t", header=None,
                         names=["chr","start","end","cluster","strand"])
        sample_size = len(df)
        obs_counts = df["cluster"].value_counts().to_dict()
        for i, cl in enumerate(sorted(cluster_counts.index)):
            K = cluster_counts.loc[cl]
            k = obs_counts.get(cl, 0)
            p = hypergeom.sf(k-1, pop_size, K, sample_size)
            cdf = hypergeom.cdf(k-1, pop_size, K, sample_size)
            p_values[i, j] = p
            cdf_p_values[i, j] = cdf
            # Bonferroni correction
            bonf = threshold / (unique_clusters * len(intersect_files))
            if p < cdf and p <= bonf:
                sig_matrix[i, j] = 1
            elif cdf <= bonf:
                sig_matrix[i, j] = -1
            else:
                sig_matrix[i, j] = 0

    # step 5: save matrices
    cols = [os.path.basename(x) for x in intersect_files]
    idx = sorted(cluster_counts.index)
    pd.DataFrame(p_values, index=idx, columns=cols).to_csv(
        os.path.join(output_directory, 'p_values_df.csv')
    )
    pd.DataFrame(cdf_p_values, index=idx, columns=cols).to_csv(
        os.path.join(output_directory, 'cdf_p_values_df.csv')
    )
    pd.DataFrame(sig_matrix, index=idx, columns=cols).to_csv(
        os.path.join(output_directory, 'sig_matrix_df.csv')
    )

    # step 6: plot significance heatmap
    plt.figure(figsize=(12, 4))
    sns.heatmap(sig_matrix, cmap='coolwarm', cbar=False,
                xticklabels=cols, yticklabels=idx,
                linewidths=1, linecolor='black')
    plt.xlabel('Cell line')
    plt.ylabel('Cluster')
    plt.title('Significance matrix')
    plt.savefig(os.path.join(output_directory, 'sig_matrix.png'))
    plt.close()

    # step 7: common intersection across all ATAC files
    common_bed = os.path.join(output_directory, "nplb_atac_common.bed")
    prev = atac_filtered_bed
    for peak in sorted(os.listdir(atac_dir)):
        if not peak.endswith('.bed'):
            continue
        bed_path = os.path.join(atac_dir, peak)
        os.system(f"bedtools intersect -a {prev} -b {bed_path} -u > {common_bed}")
        prev = common_bed

    # load common intersections for stats
    df_common = pd.read_csv(common_bed, sep="\t", header=None,
                             names=["chr","start","end","cluster","strand"])
    sample_common = len(df_common)
    obs_common = df_common["cluster"].value_counts().to_dict()

    # prepare common result arrays
    p_common = np.zeros((unique_clusters, 1))
    cdf_common = np.zeros((unique_clusters, 1))
    sig_common = np.zeros((unique_clusters, 1))

    # compute hypergeometric test for common intersection
    for i, cl in enumerate(sorted(cluster_counts.index)):
        K = cluster_counts.loc[cl]
        k = obs_common.get(cl, 0)
        p = hypergeom.sf(k-1, pop_size, K, sample_common)
        cdf = hypergeom.cdf(k-1, pop_size, K, sample_common)
        p_common[i, 0] = p
        cdf_common[i, 0] = cdf
        bonf = threshold / unique_clusters
        if p < cdf and p <= bonf:
            sig_common[i, 0] = 1
        elif cdf <= bonf:
            sig_common[i, 0] = -1
        else:
            sig_common[i, 0] = 0

    # save common significance matrix
    common_col = ["common_atac"]
    pd.DataFrame(p_common, index=sorted(cluster_counts.index), columns=common_col).to_csv(
        os.path.join(output_directory, "common_p_values_df.csv")
    )
    pd.DataFrame(cdf_common, index=sorted(cluster_counts.index), columns=common_col).to_csv(
        os.path.join(output_directory, "common_cdf_values_df.csv")
    )
    pd.DataFrame(sig_common, index=sorted(cluster_counts.index), columns=common_col).to_csv(
        os.path.join(output_directory, "common_sig_matrix_df.csv")
    )

    # plot common heatmap
    plt.figure(figsize=(4, max(5, unique_clusters/2)))
    sns.heatmap(sig_common, cmap='coolwarm', cbar=False,
                xticklabels=common_col, yticklabels=sorted(cluster_counts.index),
                linewidths=1, linecolor='black')
    plt.xlabel('Common ATAC')
    plt.ylabel('Cluster')
    plt.title('Common ATAC Intersection Significance')
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, "common_sig_matrix.png"))
    plt.close()

    print("Analysis complete. Outputs written to:", output_directory)

if __name__ == "__main__":
    main()