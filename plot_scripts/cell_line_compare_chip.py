import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats


def main():
    parser = argparse.ArgumentParser(description="CTCF neighbourhood vs ChIP-seq analysis")
    parser.add_argument("--nplb_bed", "-n", required=True,
                        help="Path to nplb_clustered.bed")
    parser.add_argument("--chip_dir", "-c", required=True,
                        help="Directory containing ChIP-seq peak BED files")
    args = parser.parse_args()

    nplb_bed = args.nplb_bed
    chip_dir = args.chip_dir

    nplb_dir = os.path.dirname(nplb_bed)
    subdir_name = "cell_line_compare_chip"
    output_directory = os.path.join(nplb_dir, subdir_name)

    # ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # find all .bed files in chip_dir
    input_chip_files = sorted([
        os.path.join(chip_dir, f)
        for f in os.listdir(chip_dir)
        if f.endswith(".bed")
    ])

    # intersect nplb clusters with each ChIP file
    intersect_files = []
    for peak_file in sorted([os.path.join(chip_dir,f) for f in os.listdir(chip_dir) if f.endswith('.bed')]):
        out_path = os.path.join(output_directory, f"nplb_intersect_{os.path.basename(peak_file)}")
        os.system(f"bedtools intersect -a {nplb_bed} -b {peak_file} -u > {out_path}")
        intersect_files.append(out_path)

    # load nplb clustered BED
    nplb_df = pd.read_csv(nplb_bed, sep="\t", header=None,
                          names=["chr","start","end","cluster","strand"])
    pop_size = len(nplb_df)
    cluster_counts = nplb_df["cluster"].value_counts().sort_index()
    unique_clusters = cluster_counts.shape[0]

    from scipy.stats import hypergeom

    # prepare result matrices
    p_values = np.zeros((unique_clusters, len(intersect_files)))
    cdf_p_values = np.zeros((unique_clusters, len(intersect_files)))
    sig_matrix = np.zeros((unique_clusters, len(intersect_files)))
    threshold = 0.001

    # for each cell line (intersect file)
    for j, file in enumerate(intersect_files):
        df = pd.read_csv(file, sep="\t", header=None,
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
            if p < cdf and p <= (threshold/(unique_clusters*len(intersect_files))):
                sig_matrix[i, j] = 1
            elif cdf <= (threshold/(unique_clusters*len(intersect_files))):
                sig_matrix[i, j] = -1
            else:
                sig_matrix[i, j] = 0

    pd.DataFrame(p_values, index=sorted(cluster_counts.index),
                 columns=[os.path.basename(f) for f in intersect_files]) \
         .to_csv(os.path.join(output_directory, 'p_values_df.csv'))
    pd.DataFrame(cdf_p_values, index=sorted(cluster_counts.index),
                 columns=[os.path.basename(f) for f in intersect_files]) \
         .to_csv(os.path.join(output_directory, 'cdf_p_values_df.csv'))
    pd.DataFrame(sig_matrix, index=sorted(cluster_counts.index),
                 columns=[os.path.basename(f) for f in intersect_files]) \
         .to_csv(os.path.join(output_directory, 'sig_matrix_df.csv'))

    plt.figure(figsize=(12,4))
    sns.heatmap(sig_matrix, cmap='coolwarm', cbar=False,
                xticklabels=[os.path.basename(f) for f in intersect_files],
                yticklabels=sorted(cluster_counts.index), linewidths=1, linecolor='black')
    plt.xlabel('Cell line')
    plt.ylabel('Cluster')
    plt.title('Significance matrix')
    plt.savefig(os.path.join(output_directory, 'sig_matrix.png'))
    plt.close()

    # common intersection across all ChIP files
    common_chip_bed = os.path.join(output_directory, "nplb_chip_common.bed")
    prev = nplb_bed
    for peak in sorted(os.listdir(chip_dir)):
        if not peak.endswith('.bed'):
            continue
        chip_path = os.path.join(chip_dir, peak)
        os.system(f"bedtools intersect -a {prev} -b {chip_path} -u > {common_chip_bed}")
        prev = common_chip_bed

    # load common chip intersections for stats
    df_common_chip = pd.read_csv(common_chip_bed, sep="\t", header=None,
                                 names=["chr","start","end","cluster","strand"])
    sample_common_chip = len(df_common_chip)
    obs_common_chip = df_common_chip["cluster"].value_counts().to_dict()

    # prepare common chip result arrays
    p_common_chip = np.zeros((unique_clusters, 1))
    cdf_common_chip = np.zeros((unique_clusters, 1))
    sig_common_chip = np.zeros((unique_clusters, 1))
    bonf_chip = threshold / unique_clusters

    # compute hypergeometric test for common chip intersection
    for i, cl in enumerate(sorted(cluster_counts.index)):
        K = cluster_counts.loc[cl]
        k = obs_common_chip.get(cl, 0)
        p = hypergeom.sf(k-1, pop_size, K, sample_common_chip)
        cdf = hypergeom.cdf(k-1, pop_size, K, sample_common_chip)
        p_common_chip[i, 0] = p
        cdf_common_chip[i, 0] = cdf
        if p < cdf and p <= bonf_chip:
            sig_common_chip[i, 0] = 1
        elif cdf <= bonf_chip:
            sig_common_chip[i, 0] = -1
        else:
            sig_common_chip[i, 0] = 0

    # save common chip significance matrix
    cols_chip = ["common_chip"]
    pd.DataFrame(p_common_chip, index=sorted(cluster_counts.index), columns=cols_chip) \
      .to_csv(os.path.join(output_directory, "common_chip_p_values_df.csv"))
    pd.DataFrame(cdf_common_chip, index=sorted(cluster_counts.index), columns=cols_chip) \
      .to_csv(os.path.join(output_directory, "common_chip_cdf_values_df.csv"))
    pd.DataFrame(sig_common_chip, index=sorted(cluster_counts.index), columns=cols_chip) \
      .to_csv(os.path.join(output_directory, "common_chip_sig_matrix_df.csv"))

    # plot common chip heatmap
    plt.figure(figsize=(4, max(5, unique_clusters/2)))
    sns.heatmap(sig_common_chip, cmap='coolwarm', cbar=False,
                xticklabels=cols_chip, yticklabels=sorted(cluster_counts.index),
                linewidths=1, linecolor='black')
    plt.xlabel('Common ChIP')
    plt.ylabel('Cluster')
    plt.title('Common ChIP Intersection Significance')
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, "common_chip_sig_matrix.png"))
    plt.close()

    print("Analysis complete. Outputs written to:", output_directory)

if __name__ == "__main__":
    main()