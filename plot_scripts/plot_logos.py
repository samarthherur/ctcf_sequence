import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logomaker
import multiprocessing as mp

# ---------------
# Helper Functions
# ---------------

def calculate_pwm_df(sequences):
    """
    Compute a PWM (Position Weight Matrix) from a list of sequences.
    Each sequence must be of equal length.
    """
    sequence_length = len(sequences[0])
    nucleotides = ['A', 'C', 'G', 'T']
    num_sequences = len(sequences)
    sequence_array = np.array([list(seq) for seq in sequences])
    counts = np.zeros((sequence_length, len(nucleotides)), dtype=int)
    for i, nuc in enumerate(nucleotides):
        counts[:, i] = np.sum(sequence_array == nuc, axis=0)
    probabilities = counts / num_sequences
    pwm_df = pd.DataFrame(probabilities, columns=nucleotides)
    return pwm_df

def calculate_pwm_for_cluster(args):
    """
    Helper for multiprocessing: computes the PWM for a given cluster.
    """
    cluster, sequences = args
    pwm_df = calculate_pwm_df(sequences)
    return cluster, pwm_df

def read_logo_file(file_path):
    """
    Read the logo file (architecture details) and return a DataFrame with columns ['Cluster','Sequence'].
    Assumes no header and that the first column is cluster ID and the third column is sequence.
    """
    data = pd.read_csv(file_path, sep='\t', header=None)
    data = data.iloc[:, [0, 2]]
    data.columns = ['Cluster', 'Sequence']
    data['Cluster'] = data['Cluster'].astype(int)
    return data

def read_binding_file(file_path):
    """
    Read the binding file and extract 'cluster' and 'binding' columns.
    Assumes cluster in column 4 and binding in column 6.
    """
    df = pd.read_csv(file_path, sep='\t', header=None)
    binding_df = pd.DataFrame({
        'cluster': df.iloc[:, 3].astype(int),
        'binding': pd.to_numeric(df.iloc[:, 5], errors='coerce')
    })
    return binding_df

def read_tss_file(file_path):
    """
    Read TSS file, extract cluster (col4), strand (col5), distance (col12),
    invert distance for negative strand, and return abs_distance.
    """
    tss_df = pd.read_csv(file_path, sep='\t', header=None, usecols=[3,4,11])
    tss_df.columns = ['cluster', 'strand', 'distance']
    tss_df['cluster'] = tss_df['cluster'].astype(int)
    tss_df['distance'] = pd.to_numeric(tss_df['distance'], errors='coerce')
    tss_df['abs_distance'] = tss_df.apply(
        lambda row: abs(row['distance']) if row['strand'] == '+' else abs(-row['distance']), axis=1
    )
    return tss_df[['cluster', 'abs_distance']]

# ---------------
# Plotting Function
# ---------------

def plot_combined_figure(logo_pwm_dict, binding_df, tss_df, out_path, cluster_order):
    """
    Create a figure with one row per cluster and three columns:
      1) Sequence logo (PWM)
      2) Horizontal bar plot of average binding
      3) Box or scatter plot of absolute TSS distance
    """
    num_clusters = len(cluster_order)
    fig, axs = plt.subplots(
        nrows=num_clusters, ncols=3,
        figsize=(15, 1.5 * num_clusters),
        constrained_layout=True
    )
    if num_clusters == 1:
        axs = np.expand_dims(axs, axis=0)

    for i, cluster in enumerate(cluster_order):
        # Logo
        ax0 = axs[i, 0]
        pwm = logo_pwm_dict.get(cluster)
        if pwm is not None:
            info = logomaker.transform_matrix(pwm, from_type='probability', to_type='information')
            logomaker.Logo(info, ax=ax0)
        ax0.set_title(f"Cluster {cluster}")

        # Binding bar
        ax1 = axs[i, 1]
        subb = binding_df[binding_df['cluster'] == cluster]
        if not subb.empty:
            ax1.barh([0], [subb['binding'].mean()], color='skyblue')
        ax1.set_xlim(0, 3)
        ax1.set_yticks([])
        ax1.set_title("Avg Binding")

        # TSS distance
        ax2 = axs[i, 2]
        subt = tss_df[tss_df['cluster'] == cluster]
        if len(subt) > 1:
            sns.boxplot(x='abs_distance', data=subt, ax=ax2, orient='h', showfliers=False)
        elif len(subt) == 1:
            ax2.scatter(subt['abs_distance'], [0], color='b')
        ax2.set_xlim(0, 240000)
        ax2.set_yticks([])
        ax2.set_title("TSS Distance")

    plt.savefig(out_path)
    plt.close()

# ---------------
# Main Function
# ---------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot logos, binding, and TSS for NPLB clusters"
    )
    parser.add_argument('-l', '--logo_file',    required=True, help='Path to architectureDetails file')
    parser.add_argument('-b', '--binding_file', required=True, help='Path to binding counts file')
    parser.add_argument('-t', '--tss_file',     required=True, help='Path to closest TSS file')
    parser.add_argument('-o', '--out_path',     required=True, help='Output PNG path')
    args = parser.parse_args()

    # Read data
    logo_data = read_logo_file(args.logo_file)
    cluster_groups = logo_data.groupby('Cluster')['Sequence'].apply(list).to_dict()
    items = list(cluster_groups.items())
    with mp.Pool() as pool:
        pwm_results = pool.map(calculate_pwm_for_cluster, items)
    logo_pwm_dict = dict(pwm_results)

    binding_df = read_binding_file(args.binding_file)
    # Optionally save binding_df to TSV next to output
    binding_output = os.path.join(os.path.dirname(args.out_path), 'binding_df.tsv')
    binding_df.to_csv(binding_output, sep='\t', index=False)

    tss_df = read_tss_file(args.tss_file)
    tss_output = os.path.join(os.path.dirname(args.out_path), 'tss_df.tsv')
    tss_df.to_csv(tss_output, sep='\t', index=False)

    cluster_order = sorted(logo_data['Cluster'].unique())
    plot_combined_figure(logo_pwm_dict, binding_df, tss_df, args.out_path, cluster_order)

if __name__ == '__main__':
    main()