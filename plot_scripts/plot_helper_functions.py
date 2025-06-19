import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os 
from scipy.stats import wilcoxon
from scipy.stats import ttest_ind
import sys


def plot_score(output_dir, column_number, cluster_values, df_names, p_value_thr, neighbourhood_df, ylim=None):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize DataFrames for each cluster
    new_dfs = []
    
    # Calculate average scores and store results
    score_dict = {}
    
    for i, cluster_value in enumerate(cluster_values):
        new_df = pd.DataFrame()
        
        for j in range(df_names[i].shape[0]):
            chr_val = df_names[i].iloc[j, 0]
            start = df_names[i].iloc[j, 1]
            end = df_names[i].iloc[j, 2]
            
            # Filter neighbourhood data
            filtered_df = neighbourhood_df.loc[
                (neighbourhood_df[0] == chr_val) & 
                (neighbourhood_df[1] == int(start)) & 
                (neighbourhood_df[2] == int(end))
            ]
            new_df = pd.concat([new_df, filtered_df])
        
        new_dfs.append(new_df)
        
        # Calculate and log average score for each cluster
        average_score = new_df[column_number].mean()
        score_dict[cluster_value] = average_score
        print(f'Average score for cluster {cluster_value}: {round(average_score, 4)}')
        
        with open(f'{output_dir}/average_score_{column_number}.txt', 'a') as f:
            f.write(f'Average score for cluster {cluster_value}: {round(average_score, 4)}\n')
    
    
    # Concatenate all DataFrames for overall statistics
    concatenated_df = pd.concat(new_dfs)
    average_score_overall = concatenated_df[column_number].mean()
    print(f'Average score overall: {round(average_score_overall, 4)}')
    
    with open(f'{output_dir}/average_score_{column_number}.txt', 'a') as f:
        f.write(f'Average score overall: {round(average_score_overall, 4)}\n')
    
    # Wilcoxon test and significance detection
    sig_p_value_dict = {}
    shuffled_df = concatenated_df[column_number].sample(frac=1).reset_index(drop=True)
    
    for i, cluster_value in enumerate(cluster_values):
        stat, p_value = ttest_ind(new_dfs[i][column_number], shuffled_df)
        print(f't-test for cluster {cluster_value}: {round(p_value, 4)}')
        
        with open(f'{output_dir}/average_score_{column_number}.txt', 'a') as f:
            f.write(f't-test for cluster {cluster_value}: {round(p_value, 4)}\n')
        
        sig_p_value_dict[cluster_value] = p_value < p_value_thr
    
    # Plotting
    fig, axs = plt.subplots(2, 2, figsize=(20, 16))
    
    # Bar plot with significance highlighting
    for i, cluster_value in enumerate(cluster_values):
        color = 'red' if sig_p_value_dict[cluster_value] else 'blue'
        label = 'Significant' if sig_p_value_dict[cluster_value] else 'Not significant'
        axs[0, 0].bar(cluster_value, score_dict[cluster_value], color=color, label=label)
    
    
    axs[0, 0].set_title('Average score per cluster')
    axs[0, 0].set_xlabel('Cluster')
    axs[0, 0].set_ylabel('Average score')
    axs[0, 0].set_xticks(range(1, len(cluster_values) + 1))
    axs[0, 0].set_xticklabels(cluster_values)
    axs[0, 0].axhline(y=average_score_overall, color='black', linestyle='--')
    if ylim is not None:
        axs[0, 0].set_ylim(ylim[0], ylim[1])
    
    # Boxplot of scores
    axs[0, 1].boxplot([df[column_number] for df in new_dfs])
    axs[0, 1].set_title('Boxplot of Scores')
    axs[0, 1].set_xlabel('Cluster')
    axs[0, 1].set_ylabel('Score')
    axs[0, 1].set_xticks(range(1, len(cluster_values) + 1))
    axs[0, 1].set_xticklabels(cluster_values)
    
    # Violin plot of scores
    sns.violinplot(data=[df[column_number] for df in new_dfs], ax=axs[1, 0])
    axs[1, 0].set_title('Violin plot of Scores')
    axs[1, 0].set_xlabel('Cluster')
    axs[1, 0].set_ylabel('Score')
    axs[1, 0].set_xticks(range(len(cluster_values)))
    axs[1, 0].set_xticklabels(cluster_values)
    
    # Density plot of scores
    for df in new_dfs:
        sns.kdeplot(df[column_number], ax=axs[1, 1])
    axs[1, 1].set_title('Density plot of Scores')
    axs[1, 1].set_xlabel('Score')
    axs[1, 1].set_ylabel('Density')
    axs[1, 1].legend(cluster_values)
    
    plt.savefig(f'{output_dir}/score_plot_{column_number}.png')
    plt.show()  
    
    ##print order of clusters by average chip-seq score
    list_of_clusters = []
    for i in range(len(cluster_values)):
        list_of_clusters.append(cluster_values[i])
    list_of_clusters.sort(key=score_dict.get, reverse=True)
    print(f'Order of clusters by average chip-seq score: {list_of_clusters}')


#parse the architectural details file to get the chr, start and end coordinates
def parse_architectural_details(architectural_details_file):
    
    parse_architectural_details = pd.DataFrame(columns=['chr', 'start', 'end', 'cluster'])
    for i in range(architectural_details_file.shape[0]):
        chr = architectural_details_file.iloc[i, 1].split(':')[0]
        start = architectural_details_file.iloc[i, 1].split(':')[1].split('-')[0]
        
        try:
            end = architectural_details_file.iloc[i, 1].split(':')[1].split('-')[2]
        except IndexError:
            end = architectural_details_file.iloc[i, 1].split(':')[1].split('-')[1]
            end = end[:-3]
        else:
            end = architectural_details_file.iloc[i, 1].split(':')[1].split('-')[1]
            end = end[:-1] 
            
        cluster = architectural_details_file.iloc[i, 0]
        parse_architectural_details.loc[i] = [chr, start, end, cluster]
        
    return parse_architectural_details
              
        

    
    
    
    
    
    
    
    
    
    
    
    
    
        
        
        
    
    
        
    
    

    
    
    
    
    
    
    
    
    
        
        
    

    
    
    
    