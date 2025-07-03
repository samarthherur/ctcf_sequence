import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import subprocess 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import scipy.stats as stats

from matplotlib.colors import ListedColormap
import logomaker as lm

#####################################################################################################################
def fimo_to_neighbourhood(fimo_filepath, output_dir, genome_fasta, n_a = 50, n_b = 50):
    # Read in FIMO output
    fimo = pd.read_csv(fimo_filepath, sep="\t")
        
    # Remove rows with negative start values or NaNs
    fimo = fimo[fimo['start'] > 0].dropna(subset=['start', 'stop']).reset_index(drop=True)
    
    new_fimo = pd.DataFrame()
    
    # Iterate over each chromosome
    for chrom in fimo['sequence_name'].unique():
        # Sort and reset index
        fimo_chrom = (
            fimo[fimo['sequence_name'] == chrom]
            .sort_values('start')
            .reset_index(drop=True)
        )
        
        fimo_chrom_without_overlap = pd.DataFrame()
        last_stop = -1
        
        for i in range(len(fimo_chrom)):
            if fimo_chrom.iloc[i]['start'] > last_stop:
                # concat the single-row slice back onto fimo_chrom_without_overlap
                fimo_chrom_without_overlap = pd.concat(
                    [fimo_chrom_without_overlap, fimo_chrom.iloc[[i]]],
                    ignore_index=True
                )
                last_stop = fimo_chrom.iloc[i]['stop']
        
        new_fimo = pd.concat([new_fimo, fimo_chrom_without_overlap], ignore_index=True)
    
    fimo = new_fimo
    
    # Convert start/stop to int
    fimo['start'] = fimo['start'].astype(int)
    fimo['stop']  = fimo['stop'].astype(int)
    
    # Split by strand
    fimo_plus  = fimo[fimo['strand'] == '+'].copy()
    fimo_minus = fimo[fimo['strand'] == '-'].copy()
    
    # Shift minus-strand coordinates by -1
    fimo_minus['start'] -= 1
    fimo_minus['stop']  -= 1

    # Apply asymmetric extension per strand
    # Plus strand: upstream by n_a, downstream by n_b
    fimo_plus['start'] = fimo_plus['start'] - n_a
    fimo_plus['stop']  = fimo_plus['stop'] + n_b
    # Minus strand: upstream (genomic downstream) by n_a, downstream by n_b
    fimo_minus['start'] = fimo_minus['start'] - n_b
    fimo_minus['stop']  = fimo_minus['stop'] + n_a
    
    # Build 'name' column
    for df in (fimo_plus, fimo_minus):
        df['name'] = (
            df['sequence_name']
            + ':'
            + df['start'].astype(str)
            + '-'
            + df['stop'].astype(str)
        )
    
    # Columns needed for BED
    fimo_plus_bed  = fimo_plus[['sequence_name','start','stop','name','score','strand']]
    fimo_minus_bed = fimo_minus[['sequence_name','start','stop','name','score','strand']]
    
    # Write plus/minus BED files
    fimo_plus_bed.to_csv(os.path.join(output_dir, 'fimo_plus.bed'),
                          sep="\t", header=False, index=False)
    fimo_minus_bed.to_csv(os.path.join(output_dir, 'fimo_minus.bed'),
                           sep="\t", header=False, index=False)
    
    # Concatenate to get final neighbourhood BED
    fimo_neighbourhood_bed = pd.concat([fimo_plus_bed, fimo_minus_bed], ignore_index=True)
    fimo_neighbourhood_bed.to_csv(
        os.path.join(output_dir, f'fimo_neighbourhood_{n_a}_{n_b}.bed'),
        sep="\t", header=False, index=False
    )
    
    # Use bedtools to fetch FASTA for both strands
    os.system(
        f'bedtools getfasta -s -fi {genome_fasta} '
        f'-bed {os.path.join(output_dir,"fimo_minus.bed")} '
        f'-fo {os.path.join(output_dir,"fimo_minus.fasta")}'
    )
    os.system(
        f'bedtools getfasta -s -fi {genome_fasta} '
        f'-bed {os.path.join(output_dir,"fimo_plus.bed")} '
        f'-fo {os.path.join(output_dir,"fimo_plus.fasta")}'
    )
    
    # Combine FASTA records
    fimo_minus_seqs = list(SeqIO.parse(os.path.join(output_dir, 'fimo_minus.fasta'), 'fasta'))
    fimo_plus_seqs  = list(SeqIO.parse(os.path.join(output_dir, 'fimo_plus.fasta'),  'fasta'))
    fimo_neighbourhood_seqs = fimo_minus_seqs + fimo_plus_seqs
    
    SeqIO.write(
        fimo_neighbourhood_seqs,
        os.path.join(output_dir, f'fimo_neighbourhood_{n_a}_{n_b}.fasta'),
        'fasta'
    )
    
    # Print summary
    print('FIMO neighbourhood files created:')
    print('  fimo_plus.bed  - neighbourhood file req for FASTA output')
    print('  fimo_minus.bed - neighbourhood file req for FASTA output')
    print(f'  fimo_neighbourhood_{n_a}_{n_b}.bed')
    print(f'  fimo_neighbourhood_{n_a}_{n_b}.fasta')
    
    
#####################################################################################################################

###function to plot fasta into heatmap


def plot_fasta_heatmap(fasta_file):
    #parse the fasta file
    list_seq = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        seq = record.seq
        list_seq.append(str(seq))
    
    #check if all sequences are of the same length
    seqlength = len(list_seq[0])
    for seq in list_seq:
        if len(seq) != seqlength:
            raise ValueError('Sequences are not of the same length')
        else:
            continue
   
    
    #nucleotide dictionary
    nucleotide_dict = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3, 'N':4, 'n':4}
    
    #convert list of sequences to numpy array
    numberofseq = len(list_seq) 
    seqlength = len(list_seq[0])
    
    seq_array = np.zeros((numberofseq, seqlength))
    for i, seq in enumerate(list_seq):
        for j, nucleotide in enumerate(seq):
            seq_array[i,j] = nucleotide_dict[nucleotide]
    
    ##plot the heatmap
    colors = ['green', 'blue', 'yellow', 'red' , 'black']
    cmap = ListedColormap(colors)    
    norm = plt.Normalize(0,4)
    
    #plot the matrix
    fig, ax = plt.subplots(figsize=(20,10))
    cax = ax.matshow(seq_array, cmap = cmap, norm = norm, aspect='auto')
    #show the color bar
    fig.colorbar(cax)
    #set x and y axis labels
    plt.xlabel('Sequence length')
    plt.ylabel('Sequence number')
    #title
    plt.title('Neighbourhood Heatmap')
    # Save heatmap to output directory
    output_dir = os.path.dirname(fasta_file)
    heatmap_fname = os.path.join(output_dir, f"{os.path.basename(fasta_file).rsplit('.',1)[0]}_heatmap.png")
    fig.savefig(heatmap_fname)
    plt.close(fig)
    # plt.show()
    
#####################################################################################################################

# ###function to plot sequence logo from fasta file
def plot_logo_from_fasta(fasta_input):
    #parse the fasta file
    list_seq = []
    for record in SeqIO.parse(fasta_input, 'fasta'):
        seq = record.seq
        seq1 = str(seq)
        seq1 = seq1.upper()
        list_seq.append(seq1)
    
    #check if all sequences are of the same length
    seqlength = len(list_seq[0])
    for seq in list_seq:
        if len(seq) != seqlength:
            raise ValueError('Sequences are not of the same length')
        else:
            continue
    
    #create a dataframe
    df = lm.alignment_to_matrix(list_seq)
    #transform the matrix
    df = lm.transform_matrix(df, from_type='counts', to_type='information')
    #remove column 'N'
    if 'N' in df.columns:
        df = df.drop(columns = ['N'])    
    #plot the logo
    lm.Logo(df)
    # Save sequence logo to output directory
    output_dir = os.path.dirname(fasta_input)
    logo_fname = os.path.join(output_dir, f"{os.path.basename(fasta_input).rsplit('.',1)[0]}_logo.png")
    plt.gcf().savefig(logo_fname)
    plt.close()
    # plt.show()

#####################################################################################################################

### function to plot sequence logo from motif file
def chip_filtered_neighbourhood(neighbourhood_bed_file, chipseq_bed_file, genome_fasta):
    #read in the bed files
    neighbourhood_bed = pd.read_csv(neighbourhood_bed_file, sep="\t", header=None)
    chipseq_bed = pd.read_csv(chipseq_bed_file, sep="\t", header=None)
    
    #create a list of neighbourhoods that overlap with chipseq peaks
    os.system('bedtools intersect -a ' + neighbourhood_bed_file + ' -b ' + chipseq_bed_file + ' -u > ' + neighbourhood_bed_file[:-4] + '_filtered.bed')
    #get fasta file for the filtered neighbourhood + and - strands
    os.system('bedtools getfasta -s -fi ' + genome_fasta + ' -bed ' + neighbourhood_bed_file[:-4] + '_filtered.bed -fo ' + neighbourhood_bed_file[:-4] + '_filtered.fasta')
    
    #print number of sequences in the filtered neighbourhood and the original neighbourhood
    print(f'Number of sequences in the original neighbourhood: {len(neighbourhood_bed)}')
    filtered_neighbourhood_bed  = pd.read_csv(neighbourhood_bed_file[:-4] + '_filtered.bed', sep="\t", header=None)
    print(f'Number of sequences in the filtered neighbourhood: {len(filtered_neighbourhood_bed)}')
    
    #print the saved bed and fasta files
    print(f'Filtered neighbourhood bed file saved as {neighbourhood_bed_file[:-4]}_filtered.bed')
    print(f'Filtered neighbourhood fasta file saved as {neighbourhood_bed_file[:-4]}_filtered.fasta')
    
#####################################################################################################################


###function to intersect filtered neighbourhood with boundary bed file
def intersect_filtered_neighbourhood_boundary(filtered_neighbourhood_bed_file, boundary_bed_file, genome_fasta):
    #read in the bed files
    filtered_neighbourhood_bed = pd.read_csv(filtered_neighbourhood_bed_file, sep="\t", header=None)
    boundary_bed = pd.read_csv(boundary_bed_file, sep="\t", header=None)
    
    #create a list of neighbourhoods that overlap with boundaries
    os.system('bedtools intersect -a ' + filtered_neighbourhood_bed_file + ' -b ' + boundary_bed_file + ' -u > ' + filtered_neighbourhood_bed_file[:-4] + '_int_boundary.bed')
    #get fasta file
    os.system('bedtools getfasta -s -fi ' + genome_fasta + ' -bed ' + filtered_neighbourhood_bed_file[:-4] + '_int_boundary.bed -fo ' + filtered_neighbourhood_bed_file[:-4] + '_int_boundary.fasta')
    
    #report entries in a that have no overlap with b
    os.system('bedtools intersect -v -a ' + filtered_neighbourhood_bed_file + ' -b ' + boundary_bed_file + ' > ' + filtered_neighbourhood_bed_file[:-4] + '_no_int_boundary.bed')
    #get fasta file
    os.system('bedtools getfasta -s -fi ' + genome_fasta + ' -bed ' + filtered_neighbourhood_bed_file[:-4] + '_no_int_boundary.bed -fo ' + filtered_neighbourhood_bed_file[:-4] + '_no_int_boundary.fasta')
    
    
    #print number of sequences in the filtered neighbourhood and the original neighbourhood
    print(f'Number of sequences in the filtered neighbourhood: {len(filtered_neighbourhood_bed)}')
    int_boundary_bed  = pd.read_csv(filtered_neighbourhood_bed_file[:-4] + '_int_boundary.bed', sep="\t", header=None)
    print(f'Number of neighbourhood sequences that intersect with boundaries: {len(int_boundary_bed)}')
    no_int_boundary_bed  = pd.read_csv(filtered_neighbourhood_bed_file[:-4] + '_no_int_boundary.bed', sep="\t", header=None)
    print(f'Number of neighbourhood sequences that do not intersect with boundaries: {len(no_int_boundary_bed)}')
    
    #print the saved bed and fasta files
    print(f'Neighbourhood intersecting with boundary bed file saved as {filtered_neighbourhood_bed_file[:-4]}_int_boundary.bed')
    print(f'Neighbourhood that does not intersect with boundary bed file saved as {filtered_neighbourhood_bed_file[:-4]}_no_int_boundary.bed')
    print(f'Neighbourhood intersecting with boundary fasta file saved as {filtered_neighbourhood_bed_file[:-4]}_int_boundary.fasta')
    print(f'Neighbourhood that does not intersect with boundary fasta file saved as {filtered_neighbourhood_bed_file[:-4]}_no_int_boundary.fasta')


#####################################################################################################################

###function that creates a matrix of all neighbourhoods
def matrix_neighbourhood(neighbourhood_bed_file, chipseq_bed, boundary_bed, output_dir):
    import os  # Ensure os is imported inside the function
    # Define output paths using os.path.join
    chip_out = os.path.join(output_dir, 'neighbourhood_ifchip.bed')
    boundary_out = os.path.join(output_dir, 'neighbourhood_ifboundary.bed')

    #create a list of neighbourhoods that overlap with chipseq peaks 
    #bedtools with -c option will report the number of times that the feature in A overlaps with features in B
    os.system(f'bedtools intersect -a {neighbourhood_bed_file} -b {chipseq_bed} -c > {chip_out}')

    #create a list of neighbourhoods that overlap with boundaries
    os.system(f'bedtools intersect -a {neighbourhood_bed_file} -b {boundary_bed} -c > {boundary_out}')

    #read in the bed files
    neighbourhood_ifchip = pd.read_csv(chip_out, sep="\t", header=None)
    neighbourhood_ifboundary = pd.read_csv(boundary_out, sep="\t", header=None)

    if len(neighbourhood_ifboundary) != len(neighbourhood_ifchip):
        raise ValueError('Length of neighbourhood_ifchip and neighbourhood_ifboundary are not equal')
    else:
        pass

    #apply step function both neighbourhood_chip and neighbourhood_boundary
    neighbourhood_ifchip[6] = np.where(neighbourhood_ifchip[6] > 0, 1, 0)
    neighbourhood_ifboundary[6] = np.where(neighbourhood_ifboundary[6] > 0, 1, 0)

    #create a df with 6th column of neighbourhood_chip and neighbourhood_boundary
    df = pd.DataFrame()
    df['chip'] = neighbourhood_ifchip[6]
    df['boundary'] = neighbourhood_ifboundary[6]
    df['score'] = neighbourhood_ifchip[4]

    df.head()

    # save the df to tsv file
    matrix_path = os.path.join(output_dir, 'neighbourhood_matrix.tsv')
    df.to_csv(matrix_path, sep="\t", index=False)
    print(f'Neighbourhood matrix saved as {matrix_path}')
    
def fasta_without_repeats(fasta_path):
    """
    Read a FASTA, filter out sequences containing N/n or lowercase letters,
    and return the filtered SeqRecord list and count of retained sequences.
    """
    from Bio import SeqIO
    seqs = []
    total = 0
    for rec in SeqIO.parse(fasta_path, "fasta"):
        total += 1
        s = str(rec.seq)
        if "N" in s or "n" in s or not s.isupper():
            continue
        seqs.append(rec)
    return seqs, len(seqs)
    

# -------------------------------------------------------------------------------------------------
# Plot intersect vs. non-intersect score distributions
def plot_intersect_distribution(int_scores, no_int_scores, output_dir, n_a, n_b, bins=30):
    """
    Plot & save histogram + KDE comparing intersect vs non-intersect scores.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os

    # define bins across combined range
    min_score = min(int_scores.min(), no_int_scores.min())
    max_score = max(int_scores.max(), no_int_scores.max())
    bin_edges = np.linspace(min_score, max_score, bins)

    colors = {"Intersect": "tab:blue", "No Intersect": "tab:orange"}
    plt.figure(figsize=(8, 6))
    sns.histplot(int_scores, bins=bin_edges, stat="density", element="step",
                 fill=True, color=colors["Intersect"], alpha=0.3, label="Intersect")
    sns.histplot(no_int_scores, bins=bin_edges, stat="density", element="step",
                 fill=True, color=colors["No Intersect"], alpha=0.3, label="No Intersect")
    sns.kdeplot(int_scores, bw_adjust=1, linewidth=2, color=colors["Intersect"], label="Intersect KDE")
    sns.kdeplot(no_int_scores, bw_adjust=1, linewidth=2, color=colors["No Intersect"], label="No Intersect KDE")
    plt.xlabel("Scores")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    outfile = os.path.join(output_dir, f"score_distribution_{n_a}_{n_b}.png")
    plt.savefig(outfile)
    plt.close()


# -------------------------------------------------------------------------------------------------
# Count windows for various neighbourhood stages and repeat-free sequences
def count_neighbourhood_windows(output_dir, n_a, n_b):
    """
    Count windows for raw, filtered, intersecting, non-intersecting stages,
    plus repeat-free sequences, writing counts to report.
    """
    import os
    from Bio import SeqIO

    stages = {
        'raw': f"fimo_neighbourhood_{n_a}_{n_b}.bed",
        'filtered': f"fimo_neighbourhood_{n_a}_{n_b}_filtered.bed",
        'int_boundary': f"fimo_neighbourhood_{n_a}_{n_b}_filtered_int_boundary.bed",
        'no_int_boundary': f"fimo_neighbourhood_{n_a}_{n_b}_filtered_no_int_boundary.bed"
    }
    for name, fname in stages.items():
        path = os.path.join(output_dir, fname)
        if os.path.exists(path):
            with open(path) as bf:
                cnt = sum(1 for _ in bf)
            print(f"Number of {name} neighbourhood windows: {cnt}")

    wr_fasta = os.path.join(output_dir, f"fimo_neighbourhood_{n_a}_{n_b}_filtered_wr.fasta")
    if os.path.exists(wr_fasta):
        cnt = sum(1 for _ in SeqIO.parse(wr_fasta, "fasta"))
        print(f"Number of repeat-free sequences: {cnt}")
    
