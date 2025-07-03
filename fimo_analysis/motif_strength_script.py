import bisect
from collections import defaultdict

# Function to read BED file and return a sorted list of intervals
def read_bed(file_path):
    intervals = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip() == '':
                continue
            cols = line.strip().split('\t')
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            intervals.append((chrom, start, end))
    # Sort intervals by chromosome and start position
    intervals.sort()
    return intervals

# Function to read motifs from FIMO TSV output
def read_fimo(file_path):
    motifs = []
    with open(file_path, 'r') as f:
        i = 0
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue  # Skip lines that don't have enough columns
            if i == 0:
                i += 1
                continue
            
            motif_id = cols[0]
            motif_alt_id = cols[1]
            sequence_name = cols[2]
            start = int(cols[3]) - 1  # Convert to 0-based coordinate
            end = int(cols[4])        # FIMO uses 1-based inclusive coordinates
            strand = cols[5]
            score = cols[6]
            p_value = cols[7]
            q_value = cols[8]
            matched_seq = cols[9] if len(cols) > 9 else ''
            motifs.append({
                'motif_id': motif_id,
                'motif_alt_id': motif_alt_id,
                'chrom': sequence_name,
                'start': start,
                'end': end,
                'strand': strand,
                'score': score,
                'p_value': p_value,
                'q_value': q_value,
                'matched_seq': matched_seq,
                'boundary_index': -1,  # To be assigned later
                'classification': ''
            })
    return motifs

# Assign motifs to boundaries
def assign_motifs_to_boundaries(boundaries, motifs):
    # Create an index for boundaries for quick lookup
    boundary_dict = defaultdict(list)
    for idx, (chrom, start, end) in enumerate(boundaries):
        boundary_dict[chrom].append((start, end, idx))
    for chrom in boundary_dict:
        boundary_dict[chrom].sort()
    # Assign motifs to boundaries
    for motif in motifs:
        chrom = motif['chrom']
        if chrom not in boundary_dict:
            motif['classification'] = 'Class 1'
            continue
        positions = boundary_dict[chrom]
        # Binary search to find the correct interval
        idx = bisect.bisect_left(positions, (motif['start'], motif['end'], 0))
        found = False
        for i in [idx - 1, idx]:
            if 0 <= i < len(positions):
                start, end, boundary_idx = positions[i]
                if motif['start'] >= start and motif['end'] <= end:
                    motif['boundary_index'] = boundary_idx
                    found = True
                    break
        if not found:
            motif['classification'] = 'Class 1'
    return motifs


def classify_motifs(boundaries, motifs):
    """
    Classifies motifs into three classes:
    1) Class 1: Motifs outside any boundary.
    2) Class 3: '+' strand motifs that have '-' strand motifs in the next boundary.
       The '-' strand motifs in the next boundary are also classified as Class 3.
    3) Class 2: All other motifs within boundaries.
    """
    from collections import defaultdict

    # Group motifs by boundary index
    motifs_by_boundary = defaultdict(list)
    for motif in motifs:
        boundary_idx = motif['boundary_index']
        motifs_by_boundary[boundary_idx].append(motif)

    # Get the sorted list of boundary indices, excluding -1 (outside boundaries)
    boundary_indices = sorted(idx for idx in motifs_by_boundary.keys() if idx != -1)

    # Map boundary index to its position in the sorted order
    boundary_position = {idx: pos for pos, idx in enumerate(boundary_indices)}

    # For each boundary, process '+' strand motifs
    for current_boundary_idx in boundary_indices:
        current_motifs = motifs_by_boundary[current_boundary_idx]
        plus_strand_motifs = [motif for motif in current_motifs if motif['strand'] == '+']

        # Get next boundary
        current_pos = boundary_position[current_boundary_idx]
        next_pos = current_pos + 1
        if next_pos < len(boundary_indices):
            next_boundary_idx = boundary_indices[next_pos]
            next_motifs = motifs_by_boundary[next_boundary_idx]
            minus_strand_motifs = [motif for motif in next_motifs if motif['strand'] == '-']

            if plus_strand_motifs and minus_strand_motifs:
                # Classify '+' strand motifs in current boundary as Class 3
                for motif in plus_strand_motifs:
                    motif['classification'] = 'Class 3'
                # Classify '-' strand motifs in next boundary as Class 3
                for motif in minus_strand_motifs:
                    motif['classification'] = 'Class 3'

    # After processing, set any unclassified motifs within boundaries to Class 2
    for motif in motifs:
        if motif['classification'] == '':
            if motif['boundary_index'] != -1:
                motif['classification'] = 'Class 2'
            else:
                # Motifs outside boundaries are already classified as Class 1
                pass

    return motifs



# Main function
def main(boundary_file, motif_file, output_file):
    boundaries = read_bed(boundary_file)
    motifs = read_fimo(motif_file)
    motifs = assign_motifs_to_boundaries(boundaries, motifs)
    motifs = classify_motifs(boundaries, motifs)
    # Write output
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(['chrom', 'start', 'end', 'motif_id', 'strand', 'score', 'classification']) + '\n')
        for motif in motifs:
            if motif['classification'] == '':
                motif['classification'] = 'Class 2'  # Default to Class 2 if not assigned
            f.write('\t'.join([
                motif['chrom'],
                str(motif['start']),
                str(motif['end']),
                motif['motif_id'],
                motif['strand'],
                motif['score'],
                motif['classification']
            ]) + '\n')
            
##############################################################################################################
            
import csv
import os

def filter_fimo_by_ctcf_peaks(fimo_tsv, ctcf_peaks_bed, output_tsv):
    """
    Filters the FIMO output TSV file to include only motifs that overlap with CTCF peaks.
    """
    # Step 1: Convert fimo.tsv to BED format with unique IDs
    fimo_bed_with_id = 'fimo_with_id.bed'
    with open(fimo_tsv, 'r') as tsv_file, open(fimo_bed_with_id, 'w', newline='') as bed_file:
        reader = csv.DictReader(tsv_file, delimiter='\t')
        for idx, row in enumerate(reader):
            # Skip rows with missing 'start' or 'stop'
            if not row.get('start') or not row.get('stop'):
                print(f"Skipping row with missing 'start' or 'stop': {row}")
                continue
            try:
                chrom = row['sequence_name']
                start = int(row['start']) - 1  # Convert to 0-based BED coordinate
                end = int(row['stop'])         # End remains the same
                name = row['motif_id']
                score = row['score']
                strand = row['strand']
                unique_id = str(idx)
                # BED format: chrom, start, end, name, score, strand, unique_id
                bed_fields = [chrom, str(start), str(end), name, score, strand, unique_id]
                bed_file.write('\t'.join(bed_fields) + '\n')
            except ValueError as ve:
                print(f"ValueError encountered: {ve} in row: {row}")
                continue

    # Step 2: Intersect with CTCF peaks using bedtools
    intersected_bed = 'fimo_intersected_with_id.bed'
    intersect_command = f'bedtools intersect -a {fimo_bed_with_id} -b {ctcf_peaks_bed} -wa > {intersected_bed}'
    os.system(intersect_command)

    # Step 3: Extract corresponding lines from fimo.tsv
    # Read IDs from intersected BED file
    intersected_ids = set()
    with open(intersected_bed, 'r') as bed_file:
        for line in bed_file:
            cols = line.strip().split('\t')
            unique_id = cols[6]  # unique_id is the 7th column (0-based index 6)
            intersected_ids.add(unique_id)

    # Read fimo.tsv and write entries with matching IDs to output_tsv
    with open(fimo_tsv, 'r') as tsv_file, open(output_tsv, 'w', newline='') as out_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        writer = csv.writer(out_file, delimiter='\t')
        header = next(reader)
        writer.writerow(header)  # Write header to output
        for idx, row in enumerate(reader):
            if str(idx) in intersected_ids:
                writer.writerow(row)

    # Cleanup intermediate files
    os.remove(fimo_bed_with_id)
    os.remove(intersected_bed)
    
##############################################################################################################

##############################################################################################################
# -------------------------------------------------------------------------------------------------
def summarize_classified_motifs(classified_tsv, output_dir):
    """
    Summarize classified motifs: counts, medians, KS tests, and density plot.
    """
    import pandas as pd
    import numpy as np
    import scipy.stats as stats
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    from itertools import combinations

    df = pd.read_csv(classified_tsv, sep='\t')
    total = len(df)
    print(f"Total classified motifs: {total}")

    classes = np.sort(df['classification'].unique())
    for cls in classes:
        cnt = (df['classification'] == cls).sum()
        print(f"Count {cls}: {cnt}")

    for cls in classes:
        med = df.loc[df['classification']==cls, 'score'].median()
        print(f"Median score {cls}: {med}")

    for c1, c2 in combinations(classes, 2):
        stat, pval = stats.ks_2samp(
            df.loc[df['classification']==c1, 'score'],
            df.loc[df['classification']==c2, 'score']
        )
        print(f"KS {c1} vs {c2}: stat={stat}, p={pval}")

    plt.figure(figsize=(8, 6))
    for cls in classes:
        sns.kdeplot(df.loc[df['classification']==cls, 'score'], shade=True, label=cls)
    plt.xlabel('Scores')
    plt.ylabel('Density')
    plt.legend()
    plt.tight_layout()
    outfile = os.path.join(output_dir, 'classification_density.png')
    plt.savefig(outfile)
    plt.close()
    print(f"Density plot saved to: {outfile}")

# Example usage
if __name__ == '__main__':
    boundary_file = '/home/samarth/dekker_microc/H1ESC/4DNFIED5HLDC.bed'
    motif_file = '/home/samarth/hg38/fimo_ctcf_hg38_new/fimo.tsv'
    output_file = '/home/samarth/hg38/fimo_ctcf_hg38_new/classified_motifs.tsv'
    chip_file = '/home/samarth/ChIP_Seq/chip_h1_hg38_narrowpeak.bed'
    
    filtered_fimo = motif_file.replace('.tsv', '_filtered.tsv')
    filter_fimo_by_ctcf_peaks(motif_file, chip_file, filtered_fimo)
    main(boundary_file, filtered_fimo, output_file)
