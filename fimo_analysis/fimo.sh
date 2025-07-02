#!/bin/bash
#SBATCH --job-name=fimo
#SBATCH --output=fimo.out
#SBATCH --error=fimo.err
#SBATCH --time=10:00:00
#SBATCH --partition=gpu
#SBATCH --ntasks=6
#SBATCH --nodelist=cn1

python /home/samarth/ctcf_sequence/fimo_analysis/fimo_neighbourhood_analysis.py \
  --fimo_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/fimo.tsv \
  --genome_fasta    /home/samarth/ctcf_sequence_data/input_data/genome_files/hg38.fa \
  --chip_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/chip_h1_hg38_narrowpeak.bed \
  --boundary        /home/samarth/ctcf_sequence_data/input_data/h1_human_data/4DNFIED5HLDC_sorted.bed \
  --output_dir      /home/samarth/ctcf_sequence_data/output/fimo/H1esc_CTCF_50_50 \
  --n_a 50 \
  --n_b 50 \
  --ignore_repeats
    
python /home/samarth/ctcf_sequence/fimo_analysis/fimo_neighbourhood_analysis.py \
  --fimo_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/fimo.tsv \
  --genome_fasta    /home/samarth/ctcf_sequence_data/input_data/genome_files/hg38.fa \
  --chip_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/chip_h1_hg38_narrowpeak.bed \
  --boundary        /home/samarth/ctcf_sequence_data/input_data/h1_human_data/4DNFIED5HLDC_sorted.bed \
  --output_dir      /home/samarth/ctcf_sequence_data/output/fimo/HFF_CTCF_50_50 \
  --n_a 50 \
  --n_b 50 \
  --ignore_repeats

python /home/samarth/ctcf_sequence/fimo_analysis/fimo_neighbourhood_analysis.py \
  --fimo_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/fimo.tsv \
  --genome_fasta    /home/samarth/ctcf_sequence_data/input_data/genome_files/hg38.fa \
  --chip_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/chip_h1_hg38_narrowpeak.bed \
  --boundary        /home/samarth/ctcf_sequence_data/input_data/h1_human_data/4DNFIED5HLDC_sorted.bed \
  --output_dir      /home/samarth/ctcf_sequence_data/output/fimo/H1esc_CTCF_0_0 \
  --n_a 0 \
  --n_b 0 \
  --ignore_repeats

python /home/samarth/ctcf_sequence/fimo_analysis/fimo_neighbourhood_analysis.py \
  --fimo_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/fimo.tsv \
  --genome_fasta    /home/samarth/ctcf_sequence_data/input_data/genome_files/hg38.fa \
  --chip_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/chip_h1_hg38_narrowpeak.bed \
  --boundary        /home/samarth/ctcf_sequence_data/input_data/h1_human_data/4DNFIED5HLDC_sorted.bed \
  --output_dir      /home/samarth/ctcf_sequence_data/output/fimo/HFF_CTCF_0_0 \
  --n_a 0 \
  --n_b 0 \
  --ignore_repeats

python /home/samarth/ctcf_sequence/fimo_analysis/fimo_neighbourhood_analysis.py \
  --fimo_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/fimo.tsv \
  --genome_fasta    /home/samarth/ctcf_sequence_data/input_data/genome_files/hg38.fa \
  --chip_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/chip_h1_hg38_narrowpeak.bed \
  --boundary        /home/samarth/ctcf_sequence_data/input_data/h1_human_data/4DNFIED5HLDC_sorted.bed \
  --output_dir      /home/samarth/ctcf_sequence_data/output/fimo/H1esc_CTCF_20_12 \
  --n_a 20 \
  --n_b 12 \
  --ignore_repeats

python /home/samarth/ctcf_sequence/fimo_analysis/fimo_neighbourhood_analysis.py \
  --fimo_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/fimo.tsv \
  --genome_fasta    /home/samarth/ctcf_sequence_data/input_data/genome_files/hg38.fa \
  --chip_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/chip_h1_hg38_narrowpeak.bed \
  --boundary        /home/samarth/ctcf_sequence_data/input_data/h1_human_data/4DNFIED5HLDC_sorted.bed \
  --output_dir      /home/samarth/ctcf_sequence_data/output/fimo/HFF_CTCF_20_12 \
  --n_a 20 \
  --n_b 12 \
  --ignore_repeats

python /home/samarth/ctcf_sequence/fimo_analysis/fimo_neighbourhood_analysis.py \
  --fimo_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/fimo.tsv \
  --genome_fasta    /home/samarth/ctcf_sequence_data/input_data/genome_files/hg38.fa \
  --chip_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/chip_h1_hg38_narrowpeak.bed \
  --boundary        /home/samarth/ctcf_sequence_data/input_data/h1_human_data/4DNFIED5HLDC_sorted.bed \
  --output_dir      /home/samarth/ctcf_sequence_data/output/fimo/H1esc_CTCF_20_22 \
  --n_a 20 \
  --n_b 22 \
  --ignore_repeats

python /home/samarth/ctcf_sequence/fimo_analysis/fimo_neighbourhood_analysis.py \
  --fimo_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/fimo.tsv \
  --genome_fasta    /home/samarth/ctcf_sequence_data/input_data/genome_files/hg38.fa \
  --chip_filepath   /home/samarth/ctcf_sequence_data/input_data/h1_human_data/chip_h1_hg38_narrowpeak.bed \
  --boundary        /home/samarth/ctcf_sequence_data/input_data/h1_human_data/4DNFIED5HLDC_sorted.bed \
  --output_dir      /home/samarth/ctcf_sequence_data/output/fimo/HFF_CTCF_20_22 \
  --n_a 20 \
  --n_b 22 \
  --ignore_repeats






