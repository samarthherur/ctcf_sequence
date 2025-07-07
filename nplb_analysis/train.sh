#!/bin/bash
#SBATCH --job-name=nplb_train_all
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=7-23:00:00          # 7 days, 23 hours
#SBATCH --partition=gpu
#SBATCH --nodes=1                  # <-- one node only
#SBATCH --ntasks=80                 # <-- eighty “slots” for our 80 jobs
#SBATCH --nodelist=cn1             # <-- pin to the node that works

cd $SLURM_SUBMIT_DIR
module load NPLB-1.0.0 python-2.7 ghostscript-10.03.1

############################################
# Define your eight FASTA→output pairs in order:
############################################
declare -a FASTA_FILES=(
  "/home/samarth/ctcf_sequence_data/output/fimo/H1esc_CTCF_0_0/fimo_neighbourhood_0_0_filtered_wr.fasta"
  "/home/samarth/ctcf_sequence_data/output/fimo/H1esc_CTCF_20_12/fimo_neighbourhood_20_12_filtered_wr.fasta"
  "/home/samarth/ctcf_sequence_data/output/fimo/H1esc_CTCF_20_22/fimo_neighbourhood_20_22_filtered_wr.fasta"
  "/home/samarth/ctcf_sequence_data/output/fimo/H1esc_CTCF_50_50/fimo_neighbourhood_50_50_filtered_wr.fasta"
  "/home/samarth/ctcf_sequence_data/output/fimo/HFF_CTCF_0_0/fimo_neighbourhood_0_0_filtered_wr.fasta"
  "/home/samarth/ctcf_sequence_data/output/fimo/HFF_CTCF_20_12/fimo_neighbourhood_20_12_filtered_wr.fasta"
  "/home/samarth/ctcf_sequence_data/output/fimo/HFF_CTCF_20_22/fimo_neighbourhood_20_22_filtered_wr.fasta"
  "/home/samarth/ctcf_sequence_data/output/fimo/HFF_CTCF_50_50/fimo_neighbourhood_50_50_filtered_wr.fasta"
)

declare -a OUT_DIRS=(
  "/home/samarth/ctcf_sequence_data/output/NPLB/H1esc_0_0_train/"
  "/home/samarth/ctcf_sequence_data/output/NPLB/H1esc_20_12_train/"
  "/home/samarth/ctcf_sequence_data/output/NPLB/H1esc_20_22_train/"
  "/home/samarth/ctcf_sequence_data/output/NPLB/H1esc_50_50_train/"
  "/home/samarth/ctcf_sequence_data/output/NPLB/HFF_0_0_train/"
  "/home/samarth/ctcf_sequence_data/output/NPLB/HFF_20_12_train/"
  "/home/samarth/ctcf_sequence_data/output/NPLB/HFF_20_22_train/"
  "/home/samarth/ctcf_sequence_data/output/NPLB/HFF_50_50_train/"
)

# Launch them all in parallel
for i in "${!FASTA_FILES[@]}"; do
  promoterLearn \
    -f "${FASTA_FILES[$i]}" \
    -o "${OUT_DIRS[$i]}" & 
done

# Wait for all 8 to finish
wait