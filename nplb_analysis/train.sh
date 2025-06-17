#!/bin/bash
#SBATCH --job-name=nplb_train
#SBATCH --output=nplb_train.out
#SBATCH --error=nplb_train.err
#SBATCH --time=5-10:00:00
#SBATCH --partition=gpu
#SBATCH --ntasks=40
#SBATCH --nodelist=cn1

cd $SLURM_SUBMIT_DIR
module load NPLB-1.0.0 python-2.7 ghostscript-10.03.1
promoterLearn -f /home/samarth/ctcf_sequence_data/output/fimo/HFF/fimo_neighbourhood_50_filtered_wr.fasta -o /home/samarth/ctcf_sequence_data/NPLB/HFF/ 

