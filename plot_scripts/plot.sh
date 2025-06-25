#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --output=plot.out
#SBATCH --error=plot.err
#SBATCH --time=5-10:00:00
#SBATCH --partition=gpu
#SBATCH --ntasks=40
#SBATCH --nodelist=cn1

cd $SLURM_SUBMIT_DIR
python plot_heatmaps.py \
  -n /home/samarth/final_plots_thesis/preliminary_processing/nplb_clustered.bed \
  -r /home/samarth/hg38/hg38.fa \
  -i /home/samarth/dekker_microc/H1ESC/4DNFI6YRK53R_insulation.bw \
  -c /home/samarth/ChIP_Seq/ENCFF041TAK_chip_h1_hg38_signal.bw \
  -p /home/samarth/phastcons/hg38.phastCons100way.bw \
  -e /home/samarth/ChIP_Seq/epi_h1esc_list