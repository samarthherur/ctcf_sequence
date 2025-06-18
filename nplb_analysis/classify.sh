
#!/bin/bash
#SBATCH --job-name=nplb_classify
#SBATCH --output=nplb_classify.out
#SBATCH --error=nplb_classify.err
#SBATCH --time=10:00:00
#SBATCH --partition=gpu
#SBATCH --ntasks=40
#SBATCH --nodelist=cn12

cd $SLURM_SUBMIT_DIR
module load NPLB-1.0.0 python-2.7 ghostscript-10.03.1
promoterClassify -f /home/samarth/hg38/fimo_ctcf_hg38_hFFc6/fimo_neighbourhood_50_filtered_wr.fasta -m /home/samarth/NPLB/output_boundary_ctcf_concat_corrected_wr/bestModel.p -o /home/samarth/NPLB/output_concat_corrected_hff_wr/ 
