#!/bin/bash
#SBATCH --job-name=NI_for_visualization_b300
#SBATCH --array=1-2
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sgao57@jhu.edu
#SBATCH -e eofiles/NI/NI_for_visualization_b300_%A_%a.err
#SBATCH -o eofiles/NI/NI_for_visualization_b300_%A_%a.out

module load R
Rscript NI_for_visualization.R $SLURM_ARRAY_TASK_ID