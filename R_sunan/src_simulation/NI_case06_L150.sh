#!/bin/bash
#SBATCH --job-name=NI_case06_L150
#SBATCH --array=1-200
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sgao57@jhu.edu
#SBATCH -e eofiles/NI/NI_case06_L150_%A_%a.err
#SBATCH -o eofiles/NI/NI_case06_L150_%A_%a.out

module load R
Rscript NI_case06_L150.R $SLURM_ARRAY_TASK_ID