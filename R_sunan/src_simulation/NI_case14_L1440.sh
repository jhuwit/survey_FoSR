#!/bin/bash
#SBATCH --job-name=NI_case14_L1440
#SBATCH --array=166-200
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sgao57@jhu.edu
#SBATCH -e eofiles/NI/NI_case14_L1440_%A_%a.err
#SBATCH -o eofiles/NI/NI_case14_L1440_%A_%a.out

module load R
Rscript NI_case14_L1440.R $SLURM_ARRAY_TASK_ID