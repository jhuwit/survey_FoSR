#!/bin/bash
#SBATCH --job-name=NI_Origial_survey
#SBATCH --array=1-200
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sgao57@jhu.edu
#SBATCH -e eofiles/NI/Origial_survey_%A_%a.err
#SBATCH -o eofiles/NI/Origial_survey_%A_%a.out

module load R
Rscript NI_Origial_survey.R $SLURM_ARRAY_TASK_ID