#!/bin/bash
#SBATCH --job-name=NI_noOrigial_survey_b500
#SBATCH --array=1-200
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sgao57@jhu.edu
#SBATCH -e eofiles/NI/noOrigial_survey_b500_%A_%a.err
#SBATCH -o eofiles/NI/noOrigial_survey_b500_%A_%a.out

module load R
Rscript NI_noOrigial_survey_b500.R $SLURM_ARRAY_TASK_ID