#!/bin/bash
#SBATCH --job-name=counts_xplr3
#SBATCH --time=06:30:00
#SBATCH --mem=70G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl
#SBATCH --error=slurm-%j.err
#SBATCH --output=slurm-%j.out

# Correct path to conda.sh
source /hpc/shared/prekovic/dhaynessimmons/miniconda3/etc/profile.d/conda.sh

# # Activate the Conda environment containing R
# conda init
conda activate r_env

# Now, run the R script. R is "activated" in the sense that
# the version of R and the packages available are those installed in `r_env`.
Rscript tcga_cancer_counts_xplr_v4.r