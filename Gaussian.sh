#!/bin/bash

#SBATCH -c 1
#SBATCH -t 10-00:00:00
#SBATCH -p general
#SBATCH -q grp_mlynch11
#SBATCH -J GaussianBatch
#SBATCH --mem-per-cpu=5G
#SBATCH --exclude=pcc003

echo "Running the script with $SLURM_ARRAY_TASK_ID"

module load gsl-2.7.1-dm
module load intel-oneapi-compilers-2023.2.1-5l

./Gaussian $SLURM_ARRAY_TASK_ID
