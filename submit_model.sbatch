#!/bin/bash

# Job name:
#SBATCH --job-name=h18v01
#
# Project:
#SBATCH --account=nn9606k
#
# Wall clock limit:
#SBATCH --time=50:0:0
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=4G
#
# Number of tasks (cores):
#SBATCH  --ntasks=20

##SBATCH --qos=devel


module purge
module load intel/2019a
module load MATLAB/2019a
export MLM_LICENSE_FILE=1719@lisens-prod01.uio.no


srun --mpi=pmi2 matlab -nodisplay -nodesktop < run_CG_ESA_CCI.m
