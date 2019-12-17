#!/bin/bash

#SBATCH -D /home/sstuenzi/CryoVeg/CryoVegModel
#SBATCH -o slurm_output_%A_%a.out
#SBATCH -e slurm_error_%A_%a.err
#SBATCH -n 5 
#SBATCH -N 1 
#SBATCH -p PermaRisk
#SBATCH -J Simone
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=simone.stuenzi@awi.de

#srun matlab -nodesktop -r main
#matlab -nodisplay -nosplash -nojvm -r main.m
matlab -nodisplay -nosplash -nodesktop -r "run('main.m');exit;"

matlab -nosplash -nodisplay -nodesktop -r "run('main.m');exit;"
