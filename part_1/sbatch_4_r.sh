#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50g
#SBATCH --time=24:00:00
#SBATCH --job-name=list2dist
#SBATCH --output=/gpfs01/home/mbxjk6/slurm-%x-%j.out

module load rstudio-uon/gcc11.3.0/2022.07.1-554
Rscript Path/to/list2dist.R 
