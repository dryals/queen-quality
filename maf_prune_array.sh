#!/bin/bash

# FILENAME: maf_prune_array.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=4
#SBATCH --partition cpu
#SBATCH --time=2:00:00
#SBATCH --job-name maf_prune_array.sh

#SBATCH --output=/home/dryals/ryals/queen-quality/outputs/dump_mafprune.out
#SBATCH --error=/home/dryals/ryals/queen-quality/outputs/dump_mafprune.out

#Dylan Ryals 21 JAN 2025
#last edit   09 FEB 2026

#fixed partition

date
echo "---------------"
module load r

n=$( echo $SLURM_ARRAY_TASK_ID )
log=/depot/bharpur/data/projects/ryals/queen-quality/outputs/mafprune.out

cd $CLUSTER_SCRATCH/queen-quality/ld

mkdir -p chr${n}

cd ~/ryals/queen-quality
echo "starting chr $n..." >> $log

Rscript --vanilla --silent maf_prune.R $n

echo "FINISHED CHR $n" >> $log

echo "---------------"
date
