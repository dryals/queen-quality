#!/bin/bash

# FILENAME: prune_array.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --partition=cpu
#SBATCH --time=2:00:00
#SBATCH --job-name prune_array.sh

#SBATCH --output=/home/dryals/ryals/queen-quality/outputs/dump_prune.out
#SBATCH --error=/home/dryals/ryals/queen-quality/outputs/dump_prune.out

#Dylan Ryals 21 JAN 2025
#last edited 16 FEB 2026

date
echo "---------------"
module load r

cd $CLUSTER_SCRATCH/queen-quality/aim
n=$( echo $SLURM_ARRAY_TASK_ID )
log=/depot/bharpur/data/projects/ryals/queen-quality/outputs/prune.out

cd ~/ryals/queen-quality
echo "starting chr $n..." >> $log

Rscript --vanilla --silent scripts/LDprune_p_v2.R $n

( flock -x 9 
    echo "FINISHED CHR $n" >> $log
    echo "    FINISHED CHR $n" >> ~/ryals/queen-quality/outputs/prep.out
) 9> ~/ryals/queen-quality/locks/.PRUNEwritelock

echo "---------------"
date
