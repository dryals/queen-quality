#!/bin/bash

# FILENAME: supervised_admix_v3.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --partition=cpu
#SBATCH --time=3:00:00
#SBATCH --job-name supadmix_v3

#SBATCH --output=/home/dryals/ryals/queen-quality/outputs/supadmix.out
#SBATCH --error=/home/dryals/ryals/queen-quality/outputs/supadmix.out

#Dylan Ryals 27 DEC 2022
#last edited 16 FEB 2026

#this version for qq data


#admixture

date

cd ${CLUSTER_SCRATCH}/queen-quality/admix/supervised
#.pop file is needed to determine populations, created with makeAdmixPop.R script
filename=$( cat ${CLUSTER_SCRATCH}/queen-quality/plink/plink_admix_filename.txt )

    echo "running admixture..."
    
        ADMIX=/depot/bharpur/apps/admixture/admixture
        #remember to change this for the correct k value!!
        $ADMIX ../../plink/${filename}.bed 4 -j8 --cv=20 --supervised > ${filename}.out
        
date

echo "DONE"
