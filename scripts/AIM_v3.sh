#!/bin/bash

# FILENAME: AIM_v3.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --job-name AIM_v3_array
#SBATCH --partition cpu
#SBATCH --output=/home/dryals/ryals/queen-quality/outputs/dump.out
#SBATCH --error=/home/dryals/ryals/queen-quality/outputs/dump.out

#Dylan Ryals 26 OCT 2023
#last edit   13 FEB 2026

#this version formatted for quaan-quality project

#this version reads a filename from a file to be more modular
    #this could be some sort of command line argument or named pipe, but eh
    
#updating for response, new slurm

#launch with sbatch --array=1-16 AIM_v3.sh

date

module load biocontainers bcftools r

NCHR=$( echo $SLURM_ARRAY_TASK_ID )
log=/depot/bharpur/data/projects/ryals/queen-quality/outputs/aim.out

cd ${CLUSTER_SCRATCH}/queen-quality/aim
filename=$( cat ref_filename.txt )
echo "reading $filename ..."

cd ${CLUSTER_SCRATCH}/queen-quality/aim
mkdir -p chr${NCHR}
cd chr${NCHR}

echo "starting chr $NCHR" >> $log

bcftools view ../../$filename -r $NCHR -Ob -o chr${NCHR}refs.bcf.gz
bcftools index -c chr${NCHR}refs.bcf.gz

for pop in A C M O
#for pop in A C L M O U Y
do
    echo "starting $NCHR $pop ..."
    bcftools view chr${NCHR}refs.bcf.gz -S ~/ryals/queen-quality/references/${pop}.txt -Ou | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%AF\n' -o ${pop}.frq
    
    awk '{print $3}' ${pop}.frq > ${pop}.tmp
done

paste A.frq C.tmp M.tmp O.tmp > chr${NCHR}.popfrq
#paste A.frq C.tmp L.tmp M.tmp O.tmp U.tmp Y.tmp > chr${NCHR}.popfrq
rm *.tmp *.frq

echo "calculating Ia for chr $NCHR" >> $log

Rscript --vanilla --silent ~/ryals/queen-quality/scripts/aimIa_v2.R $NCHR
    #v2 runtime: ~3:40
    #v1 runtime: ~9:30

#report to logs, using flock to avoid conflicts
( flock -x 9 
    echo "FINISHED CHR $NCHR" >> $log
    echo "    FINISHED CHR $NCHR" >> ~/ryals/queen-quality/outputs/prep.out
) 9> ~/ryals/queen-quality/locks/.AIMwritelock
echo "done"
date







