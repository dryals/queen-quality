#!/bin/bash

# FILENAME: GSprep.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --time=1-00:00:00
#SBATCH --job-name GSprep

#SBATCH --output=/home/dryals/ryals/blup/prep.out
#SBATCH --error=/home/dryals/ryals/blup/prep.out

#Dylan Ryals 27 MAY 2025
#last edit   
    
    #creating a GRM from vcf data...
    
#modules
module load biocontainers bcftools plink
    
#pull most recent plink file...
    cd $CLUSTER_SCRATCH/pipeline/plink
    cp US1pc.* $CLUSTER_SCRATCH/blup/plink
    
    #pull samples
    cd $CLUSTER_SCRATCH/blup/plink
    cat US1pc.fam | awk '{print $1, $2}' | grep 'QC' > tarpy.samps
    
    #TODO: consider filtering options, try to get more sites to start, try real LD pruning 
    
    plink2 --bfile US1pc --make-bed --keep tarpy.samps \
    --maf 0.01 --make-king square --out tarpy

    plink2 --bfile tarpy --make-bed \
    --bp-space 5000 --make-king square --pca 60 --out tarpy.small
    
    
#GWAS time
    #TODO: incorporate ancestry, ensure that order is the same in all files
    echo "gcta..."
    cd $CLUSTER_SCRATCH/blup
    mkdir -p gwas
    cd gwas
    gcta=/depot/bharpur/apps/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
    
        $gcta --bfile ../plink/tarpy.small --make-grm --thread-num $SLURM_NTASKS --autosome-num 16 --out tarpy
        echo "mlma..."
        #adjusted phenotypes generated in R script...
        $gcta --mlma --bfile ../plink/tarpy.small --grm tarpy --pheno tarpy_viability.pheno --autosome-num 16 \
            --out tarpy_viability --thread-num $SLURM_NTASKS
            
        $gcta --mlma --bfile ../plink/tarpy.small --grm tarpy --pheno tarpy_weight.pheno --autosome-num 16 \
            --out tarpy_weight --thread-num $SLURM_NTASKS


        
