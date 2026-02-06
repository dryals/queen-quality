#!/bin/bash

# FILENAME: 1_data_prep.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --partition cpu
#SBATCH --time=1-00:00:00
#SBATCH --job-name data_prep

#SBATCH --output=/home/dryals/ryals/queen-quality/outputs/prep.out
#SBATCH --error=/home/dryals/ryals/queen-quality/outputs/prep.out

#Dylan Ryals 06 FEB 2026
#last edit

echo "date"
echo "-----------------------"

#modules
    module purge
    module load biocontainers bcftools plink

#directory setup
    mkdir -p outputs
    mkdir -p $CLUSTER_SCRATCH/queen-quality
    
#filter and prepare vcf
    #vcf location
    vcf="/depot/bharpur/data/popgenomes/gencove/NCstate/NCstate_final.bcf.gz"
    rename=/home/dryals/ryals/queen-quality/chrsrename.txt
    chrs=$( awk '{print $1}' $rename | tr '\n' ',' )
    chrsShort=$( awk '{print $2}' $rename | tr '\n' ' ' )
    
    cd $CLUSTER_SCRATCH/queen-quality
    
    #TODO: investigate duplicated samples: QC2573, QC3371
    #TODO: ensure all swamps and incorrect names are corrected!
    
    #keep no contigs, only bialleleci snps, remove duplicats (norm), rename chrs
    echo "filtering sample vcf..."
    bcftools view $vcf -v snps -r $chrs -Ou | bcftools norm -m +snps -Ou | \
        bcftools view -m2 -M2 -Ou | \
        bcftools annotate --rename-chrs $rename --threads $SLURM_NTASKS -Ob -o samples.bcf.gz
    
    bcftools index -c samples.bcf.gz
    
    #mark low propability as missing 
    echo "    mark missing..."
    bcftools filter samples.bcf.gz -S . -i 'GP[:0] > 0.99 | GP[:1] > 0.99 | GP[:2] > 0.99' \
        --threads $SLURM_NTASKS -Ob -o samples.missing.bcf.gz
        
    #filter in plink
    echo "    mind, geno, and maf filters..."
    mkdir -p plink
    #takes LONG time
    plink --bcf samples.missing.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
        --set-missing-var-ids @:# \
        --mind 0.2 --geno 0.1 --maf 0.01 \
        --threads $SLURM_NTASKS --out plink/samples-filter --silent
        
    #output sites for ref filter, samples
    cd plink
    awk '{print $2}' samples-filter.bim | tr ":" "\t" > samples-filter.sites
    awk '{print $1}' samples-filter.fam > samples-filter.names
    
#PCA and GRM
    cd $CLUSTER_SCRATCH/queen-quality/plink
    echo "PCA..."
    plink --bfile samples-filter --maf 0.05 --pca 500 \
        --threads $SLURM_NTASKS --out samples-maf --silent
    
    echo "GRM..."
    #is plink the best? KING? going with basic make-rel for now
    module purge
    module load biocontainers plink2
    
    plink2 --bfile samples-filter -make-rel square \
    --threads $SLURM_NTASKS --out samples-filter --silent
    
    
#TODO: admixture components


        

# #pull most recent plink file...
#     cd $CLUSTER_SCRATCH/pipeline/plink
#     cp US1pc.* $CLUSTER_SCRATCH/blup/plink
#     
#     #pull samples
#     cd $CLUSTER_SCRATCH/blup/plink
#     cat US1pc.fam | awk '{print $1, $2}' | grep 'QC' > tarpy.samps
#     
#     #TODO: consider filtering options, try to get more sites to start, try real LD pruning 
#     
#     plink2 --bfile US1pc --make-bed --keep tarpy.samps \
#     --maf 0.01 --make-king square --out tarpy
# 
#     plink2 --bfile tarpy --make-bed \
#     --bp-space 5000 --make-king square --pca 60 --out tarpy.small
#     
#     
# #GWAS time
#     #TODO: incorporate ancestry, ensure that order is the same in all files
#     echo "gcta..."
#     cd $CLUSTER_SCRATCH/blup
#     mkdir -p gwas
#     cd gwas
#     gcta=/depot/bharpur/apps/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
#     
#         $gcta --bfile ../plink/tarpy.small --make-grm --thread-num $SLURM_NTASKS --autosome-num 16 --out tarpy
#         echo "mlma..."
#         #adjusted phenotypes generated in R script...
#         $gcta --mlma --bfile ../plink/tarpy.small --grm tarpy --pheno tarpy_viability.pheno --autosome-num 16 \
#             --out tarpy_viability --thread-num $SLURM_NTASKS
#             
#         $gcta --mlma --bfile ../plink/tarpy.small --grm tarpy --pheno tarpy_weight.pheno --autosome-num 16 \
#             --out tarpy_weight --thread-num $SLURM_NTASKS
# 
# 
#         
echo "-----------------------"
echo "DONE"
echo "-----------------------"
