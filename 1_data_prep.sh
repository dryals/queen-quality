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
#last edit   09 FEB 2026

echo "date"
echo "-----------------------"

#modules
    module purge
    module load biocontainers bcftools plink
    
    #vcf location
    vcf="/depot/bharpur/data/popgenomes/gencove/NCstate/NCstate_final.bcf.gz"
    rename=/home/dryals/ryals/queen-quality/chrsrename.txt
    chrs=$( awk '{print $1}' $rename | tr '\n' ',' )
    chrsShort=$( awk '{print $2}' $rename | tr '\n' ' ' )

    
echo "-----------------------"
# #directory setup
    mkdir -p outputs
    mkdir -p $CLUSTER_SCRATCH/queen-quality
    mkdir -p blup
#     
# #filter and prepare vcf
#     
    cd $CLUSTER_SCRATCH/queen-quality
    
    #TODO: investigate duplicated samples: QC2573, QC3371
    #TODO: ensure all swaps and incorrect names are corrected!
    #TODO: try removing missing data before calling bialleleic sites, might retain more that way!
    
    #keep no contigs, only bialleleci snps, remove duplicats (norm), rename chrs
    echo "filtering sample vcf..."
    bcftools view $vcf -v snps -r $chrs -Ou | bcftools norm -m +snps -Ou | \
        bcftools view -m2 -M2 -Ou | \
        bcftools annotate --rename-chrs $rename --threads $SLURM_NTASKS -Ob -o samples.bcf.gz
    
    bcftools index -c samples.bcf.gz
    
    #remove non-QC samples
    bcftools query samples.bcf.gz -l > samples.names
    grep "QC" samples.names > keep.names
        
    #mark low propability as missing 
    echo "    mark missing..."
    bcftools view samples.bcf.gz -S keep.names -Ou | 
        bcftools filter  -S . -i 'GP[:0] > 0.99 | GP[:1] > 0.99 | GP[:2] > 0.99' \
        --threads $SLURM_NTASKS -Ob -o samples.missing.bcf.gz
    
#     #remove non-QC samples
#     bcftools query samples.missing.bcf.gz -l > samples.names
#     grep "QC" samples.names > keep.names
#     paste  keep.names  keep.names >  keep.plink
        
    #filter in plink
    echo "    mind, geno, and maf filters..."
    cd $CLUSTER_SCRATCH/queen-quality
    mkdir -p plink
    #takes LONG time
    plink --bcf samples.missing.bcf.gz --make-bed \
        --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
        --set-missing-var-ids @:# \
        #-keep keep.plink \
        --mind 0.2 --geno 0.1 --maf 0.01 \
        --threads $SLURM_NTASKS --out plink/samples-filter --silent
        
    #output sites for ref filter, samples
    cd plink
    awk '{print $2}' samples-filter.bim | tr ":" "\t" > samples-filter.sites
    awk '{print $1}' samples-filter.fam > samples-filter.names
    
# echo "-----------------------"
#     echo "LD pruning..."
#     echo "    calculating LD and af..."
#     cd ${CLUSTER_SCRATCH}/queen-quality/plink
#         plink --bfile samples-filter \
#             -r2 --ld-window 1000 --ld-window-kb 20 --ld-window-r2 0.5 \
#             --make-bed --threads $SLURM_NTASKS --out samples-preprune --silent
#             
#          plink --bfile samples-filter --freq --silent --out samples-preprune
#     
#     echo "    starting array pruning..."
#     cd $CLUSTER_SCRATCH/queen-quality/plink
#     mkdir -p ld
#     cd ~/ryals/queen-quality
#     echo -n "" > outputs/mafprune.out
#     #run R script to generate best set of sites maf
#     #start
#     sbatch --array=1-16 maf_prune_array.sh
#     #wait
#     #WARNING: collision problem, need to flock?
#     echo "    waiting for pruning (see mafprune.out)..."
#     while [ $(grep "FINISHED" outputs/mafprune.out| wc -l | awk '{print $1}') -lt 16 ] #wait for all 16 to finish
#     do
#         sleep 20 #wait between each check
#     done
#     #create full output
#     echo "    compiling results..."
#     cd ${CLUSTER_SCRATCH}/queen-quality/ld
#     #this will hold all the sites to remove
#     rm allMAFremove.txt
#     cat chr*/MAFremove.txt > allMAFremove.txt
#     count=$( wc -l allMAFremove.txt | awk '{print $1}')
#     
#         #kill script if the above fails
#         if [ $count -eq 0  ]; then
#             echo -e "\e[30;41m Filters Failed! \e[0m"
#             exit 1
#         fi
#     
#     
#     echo "    marked $count sites"
#     echo "    removing sites..."
#         cd $CLUSTER_SCRATCH/queen-quality/plink
#         plink --bfile samples-filter --make-bed --exclude ../ld/allMAFremove.txt \
#             --threads $SLURM_NTASKS --silent --out samples-pruned
#     
# echo "-----------------------"
# #PCA and GRM
#     cd $CLUSTER_SCRATCH/queen-quality/plink
#     echo "PCA..."
#     plink --bfile samples-pruned --pca 500 \
#         --threads $SLURM_NTASKS --out samples-pca --silent
#     
#     echo "GRM..."
#     #is plink the best? KING? going with basic make-rel for now
#         #TODO: try KING, compare AIC from aireml
#     module purge
#     module load biocontainers plink2
#     
#     plink2 --bfile samples-pruned -make-rel square \
#     --threads $SLURM_NTASKS --out samples-pruned --silent
#     
#     module purge
#     module load biocontainers r
#     
# echo "-----------------------"
# echo "preparing phenotypic data in R..."
#     cd ~/ryals/queen-quality
#     R --vanilla --no-save --no-echo --silent < pheno_adjust.R
# 
# 
# echo "-----------------------"
# echo "running GWAS..."
#     echo "    gcta..."
#     cd $CLUSTER_SCRATCH/queen-quality
#     mkdir -p gwas
#     cd gwas
#     gcta=/depot/bharpur/apps/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
#     
#         $gcta --bfile ../plink/samples-pruned --make-grm --thread-num $SLURM_NTASKS \
#             --autosome-num 16 --out qq
#             
#         #TODO: missing individuals here
#             
#         echo "    mlma..."
#         #adjusted phenotypes generated in R script...
#         $gcta --mlma --bfile ../plink/samples-pruned --grm qq \
#             --pheno ~/ryals/queen-quality/data/qq_vsperm.pheno \
#             --autosome-num 16 \
#             --out qq_vsperm --thread-num $SLURM_NTASKS
#             
#         $gcta --mlma --bfile ../plink/samples-pruned --grm qq \
#             --pheno ~/ryals/queen-quality/data/qq_weight.pheno \
#             --autosome-num 16 \
#             --out qq_weight --thread-num $SLURM_NTASKS
#             
#         $gcta --mlma --bfile ../plink/samples-pruned --grm qq \
#             --pheno ~/ryals/queen-quality/data/qq_lsperm.pheno \
#             --autosome-num 16 \
#             --out qq_lsperm --thread-num $SLURM_NTASKS

#TODO: do rest of traits, estimate variance explained by sig QTL
    #hit = 3 6923973
    cd $CLUSTER_SCRATCH/queen-quality/plink
    grep -n "3:6923973" samples-pruned.bim
    44939
    plink --bfile samples-pruned --recode oxford --snp '3:6923973' --out test


echo "-----------------------"  
echo "creating GRM for BLUP..."
    cd $CLUSTER_SCRATCH/queen-quality/plink
#     plink --bfile samples-filter -make-rel square \
#         --threads $SLURM_NTASKS --out samples-filter --silent
#         
    module purge
    module load biocontainers plink2
    
#     plink2 --bfile samples-filter -make-king square \
#     --threads $SLURM_NTASKS --out samples-filter --silent
    
    plink2 --bfile samples-filter -make-king square \
        --out samples-filter
    
    module purge
    module load biocontainers r




echo "-----------------------"  
echo "running BLUP..."

    #TODO: single-trait blups

    par=w.par
    cp params/${par}0 blup
    cd blup
    #aireml
    aireml=/depot/bharpur/apps/blupf90/airemlf90
    $aireml ${par}0 #&> lastrun.log
    #./airemelf90 $par
    
    #TODO
    #read G and R matricies into blup.par2
#     sed -n 16,80p file1>patch
#     sed -i 18rpatch file2
    
    cp ../params/${par}1 .
    blup=/depot/bharpur/apps/blupf90/blupf90+
    $blup ${par}1
    cp solutions ../data/sol-12feb26.txt
#     
#     
#     
#TODO: estimate CV error: scripts/cv.R

    

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
