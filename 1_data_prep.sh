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
#last edit   04 MAR 2026

echo "date"
echo "-----------------------"

#modules
    module purge
    module load biocontainers bcftools vcftools plink r
    
    #vcf location
    vcf="/depot/bharpur/data/popgenomes/gencove/NCstate/NCstate_final2.bcf.gz"
    rename=/home/dryals/ryals/queen-quality/chrsrename.txt
    chrs=$( awk '{print $1}' $rename | tr '\n' ',' )
    chrsShort=$( awk '{print $2}' $rename | tr '\n' ' ' )
    refs=/depot/bharpur/data/popgenomes/HarpurPNAS/output_snp.vcf.gz

echo "-----------------------"
#     #directory setup
#     mkdir -p outputs
#     mkdir -p locks
#     mkdir -p $CLUSTER_SCRATCH/queen-quality
#     mkdir -p blup 
#     #filter and prepare vcf  
#     cd $CLUSTER_SCRATCH/queen-quality


#TODO: Luiz Suggestions
    #filtering: 
        #fixed effect of season?
    #GWAS: 
        #fit ancestry (test against PCA)
        #clumping
        #var exp from snp eff
        #independent chromosomal segments for MHT
    #GS: 
        #dont standardize phenos
        #better estimates for reml
        #add small value to diagonal, check PD
        #check mean diag and off-diag values
        #EM-reml x 20
    #analysis: 
        #BV trend over years...
        

     
    #TODO: try removing missing data before calling bialleleic sites, might retain more that way!
    #TODO: why are there two LDprune scripts?
    
    #TODO: make sure GRM filter hapens once only, same samples used for all analyses
    
    #TODO: compare various methods of GRM construction, verify everything loaded into BLUP correctly...
    #TODO: fit PCs to BLUP as well as GWAS ... ???
    
#     
#     #keep no contigs, only bialleleci snps, remove duplicats (norm), rename chrs
#     echo "filtering sample vcf..."
#     bcftools view $vcf -v snps -r $chrs -Ou | bcftools norm -m +snps -Ou | \
#         bcftools view -m2 -M2 -Ou | \
#         bcftools annotate --rename-chrs $rename --threads $SLURM_NTASKS -Ob -o samples.bcf.gz
#     
#     bcftools index -c samples.bcf.gz
#     
#     #remove non-QC samples
#     bcftools query samples.bcf.gz -l > samples.names
#     grep "QC" samples.names > keep.names
#         
# 
# #TODO: increase post prob threshold?
#     #mark low propability as missing 
#     echo "    mark missing..."
#     bcftools view samples.bcf.gz -S keep.names -Ou | 
#         bcftools filter  -S . -i 'GP[:0] > 0.995 | GP[:1] > 0.995 | GP[:2] > 0.995' \
#         --threads $SLURM_NTASKS -Ob -o samples.missing.bcf.gz 
#      
#     #filter in plink
#     echo "    mind, geno, and maf filters..."
#     cd $CLUSTER_SCRATCH/queen-quality
#     mkdir -p plink
#     #takes LONG time
#     plink --bcf samples.missing.bcf.gz --make-bed \
#         --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
#         --set-missing-var-ids @:# \
#         --mind 0.2 --geno 0.1 --maf 0.01 \
#         --hwe 1e-8 \
#         --threads $SLURM_NTASKS --out plink/samples-filter --silent
#         
#     #output sites for ref filter, samples
#     cd plink
#     awk '{print $2}' samples-filter.bim | tr ":" "\t" > samples-filter.sites
#     awk '{print $1}' samples-filter.fam > samples-filter.names
# 
# echo "-----------------------"
# echo "preparing phenotypic data in R..."
#     cd ~/ryals/queen-quality
#     R --vanilla --no-save --no-echo --silent < scripts/clean_pheno.R     

#further sample QC
# echo "-----------------------"
#     echo "removing GRM outliers"
#     cd $CLUSTER_SCRATCH/queen-quality/plink
#     module purge
#     module load biocontainers plink2
#     
#     plink2 --bfile samples-filter -make-rel square --out samples-filter
#     #--remove ../het.remove.plink \
#         
#     
#     module purge
#     module load biocontainers bcftools vcftools plink r
            
#    
#    echo "reading GRM..."
#             
#     R
#         ReadGRMBin=function(prefix, AllN=F, size=4){
#             sum_i=function(i){
#                 return(sum(1:i))
#             }
#                 BinFileName=paste(prefix,".grm.bin",sep="")
#                 NFileName=paste(prefix,".grm.N.bin",sep="")
#                 IDFileName=paste(prefix,".grm.id",sep="")
#                 id = read.table(IDFileName)
#                 n=dim(id)[1]
#                 BinFile=file(BinFileName, "rb");
#                 grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
#                 NFile=file(NFileName, "rb");
#                 if(AllN==T){
#                     N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
#                 }
#                 else N=readBin(NFile, n=1, what=numeric(0), size=size)
#                 i=sapply(1:n, sum_i)
#                 return(list(diag=grm[i], off=grm[-i], id=id, N=N))
#             }
#         
#         G = ReadGRMBin("samples-gs")
#     
#         remove = G$id[G$diag > 1.7,1]
#         
#         remove
#         
#         write.table(remove, file = "GRM.remove",
#                     row.names = F, col.names = F, quote = F)
#                     
#     quit(save = "no")
#     
#     paste GRM.remove GRM.remove > GRM.remove.plink
#     
# 
#         
# echo "-----------------------"
#     echo "ADMIXTURE analysis"
#     cd ${CLUSTER_SCRATCH}/queen-quality
# #          echo "pulling references..."
# #         #no multiallelic sites, only snps, keep subset of references, no contigs, rename chromosomes to "1,2,3...16"
# #         bcftools view $refs -S /home/dryals/ryals/ahb/references/pureRefs.txt \
# #             -r $chrsLong -M2 -v snps -Ou | \
# #             bcftools annotate --rename-chrs $rename --threads $SLURM_NTASKS \
# #             -Ob -o reference.bcf.gz
# #             
# #         bcftools index -c reference.bcf.gz
# 
#     #copy from admix dir
#     cp ../pipeline/reference.bcf.* .
#     
#     echo "    filtering references ..."
#         bcftools view reference.bcf.gz -T plink/samples-filter.sites --threads $SLURM_NTASKS \
#         -Ob -o reference-filter.bcf.gz
# 
#      bcftools index -c reference-filter.bcf.gz
#      
#     
#     echo "launching Ia script...."
#         #count number of samples in each population
#         cd ${CLUSTER_SCRATCH}/queen-quality
#         mkdir -p aim
#         cd aim
#         #specify reference file
#         echo "reference-filter.bcf.gz" > ref_filename.txt
#         #reset logifle
#         cd /home/dryals/ryals/queen-quality
#         echo -n "" > outputs/aim.out
# 
#         #launch the Ia script array
#         sbatch --array=1-16 scripts/AIM_v3.sh
#         
#     echo "waiting for Ia results (see aim.out)..."
#     cd ~/ryals/queen-quality
#     while [ $(grep "FINISHED" outputs/aim.out | wc -l | awk '{print $1}') -lt 16 ] #wait for all 16 to finish
#     do
#         sleep 20 #wait between each check
#     done
#     
#     echo "compiling Ia results..."    
#     cd ${CLUSTER_SCRATCH}/queen-quality/aim
#     #this will hold all the aims
#     cat chr*/chr*.ia | grep -v "chr" | sort -k3 -gr > aim.ia.txt
#    
#         grep -v "NA" aim.ia.txt | awk '$3>0' | awk 'OFS=":" {print$1, $2}' > plink_aim.ia.txt
#           cat plink_aim.ia.txt | tr ":" "\t" > bct_aim.ia.txt
#           
#         
#     count=$( wc -l plink_aim.ia.txt | awk '{print $1}')
#     echo "    Calculated Ia for $count sites"
#    
#     echo "merging samples and references..."
#     cd $CLUSTER_SCRATCH/queen-quality
#     bcftools view samples.missing.bcf.gz -T plink/samples-filter.sites -S plink/samples-filter.names \
#         --threads $SLURM_NTASKS -Ob -o samples-aim.bcf.gz
#         
#         bcftools index -c samples-aim.bcf.gz
#         
#     bcftools merge reference-filter.bcf.gz samples-aim.bcf.gz -m snps -Ou | \
#         bcftools norm -m +snps -Ou | bcftools view -M2 -m2 --threads $SLURM_NTASKS -Ob -o qq-admix.bcf.gz
#     echo "  indexing..."
#         bcftools index -c qq-admix.bcf.gz
#     
#     echo "starting PLINK for supervised data" 
#         #aims, maf, and geno
#         cd ${CLUSTER_SCRATCH}/queen-quality
#             echo "    pulling informative sites..."
#             plink --bcf qq-admix.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
#                 --set-missing-var-ids @:# --silent \
#                 --extract aim/plink_aim.ia.txt --threads $SLURM_NTASKS --out plink/topaim
# 
#             #extract ld data, removing references
#             echo "    calculating ld..."
#             cd plink
#             #TODO: consider maf threshold
#             plink --bfile topaim -r2 --ld-window 1000 --ld-window-kb 50 --ld-window-r2 0.2 \
#                 --remove /home/dryals/ryals/admixPipeline/references/plink_pureRefs.txt \
#                 --silent --threads $SLURM_NTASKS --out sampleAncPrune
#                 
#             #run R script to generate best set of sites by Ia
#                 #reset logfile
#                 cd /home/dryals/ryals/queen-quality
#                 echo -n "" > outputs/prune.out
#                 #start
#                 sbatch --array=1-16 scripts/prune_array.sh
#                 #wait
#                 echo "    waiting for pruning (see prune.out)..."
#                 while [ $(grep "FINISHED" outputs/prune.out | wc -l | awk '{print $1}') -lt 16 ] #wait for all 16 to finish
#                 do
#                     sleep 20 #wait between each check
#                 done
#                 #create full output
#                 echo "    compiling results..."
#                 cd ${CLUSTER_SCRATCH}/queen-quality/aim
#                 #this will hold all the aims
#                 cat chr*/LDremove.txt > allLDremove.txt
#                 count=$( wc -l allLDremove.txt | awk '{print $1}')
#                 echo "    marked $count sites"
#                 
#             echo "    removing pruned sites..."
#             cd $CLUSTER_SCRATCH/queen-quality/plink
#             #create new admix file
#             plink --bfile topaim --make-bed --exclude ../aim/allLDremove.txt --silent \
#                 --threads $SLURM_NTASKS --out finaladmix
#         
#         #use this plink file basename for admix scripts
#         echo "finaladmix" > plink_admix_filename.txt
#         
#         #kill script if the above fails
#         if [ ! -f "finaladmix.bed" ]; then
#             echo -e "\e[30;41m Plink Failed! \e[0m"
#             exit 1
#         fi
#         
# echo "starting supervised admix..."  
#     cd $CLUSTER_SCRATCH/queen-quality
#     mkdir -p admix
#     cd admix 
#     mkdir -p supervised
#     #supervised
#     cd /home/dryals/ryals/queen-quality
#     #create pop file
#     R --vanilla --no-save --no-echo --silent < scripts/makeAdmixPop.R
#     sbatch scripts/supervised_admix_v3.sh
# 
# echo "-----------------------"
#     echo "preparing GWAS and GS"
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

#     echo "    removing GRM outliers..."
#     plink --bfile samples-pruned --make-bed \
#         --remove GRM.remove.plink \
#         --pca 500 \
#         --threads $SLURM_NTASKS --out samples-gwas
#         
#         
# 
# echo "-----------------------"
# #PCA and GRM
#     cd $CLUSTER_SCRATCH/queen-quality/plink
#     echo "PCA..."
#     plink --bfile samples-filter --pca 500 \
#         --threads $SLURM_NTASKS --out samples-gwas --silent
#      
#     echo "GRM..."
#     #is plink the best? KING? going with basic make-rel for now
#         #TODO: consider gcta GRM after all
#     module purge
#     module load biocontainers plink2
#     
#     plink2 --bfile samples-filter --keep ~/ryals/queen-quality/data/phenotyped.plink \
#         -maf 0.05 \
#         -make-rel square --out samples-gs2
# 
#     module purge
#     module load biocontainers bcftools vcftools plink r
#     
#     
#     echo "GRM in GCTA..."
#     cd ~/ryals/queen-quality/data
#     paste phenotyped.gcnames phenotyped.gcnames > phenotyped.plink
#     
#     cd $CLUSTER_SCRATCH/queen-quality/plink
#     
#     plink --bfile samples-filter --maf 0.01 \
#         --keep ~/ryals/queen-quality/data/phenotyped.plink \
#         --make-bed --threads $SLURM_NTASKS --out samples-gs
#         
#     gcta=/depot/bharpur/apps/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
#     
#         $gcta --bfile samples-gs --make-grm --thread-num $SLURM_NTASKS \
#             --autosome-num 16 --out samples-gs
# 
#                          
# echo "-----------------------"
# echo "preparing data for GWAS and GS"
#     cd ~/ryals/queen-quality
#     R --vanilla --no-save --no-echo --silent < scripts/adjust_pheno.R

echo "-----------------------"
echo "running GWAS..."
    echo "    gcta..."
    cd $CLUSTER_SCRATCH/queen-quality
    mkdir -p gwas
    cd gwas
    gcta=/depot/bharpur/apps/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
    
        $gcta --bfile ../plink/samples-filter --make-grm --thread-num $SLURM_NTASKS \
            --autosome-num 16 --out qq
            
        echo "    mlma..."
        #adjusted phenotypes generated in R script...
        $gcta --mlma --bfile ../plink/samples-filter --grm qq \
            --pheno ~/ryals/queen-quality/data/qq_vsperm.pheno \
            --autosome-num 16 \
            --out qq_vsperm --thread-num $SLURM_NTASKS
            
        $gcta --mlma --bfile ../plink/samples-filter --grm qq \
            --pheno ~/ryals/queen-quality/data/qq_weight.pheno \
            --autosome-num 16 \
            --out qq_weight --thread-num $SLURM_NTASKS
            
        $gcta --mlma --bfile ../plink/samples-filter --grm qq \
            --pheno ~/ryals/queen-quality/data/qq_lsperm.pheno \
            --autosome-num 16 \
            --out qq_lsperm --thread-num $SLURM_NTASKS
            

# #TODO: estimate variance explained by sig QTL
#     #hit = 3 6923973
#     cd $CLUSTER_SCRATCH/queen-quality/plink
#     grep -n "3:6923973" samples-pruned.bim
#     44939
#     plink --bfile samples-pruned --recode oxford --snp '3:6923973' --out test
# 
# 
# 
# 
# echo "-----------------------"  
# echo "running BLUP..."
# 
#     par=wv
# 
#     #TODO: single-trait blups
#     cd ~/ryals/queen-quality/blup
#         #create links
#         if [ ! -f  blupf90+ ]; then
#             echo "    creating links..."
#             ln -S blupf90+ /depot/bharpur/apps/blupf90/blupf90+
#             ln -S airemlf90 /depot/bharpur/apps/blupf90/airemlf90 
#         fi
# 
#     cd ~/ryals/queen-quality
#     cp params/${par}.par0 blup
#     cd blup
#     ./airemlf90 ${par}.par0
#     
#     
#     #TODO
#     #read G and R matricies into blup.par2
# #     sed -n 16,80p file1>patch
# #     sed -i 18rpatch file2
#     
#     cp ../params/${par}.par1 .
#     ./blupf90+ ${par}.par1
#     cp solutions ../data/sol-${par}.txt

# echo "-----------------------"
#     echo "  CV error: single-trait"
# #TODO: estimate CV error: scripts/cv.R
# 
#     #create -cv version which uses pheno-cv.txt
#     cd ~/ryals/queen-quality
#     cp params/${par}.par1 blup/${par}-cv.par1
#     sed -i 's/pheno.txt/pheno-cv.txt/g' blup/${par}-cv.par1
#     
#     #run cv script
#     Rscript --vanilla scripts/cv.R $par
# 
# echo "-----------------------"
#     echo "  CV error: multi-trait"
#     
#     par=wv
#  
#      #create -cv version which uses pheno-cv.txt
#     cd ~/ryals/queen-quality
#     cp params/${par}.par1 blup/${par}-cv.par1
#     sed -i 's/pheno.txt/pheno-cv.txt/g' blup/${par}-cv.par1
#     
#     #run cv script
#     Rscript --vanilla scripts/cv-multi.R $par
#     
#  
#  
echo "-----------------------"
echo "DONE"
echo "-----------------------"
