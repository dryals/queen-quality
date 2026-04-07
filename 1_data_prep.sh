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
        #clumping
        #independent chromosomal segments for MHT
        
     
    #TODO: try removing missing data before calling bialleleic sites, might retain more that way!
    #TODO: why are there two LDprune scripts?
    
    #TODO: make sure GRM filter hapens once only, same samples used for all analyses
    
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
#     #output allele frequencies
#     cd $CLUSTER_SCRATCH/queen-quality/plink
#     plink -bfile samples-filter --freq --out samples-filter
#      
#     #output sites for ref filter, samples
#     awk '{print $2}' samples-filter.bim | tr ":" "\t" > samples-filter.sites
#     awk '{print $1}' samples-filter.fam > samples-filter.names
# 
# echo "-----------------------"
# echo "preparing phenotypic data in R..."
#     cd ~/ryals/queen-quality
#     R --vanilla --no-save --no-echo --silent < scripts/clean_pheno.R     
#  
echo "-----------------------"
    echo "ADMIXTURE analysis"
    cd ${CLUSTER_SCRATCH}/queen-quality
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
    echo "    filtering references ..."
        bcftools view reference.bcf.gz -T plink/samples-filter.sites --threads $SLURM_NTASKS \
        -Ob -o reference-filter.bcf.gz

     bcftools index -c reference-filter.bcf.gz
     
    echo "launching Ia script...."
        #count number of samples in each population
        cd ${CLUSTER_SCRATCH}/queen-quality
        mkdir -p aim
        cd aim
        #specify reference file
        echo "reference-filter.bcf.gz" > ref_filename.txt
        #reset logifle
        cd /home/dryals/ryals/queen-quality
        echo -n "" > outputs/aim.out

        #launch the Ia script array
        sbatch --array=1-16 scripts/AIM_v3.sh
        
    echo "waiting for Ia results (see aim.out)..."
    cd ~/ryals/queen-quality
    while [ $(grep "FINISHED" outputs/aim.out | wc -l | awk '{print $1}') -lt 16 ] #wait for all 16 to finish
    do
        sleep 20 #wait between each check
    done
    
    echo "compiling Ia results..."    
    cd ${CLUSTER_SCRATCH}/queen-quality/aim
    #this will hold all the aims
    cat chr*/chr*.ia | grep -v "chr" | sort -k3 -gr > aim.ia.txt
   
        grep -v "NA" aim.ia.txt | awk '$3>0' | awk 'OFS=":" {print$1, $2}' > plink_aim.ia.txt
          cat plink_aim.ia.txt | tr ":" "\t" > bct_aim.ia.txt
          
        
    count=$( wc -l plink_aim.ia.txt | awk '{print $1}')
    echo "    Calculated Ia for $count sites"
   
    echo "merging samples and references..."
    cd $CLUSTER_SCRATCH/queen-quality
    bcftools view samples.missing.bcf.gz -T plink/samples-filter.sites -S plink/samples-filter.names \
        --threads $SLURM_NTASKS -Ob -o samples-aim.bcf.gz
        
        bcftools index -c samples-aim.bcf.gz
        
    bcftools merge reference-filter.bcf.gz samples-aim.bcf.gz -m snps -Ou | \
        bcftools norm -m +snps -Ou | bcftools view -M2 -m2 --threads $SLURM_NTASKS -Ob -o qq-admix.bcf.gz
    echo "  indexing..."
        bcftools index -c qq-admix.bcf.gz
    
    echo "starting PLINK for supervised data" 
        #aims, maf, and geno
        cd ${CLUSTER_SCRATCH}/queen-quality
            echo "    pulling informative sites..."
            plink --bcf qq-admix.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
                --set-missing-var-ids @:# --silent \
                --extract aim/plink_aim.ia.txt--threads $SLURM_NTASKS --out plink/topaim

            #extract ld data, removing references
            echo "    calculating ld..."
            cd plink
            #TODO: consider maf threshold
            plink --bfile topaim -r2 --ld-window 1000 --ld-window-kb 50 --ld-window-r2 0.2 \
                --remove /home/dryals/ryals/admixPipeline/references/plink_pureRefs.txt \
                --silent --threads $SLURM_NTASKS --out sampleAncPrune
                
            #run R script to generate best set of sites by Ia
                #reset logfile
                cd /home/dryals/ryals/queen-quality
                echo -n "" > outputs/prune.out
                #start
                sbatch --array=1-16 scripts/prune_array.sh
                #wait
                echo "    waiting for pruning (see prune.out)..."
                while [ $(grep "FINISHED" outputs/prune.out | wc -l | awk '{print $1}') -lt 16 ] #wait for all 16 to finish
                do
                    sleep 20 #wait between each check
                done
                #create full output
                echo "    compiling results..."
                cd ${CLUSTER_SCRATCH}/queen-quality/aim
                #this will hold all the aims
                cat chr*/LDremove.txt > allLDremove.txt
                count=$( wc -l allLDremove.txt | awk '{print $1}')
                echo "    marked $count sites"
                
            echo "    removing pruned sites..."
            cd $CLUSTER_SCRATCH/queen-quality/plink
            #create new admix file
            plink --bfile topaim --make-bed --exclude ../aim/allLDremove.txt --silent \
                --threads $SLURM_NTASKS --out finaladmix
        
        #use this plink file basename for admix scripts
        echo "finaladmix" > plink_admix_filename.txt
        
        #kill script if the above fails
        if [ ! -f "finaladmix.bed" ]; then
            echo -e "\e[30;41m Plink Failed! \e[0m"
            exit 1
        fi
        
echo "starting supervised admix..."  
    cd $CLUSTER_SCRATCH/queen-quality
    mkdir -p admix
    cd admix 
    mkdir -p supervised
    #supervised
    cd /home/dryals/ryals/queen-quality
    #create pop file
    R --vanilla --no-save --no-echo --silent < scripts/makeAdmixPop.R
    sbatch scripts/supervised_admix_v3.sh

# echo "-----------------------"
#     #PCA and GRM
#     cd $CLUSTER_SCRATCH/queen-quality/plink
#     echo "PCA..."
#     plink --bfile samples-filter --pca 500 \
#         --threads $SLURM_NTASKS --out samples-gwas --silent
#         
#     echo "    select locations..."
#     for LOC in GA CA HI
#         do
#             echo "        working $LOC ..."
#             cd $CLUSTER_SCRATCH/queen-quality/plink
#             plink --bfile samples-filter --pca 500 \
#                 --make-bed \
#                 --keep ~/ryals/queen-quality/data/phenotyped_${LOC}.gcnames \
#                 --threads $SLURM_NTASKS --out samples-gwas_${LOC} --silent
#         done
#      
#     echo "GRM..."
#     module purge
#     module load biocontainers plink2
#     
#     plink2 --bfile samples-filter --keep ~/ryals/queen-quality/data/phenotyped.plink \
#         -maf 0.05 \
#         -make-rel square --out samples-gs
# 
#     module purge
#     module load biocontainers bcftools vcftools plink r
#     
# #     echo "GRM in GCTA..."
# #     cd ~/ryals/queen-quality/data
# #     paste phenotyped.gcnames phenotyped.gcnames > phenotyped.plink
# #     
# #     cd $CLUSTER_SCRATCH/queen-quality/plink
# #     
# #     plink --bfile samples-filter --maf 0.01 \
# #         --keep ~/ryals/queen-quality/data/phenotyped.plink \
# #         --make-bed --threads $SLURM_NTASKS --out samples-gs
# #         
# #     gcta=/depot/bharpur/apps/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
# #     
# #         $gcta --bfile samples-gs --make-grm --thread-num $SLURM_NTASKS \
# #             --autosome-num 16 --out samples-gs
#                          
# echo "-----------------------"
# echo "preparing data for GWAS and GS"
#     cd ~/ryals/queen-quality
#     R --vanilla --no-save --no-echo --silent < scripts/adjust_pheno.R
# 
# echo "-----------------------"
# echo "running GWAS..."
# 
# #TODO pull GA, CA, HI and run separately
# 
#     echo "    GRM..."
#     cd $CLUSTER_SCRATCH/queen-quality
#     mkdir -p gwas
#     cd gwas
#     gcta=/depot/bharpur/apps/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
#     
#         $gcta --bfile ../plink/samples-filter --make-grm --thread-num $SLURM_NTASKS \
#             --autosome-num 16 --out qq
#             
#         echo "    mlma: all locations..."
#         #adjusted phenotypes generated in R script...
#         $gcta --mlma-loco --bfile ../plink/samples-filter --grm qq \
#             --pheno ~/ryals/queen-quality/data/qq_vsperm.pheno \
#             --autosome-num 16 \
#             --out qq_vsperm --thread-num $SLURM_NTASKS &> /dev/null
#             
#         $gcta --mlma-loco --bfile ../plink/samples-filter --grm qq \
#             --pheno ~/ryals/queen-quality/data/qq_weight.pheno \
#             --autosome-num 16 \
#             --out qq_weight --thread-num $SLURM_NTASKS &> /dev/null
#             
#        echo "    mlma: select locations..."
#             for LOC in GA CA HI
#             do
#                 echo "    working: ${LOC}..."
#                 
#                 $gcta --bfile ../plink/samples-gwas_${LOC} --make-grm --thread-num $SLURM_NTASKS \
#                     --autosome-num 16 --out qq_${LOC} &> /dev/null
#                 
#                 $gcta --mlma-loco --bfile ../plink/samples-gwas_${LOC} --grm qq_${LOC} \
#                     --pheno ~/ryals/queen-quality/data/qq_vsperm_${LOC}.pheno \
#                     --autosome-num 16 \
#                     --out qq_vsperm_${LOC} --thread-num $SLURM_NTASKS &> /dev/null
#                     
#                 $gcta --mlma-loco --bfile ../plink/samples-gwas_${LOC} --grm qq_${LOC} \
#                     --pheno ~/ryals/queen-quality/data/qq_weight_${LOC}.pheno \
#                     --autosome-num 16 \
#                     --out qq_weight_${LOC} --thread-num $SLURM_NTASKS &> /dev/null
#                     
#                 echo "        done: ${LOC}"
# 
#             done
#             
#         cp qq_*.loco.mlma ~/ryals/queen-quality/data
#             
# echo "-----------------------"
# echo "running GREML..."
#     #ld scores
#     $gcta --bfile ../plink/samples-filter --ld-score-region 200 --ld-wind 50 --ld-rsq-cutoff 0.1 --out greml
#     
#     
#     R
#     lds_seg = read.table("greml.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
#     quartiles=summary(lds_seg$ldscore_SNP)
# 
#     lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
#     lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
#     lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
#     lb4 = which(lds_seg$ldscore_SNP > quartiles[5])
# 
#     lb1_snp = lds_seg$SNP[lb1]
#     lb2_snp = lds_seg$SNP[lb2]
#     lb3_snp = lds_seg$SNP[lb3]
#     lb4_snp = lds_seg$SNP[lb4]
# 
#     write.table(lb1_snp, "snp_group1.txt", row.names=F, quote=F, col.names=F)
#     write.table(lb2_snp, "snp_group2.txt", row.names=F, quote=F, col.names=F)
#     write.table(lb3_snp, "snp_group3.txt", row.names=F, quote=F, col.names=F)
#     write.table(lb4_snp, "snp_group4.txt", row.names=F, quote=F, col.names=F)
#     quit("no")
#     
#     #multi GRM
#     $gcta --bfile ../plink/samples-filter --extract snp_group1.txt --make-grm --out greml_group1
#     $gcta --bfile ../plink/samples-filter --extract snp_group2.txt --make-grm --out greml_group2
#     $gcta --bfile ../plink/samples-filter --extract snp_group3.txt --make-grm --out greml_group3
#     $gcta --bfile ../plink/samples-filter --extract snp_group4.txt --make-grm --out greml_group4
#     
#     echo "greml_group1" > multi_grm.txt
#     echo "greml_group2" >> multi_grm.txt
#     echo "greml_group3" >> multi_grm.txt
#     echo "greml_group4" >> multi_grm.txt
#     
#     #GREML
#     for trait in weight vsperm
#     do
#         $gcta --reml --mgrm multi_grm.txt --pheno ~/ryals/queen-quality/data/qq_${trait}.pheno \
#         --reml-pred-rand --out greml_${trait}
#         
#         $gcta --bfile ../plink/samples-filter --blup-snp greml_${trait}.indi.blp --out greml_${trait}
#     done
#     
#     cp greml_*.snp.blp ~/ryals/queen-quality/data
#     
# # echo "-----------------------"  
# # echo "running BLUP..."
# # 
#     par=wv
# 
#     #TODO: single-trait blups
#     cd ~/ryals/queen-quality/blup
#         #create links
#         if [ ! -f  blupf90+ ]; then
#             echo "    creating links..."
#             ln -S blupf90+ /depot/bharpur/apps/blupf90/blupf90+
#             ln -S airemlf90 /depot/bharpur/apps/blupf90/airemlf90 
#             ln -S validationf90 /depot/bharpur/apps/blupf90/validationf90
#         fi
# 
#     cd ~/ryals/queen-quality
#     cp params/${par}.par0 blup
#     cd blup
#     ./airemlf90 ${par}.par0
#     
# #     
# #     
# #     #TODO this but automatic...
# #     #read G and R matricies into blup.par2
# # #     sed -n 16,80p file1>patch
# # #     sed -i 18rpatch file2
# #     
#     cp ../params/${par}.par1 .
#     ./blupf90+ ${par}.par1
#     cp solutions ../data/sol-${par}.txt
#     cp solutions sol-${par}
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

echo "-----------------------"
echo "DONE"
date
echo "-----------------------"
