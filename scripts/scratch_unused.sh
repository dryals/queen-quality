    echo "heterozygosity QC"
    echo "    formatting..."
    cd $CLUSTER_SCRATCH/queen-quality
    
    bcftools view samples.missing.bcf.gz -T plink/samples-filter.sites -S plink/samples-filter.names \
        --threads $SLURM_NTASKS -Ob -o samples-filter.bcf.gz
        
        bcftools index -c samples-filter.bcf.gz
        
        #convert to vcf
        bcftools view samples-filter.bcf.gz --threads $SLURM_NTASKS -Oz -o samples-filter.vcf.gz 
    echo "    calculating het..."
    vcftools --gzvcf samples-filter.vcf.gz --het --out samples-filter
    #sort by F
    cat *.het | sort -k5 -gr > sorted.het
    #list outliers in R
    R
        sorted = read.delim("sorted.het", header = F, sep = "")[,c(1,5)]
            colnames(sorted) = c("gc_id", "F")
        #remove last row (old header)
        sorted = sorted[-nrow(sorted),]
        #standardize
        sorted$Fz = scale(as.numeric(sorted$F))
        #visualize
        png(file = "hethist.png")
            hist(sorted$Fz)
        dev.off()
        
        #remove anything > 3 std dev 
        remove = sorted$gc_id[abs(sorted$Fz) > 3]
        
        nrow(remove)
        
        write.table(remove, file = "het.remove",
                    row.names = F, quote = F, col.names = F)
        
        
    quit(save = "no")
    
    #format for plink
    paste het.remove het.remove > het.remove.plink
    
