#for parallel processing
#last edit   09 FEB 2026

args = commandArgs(trailingOnly=TRUE)
chrno = args[1]

suppressMessages(library(tidyverse))

#read data
setwd('/scratch/negishi/dryals/queen-quality')
ld = read.delim('plink/samples-preprune.ld', sep = '', header = T) %>% 
  filter(CHR_A == chrno, CHR_B == chrno)
maf = read.delim('plink/samples-preprune.frq', sep = '', header = T) %>% 
  filter(CHR == chrno)

#attach maf for A
ld = ld %>% left_join(maf %>%  select(MAF_A = MAF, CHR_A = CHR, SNP_A = SNP), by = c("CHR_A", "SNP_A"))
#... for B
ld = ld %>% left_join(maf %>%  select(MAF_B = MAF, CHR_B = CHR, SNP_B = SNP), by = c("CHR_B", "SNP_B"))
#difference 
ld$diff = ld$MAF_A - ld$MAF_B
  ld$MMAF = (ld$diff > 0) * ld$MAF_A + (ld$diff < 0) * ld$MAF_B + (ld$diff == 0) * ld$MAF_A 

ld2 = ld %>%
  arrange(desc(MMAF)) %>% 
  select(SNP_A, SNP_B, MAF_A, MAF_B, diff, MMAF)
toremove = vector("character", 1)

#remove site with lower MAF
while(nrow(ld2)>0){
  
  #keep snp with best maf
  if(ld2$diff[1] > 0){
    maxsnp = ld2$SNP_A[1]
  } else{
    maxsnp = ld2$SNP_B[1]
    }
    #identify all snps in ld with this one
    tmpremove = c(ld2$SNP_A[ld2$SNP_B == maxsnp], ld2$SNP_B[ld2$SNP_A == maxsnp])
    #remove all pairs containing those snps
    ld2 = ld2[!((ld2$SNP_A %in% tmpremove) | (ld2$SNP_B %in% tmpremove)),]
    #add to remove list
    toremove = c(toremove, tmpremove)
    rm(tmpremove)
}
    
#length(toremove)/nrow(maf)
  

#write out

write.table(toremove[-1], paste0("ld/chr", chrno, "/MAFremove.txt"), quote = F, row.names = F, col.names = F)

print(paste0("chr ", chrno, ": removing ", length(toremove), " sites"))

