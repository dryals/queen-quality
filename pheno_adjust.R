library(dplyr)
library(readxl)

select=dplyr::select


#read phenotypic data
pheno = read_excel("data/phenotypes.xlsx") %>% 
  mutate(id = gsub(" ", "", QC))
#standardize locations
loc.trans = data.frame(Location = c("Hawaii", "Georgia", "Southern California", "Minnesota",
                                    "Northern California", "Washington", "West Virginia",
                                    "Michigan", "Unknown", "NCA", "NC", "USA", "MN", 
                                    "GA", "HI", "PA", "CA", "OH", "NY", "WA", "SCA",
                                    "VA", "AL", "FL", "OR"), 
                       
                       loc.fix = c("HI", "GA", "SCA", "MN",
                                   "NCA", "WA", "WV",
                                   "MI", "USA", "NCA", "NC", "USA", "MN", 
                                   "GA", "HI", "PA", "CA", "OH", "NY", "WA", "SCA",
                                   "VA", "AL", "FL", "OR"))
pheno = pheno %>% left_join(loc.trans)

#remove duplicate second entry
pheno = pheno[-(which(pheno$id == "QC2573")[2]),]

#collapse some small locs
  #TODO: try to fix these???
  small = pheno %>% group_by(loc.fix) %>% 
    summarise(n = n()) %>% 
    filter(n<5) %>% 
    ungroup() %>% 
    select(loc.fix)
  
  pheno$loc.fix[pheno$loc.fix %in% small$loc.fix] = "USA"
  

#read plink PCA
  pca.geno = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-pca.eigenvec",
                       header = F, sep = "")[,c(1, 3:5)]
    colnames(pca.geno) = c("id", "PC1", "PC2", "PC3")

#prepare gwas  
gwas = pheno %>% 
  filter(id %in% pca.geno$id) %>% 
  left_join(pca.geno %>% select(id, PC1, PC2, PC3))

sapply(gwas, function(x){sum(is.na(x))})

summary(lm(m.Body ~ loc.fix + PC1 + PC2 + PC3, data = gwas))
#PC3 is no longer significant

gwas$adj.m.Body = lm(m.Body ~ loc.fix + PC1 + PC2, data = gwas)$residuals %>% 
  round(4)
gwas$adj.v.Sperm = lm(v.Sperm ~ loc.fix + PC1 + PC2, data = gwas)$residuals %>% 
  round(4)

#write out
gwas.out = data.frame(fid = gwas$id, 
                      iid = gwas$id, 
                      weight = gwas$adj.m.Body, 
                      viability = gwas$adj.v.Sperm)
write.table(file = "data/qq_weight.pheno",
            gwas.out %>% select(-viability),
            col.names = F, row.names = F, quote = F,
            sep = "\t")
write.table(file = "data/qq_viability.pheno",
            gwas.out %>% select(-weight),
            col.names = F, row.names = F, quote = F,
            sep = "\t")
