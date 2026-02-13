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

#summary(lm(m.Body ~ loc.fix + PC1 + PC2 + PC3, data = gwas))
#summary(lm(l.Sperm ~ loc.fix + PC1 + PC2 + PC3, data = gwas))
#summary(lm(v.Spermatheca ~ loc.fix + PC1 + PC2 + PC3, data = gwas))


gwas$adj.m.Body = lm(m.Body ~ loc.fix + PC1, data = gwas)$residuals %>% 
  round(4)
  
  #TODO: v sperm should be beta distributed (or something?)
  
gwas$adj.v.Sperm = lm(v.Sperm ~ loc.fix + PC1, data = gwas)$residuals %>% 
  round(4)
gwas$adj.l.Sperm = lm(l.Sperm ~ loc.fix + PC1 + PC3, data = gwas)$residuals %>% 
  round(4)

#write out
gwas.out = data.frame(fid = gwas$id, 
                      iid = gwas$id, 
                      weight = gwas$adj.m.Body, 
                      vsperm = gwas$adj.v.Sperm,
                      lsperm = gwas$adj.l.Sperm)
                      
write.table(file = "data/qq_weight.pheno",
            gwas.out %>% select(fid, iid, weight),
            col.names = F, row.names = F, quote = F,
            sep = "\t")
write.table(file = "data/qq_vsperm.pheno",
            gwas.out %>% select(fid, iid, vsperm),
            col.names = F, row.names = F, quote = F,
            sep = "\t")
write.table(file = "data/qq_lsperm.pheno",
            gwas.out %>% select(fid, iid, lsperm),
            col.names = F, row.names = F, quote = F,
            sep = "\t")
  
  
  
  
#prepare files for BLUP
 blup_rename = function(v){
    u = unique(v[!is.na(v)])
    out = v
    for(i in 1:length(u)){
      out[v == u[i]] = i
    }
    return(as.numeric(out))
  }


preblup = pheno %>% 
  filter(id %in% pca.geno$id) 

preblup = preblup %>% 
  select(id, loc = loc.fix, 
  lsperm = l.Sperm, weight = m.Body, vsperm = v.Sperm,
  tsperm = t.Sperm) %>% 
  mutate(iid = 1:nrow(preblup),
         locid = blup_rename(loc),
         lsperm = round(scale(as.numeric(lsperm))[,1],4),
         weight = round(scale(as.numeric(weight))[,1],4),
         vsperm = round(scale(as.numeric(vsperm))[,1],4),
         tsperm = round(scale(as.numeric(tsperm))[,1],4)
         )

# sum(is.na(preblup$weight))
# sum(is.na(preblup$lsperm))
# 
# var(preblup$lsperm)
# var(preblup$weight)

#output for pheno
  blup = preblup %>% 
    select(iid, locid, lsperm, weight, vsperm, tsperm, id, loc)
  
  write.table(blup, "blup/pheno.txt", 
              col.names = F, row.names = F, quote = F)
              
#read grm
  
  G = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-filter.rel", 
  sep = "", header = F) %>% as.matrix()
  
  Gid = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-filter.rel.id", 
  sep = "", header = F)
  
  colnames(G) = rownames(G) = Gid[,1]
  
  sum(blup$id %in% colnames(G))
  
  #output relationship matrix
  
  final.mat = G[preblup$id, preblup$id]
  
  N = dim(final.mat)[1]
  covmat = matrix(ncol = 3, nrow = (N*N-N)/2 + N)
  cmr = 0
  for(i in 1:N){
    for(j in 1:i){
      cmr = cmr +1
      covmat[cmr,] = c(i,j, final.mat[i,j])
    }
  }
  
  write.table(covmat, "blup/covmat.txt", row.names = F, col.names = F, sep = " ")
  
max(covmat[,2])
  
sapply(blup, max)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
