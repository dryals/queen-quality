library(dplyr)
library(readxl)
library(lubridate)
library(Matrix)
library(matrixcalc)

select=dplyr::select

#read phenotypic data
pheno = read.csv("data/cleaned_pheno.csv")

phenotyped = read.csv("data/phenotyped.gcnames", header = F)
  colnames(phenotyped) = "gc_id"

                
#read plink PCA
  pca.geno = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-gwas.eigenvec",
  #pca.geno = read.delim("data/samples-pca.eigenvec",
                       header = F, sep = "")[,c(1, 3:5)]
    colnames(pca.geno) = c("gc_id", "PC1", "PC2", "PC3")

#retain hjust samples which have genomic data
pheno.filter = pheno %>%
filter(gc_id %in% phenotyped$gc_id)

    #collapse small levels, preserving year over location
      
      #loc.years with small count become USA.year
      small = pheno.filter %>% group_by(loc.fix, year) %>% 
        summarise(n = n()) %>% 
        filter(n<5) %>% 
        ungroup() 
      
      #consider site-years with small n as "USA"
      for(i in 1:nrow(small)){
        pheno.filter$loc.fix[pheno.filter$loc.fix == small$loc.fix[i] &
                          pheno.filter$year == small$year[i]] = "USA"
      }

    
    #overwrite loc.year
    pheno.filter$loc.year = paste0(pheno.filter$loc.fix, pheno.filter$year)
    
    #collapse remaining small USA levels
    small.usa = pheno.filter %>% filter(loc.fix == "USA") %>% 
      group_by(loc.year) %>% 
      summarise(n = n()) %>% 
      filter(n<5) %>% 
      ungroup
    
    pheno.filter$loc.year[pheno.filter$loc.year %in% small.usa$loc.year] = "USA20XX"
    
    table(pheno.filter$loc.year)
    
    
# #read admix components 
#     #read admix data
#     
#     admix = cbind(read.delim("/scratch/negishi/dryals/queen-quality/plink/finaladmix.fam", header = F, sep = "")[,1],
#                   read.delim("/scratch/negishi/dryals/queen-quality/admix/supervised/finaladmix.4.Q", header = F, sep = "")
#                   )
#       colnames(admix)[1]="gc_id"
#       
#     reffam = read.delim("references/refData.txt")
#     
#     #verify lineage identity
#     lins = admix %>% left_join(reffam %>% select(gc_id = SRR, lineage)) %>% 
#       filter(grepl ("SRR", gc_id)) %>%
#       pivot_longer(starts_with("V")) %>%
#       group_by(name, lineage) %>%
#       mutate(m = mean(value)) %>% ungroup %>% 
#       group_by(name) %>% arrange(desc(m)) %>% slice(1) %>% select(lineage, name)
#     
#     #rename
#     admix = admix %>% rename(A = V4, M = V1, C = V2, O = V3) %>% 
#       left_join(pheno %>% select(gc_id, loc.year))
    

#prepare gwas  
gwas = pheno.filter %>% 
  left_join(pca.geno %>% select(gc_id, PC1, PC2, PC3), by = 'gc_id') #%>%
  #left_join(admix %>% select(gc_id, A, M, C, O), by = 'gc_id')

#sapply(gwas, function(x){sum(is.na(x))})

#summary(lm(v.Sperm ~ loc.year + PC1 + PC2, data = gwas))


gwas$adj.m.Body = lm(m.Body ~ loc.year + PC1 + PC2 + PC3, data = gwas)$residuals %>% 
  round(4)
gwas$adj.v.Sperm = lm(v.Sperm ~ loc.year + PC1 + PC2 + PC3, data = gwas)$residuals %>% 
  round(4)
  
#write out whole data
  gwas.out = data.frame(fid = gwas$gc_id, 
                        iid = gwas$gc_id, 
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

# #adjust within select locations
#   for (LOC in c("HI", "CA", "GA")){
# 
#     #correct samples
#     phenotyped.loc = read.csv(paste0("data/phenotyped_", LOC, ".gcnames"),
#                               header = F, sep = "")
#       colnames(phenotyped.loc) = c("gc_id", "fam_id")
#     #read correct PCA
#     pca.loc = read.delim(paste0("/scratch/negishi/dryals/queen-quality/plink/samples-gwas_", 
#                           LOC, ".eigenvec"),
#                        header = F, sep = "")[,c(1, 3:5)]
#     colnames(pca.loc) = c("gc_id", "PC1", "PC2", "PC3")
#     #combine
#     gwas.loc = pheno %>% 
#     filter(gc_id %in% phenotyped.loc$gc_id) %>%
#     left_join(pca.loc %>% select(gc_id, PC1, PC2, PC3), by = 'gc_id')
#     #pool small year levels
#     smallyears = gwas.loc %>%
#     group_by(year) %>%
#     summarise(n = n()) %>%
#     filter(n<5)
#     
#     gwas.loc$year.fix = gwas.loc$year
#       gwas.loc$year.fix[gwas.loc$year.fix %in% smallyears$year] = "20XX"
#       
#     #adjust
#     gwas.loc$adj.m.Body = lm(m.Body ~ year.fix + PC1 + PC2 + PC3, data = gwas.loc)$residuals %>% 
#       round(4)
#     gwas.loc$adj.v.Sperm = lm(v.Sperm ~ year.fix + PC1 + PC2 + PC3, data = gwas.loc)$residuals %>% 
#       round(4)
#       
#     #write out
#     gwas.loc.out = data.frame(fid = gwas.loc$gc_id, 
#                       iid = gwas.loc$gc_id, 
#                       weight = gwas.loc$adj.m.Body, 
#                       vsperm = gwas.loc$adj.v.Sperm)
#                       
#     write.table(file = paste0("data/qq_weight_", LOC, ".pheno"),
#                 gwas.loc.out %>% select(fid, iid, weight),
#                 col.names = F, row.names = F, quote = F,
#                 sep = "\t")
#                 
#     write.table(file = paste0("data/qq_vsperm_", LOC, ".pheno"),
#                 gwas.loc.out %>% select(fid, iid, vsperm),
#                 col.names = F, row.names = F, quote = F,
#                 sep = "\t")
#     
#  }  
  

#Prepare GS
#read grm
  
  G.p = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-gs.rel", 
  sep = "", header = F) %>% as.matrix()
  
  G.p.id = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-gs.rel.id", 
  sep = "", header = T)
  
  colnames(G.p) = rownames(G.p) = G.p.id[,1]
  
  #force positive definite
  make.positive.definite <- function(M)
    {
      tol=1e-6
      eig <- eigen(M, symmetric=TRUE)
      rtol <- tol * eig$values[1]
      if(min(eig$values) < rtol)
      {
        vals <- eig$values
        vals[vals < rtol] <- rtol
        srev <- eig$vectors %*% (vals * t(eig$vectors))
        dimnames(srev) <- dimnames(M)
        return(srev)
      } else
      {
        return(M)
      }
    }
    
    G.p.tmp = make.positive.definite(G.p)
    G.p.tmp = round(G.p.tmp, 4)
    
    #is.positive.definite(G.p.tmp)
    
    #overwrite
    G.p = G.p.tmp
    
#prepare files for BLUP
 blup_rename = function(v){
    u = unique(v[!is.na(v)])
    out = v
    for(i in 1:length(u)){
      out[v == u[i]] = i
    }
    return(as.numeric(out))
  }

  #ensure same individuals in same order
preblup = data.frame(gc_id = colnames(G.p)) %>%
  left_join(pheno.filter) %>%
  #join pc's 
  left_join(gwas %>% select(gc_id, PC1, PC2, PC3))
  
preblup = preblup %>% 
  select(gc_id, pheno_id, loc = loc.year, PC1, PC2, PC3,
  lsperm = l.Sperm, weight = m.Body, vsperm = v.Sperm,
  tsperm = t.Sperm) %>% 

  
  #TODO: this but earlier in the pipeline!!
  #output phenotypes in correct units
  
    mutate(
         locid = blup_rename(loc),
         PC1 = round(PC1, 4),
         PC2 = round(PC2,4),
         lsperm = round(lsperm / 1e5 ,4),
         weight = round(weight,4),
         vsperm = round(vsperm * 100,4),
         tsperm = round(tsperm,4)
         )
         
  preblup[is.na(preblup)] = -999
  
      
#output phenotypes
  blup = preblup
  
    blup = blup %>%
    mutate(iid = 1:nrow(blup)) %>%
    select(iid, locid, PC1, PC2, PC3, lsperm, weight, vsperm, tsperm, loc, gc_id, pheno_id)
  
  write.table(blup, "blup/pheno.txt", 
              col.names = F, row.names = F, quote = F)
    
#output relationship matrix
  
  final.mat = G.p[blup$gc_id, blup$gc_id]
  
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
  
#print number of samples
print("n = ")
max(covmat[,2])
  
#print colnames and max values
sapply(blup, max)
