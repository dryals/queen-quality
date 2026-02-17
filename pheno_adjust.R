library(dplyr)
library(readxl)

select=dplyr::select


#read phenotypic data
pheno = read_excel("data/phenotypes.xlsx") %>% 
  mutate(pheno_id = gsub(" ", "", QC))
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
pheno = pheno %>% left_join(loc.trans, by = "Location")

  #remove duplicate entries
  pheno = pheno[-(which(pheno$pheno_id == "QC2573")[2]),]
  pheno = pheno[-(which(pheno$pheno_id == "QC2422")[2]),]

#data cleaning and prep

    #pull all names from sequencer
    allnames = read.delim("data/all.names", header = F) %>% 
      rename(gc_id = V1) %>% 
      filter(grepl("QC", gc_id))
    
    #attempt bradley fixes
    bradley = read.csv("data/bradley_edits.csv")
    
    allnames = allnames %>% left_join(bradley %>% select(gc_id = manifest_id, new_id), by = 'gc_id')
      allnames$new_id[is.na(allnames$new_id)] = allnames$gc_id[is.na(allnames$new_id)] 
    
    # #how many fail to match?
    # allnames$gc_id[!allnames$gc_id %in% pheno$id]
    # allnames$new_id[!allnames$new_id %in% pheno$id]
      
      
    #manually drop duplicated genomic sample 
      allnames = allnames[-which(allnames$new_id == "QC2573")[1],]
  
      
    #amend pheno with genetic ids
      pheno = pheno %>% 
        left_join(allnames %>% 
                    select(pheno_id = new_id, gc_id), by = 'pheno_id')
      
      #sum(!is.na(pheno$gc_id))
      #nrow(allnames)
  
  
  #check for phenotype outliers
      pheno.num = pheno %>% 
        mutate(
          m.Body = as.numeric(m.Body),
          v.Sperm = as.numeric(v.Sperm),
          l.Sperm = as.numeric(l.Sperm),
          t.Sperm = as.numeric(t.Sperm),
          w.Head = as.numeric(w.Head),
          w.Thorax = as.numeric(w.Thorax),
          d.Spermatheca = as.numeric(d.Spermatheca),
          f.Spermatheca = as.numeric(f.Spermatheca)
        )
      
      # for(trait in colnames(pheno.num)[5:12]){
      #   
      #   hist(pheno.num[,trait] %>% unlist(), main = trait)
      #   
      # }
      
      #remove outlier phenotypes
      pheno.num = pheno.num[-(which.min(pheno.num$m.Body)),]
      
      #TODO: problems with thorax and head, ask bradley
      
  #collapse some small locs
    #TODO: try to fix these???
    small = pheno.num %>% group_by(loc.fix) %>% 
      summarise(n = n()) %>% 
      filter(n<5) %>% 
      ungroup() %>% 
      select(loc.fix)
    
    pheno.num$loc.fix[pheno.num$loc.fix %in% small$loc.fix] = "USA"
    
    #write cleaned phenotypes
    write.csv(pheno.num, "data/cleaned_pheno.csv",
              row.names = F, quote= F)
  

#read plink PCA
  pca.geno = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-pca.eigenvec",
  #pca.geno = read.delim("data/samples-pca.eigenvec",
                       header = F, sep = "")[,c(1, 3:5)]
    colnames(pca.geno) = c("gc_id", "PC1", "PC2", "PC3")

#prepare gwas  
gwas = pheno.num %>% 
  filter(gc_id %in% pca.geno$gc_id) %>% 
  left_join(pca.geno %>% select(gc_id, PC1, PC2, PC3), by = 'gc_id')

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


preblup = pheno.num %>% 
  filter(gc_id %in% pca.geno$gc_id) 

preblup = preblup %>% 
  select(gc_id, loc = loc.fix, 
  lsperm = l.Sperm, weight = m.Body, vsperm = v.Sperm,
  tsperm = t.Sperm) %>% 
  mutate(iid = 1:nrow(preblup),
         locid = blup_rename(loc),
         lsperm = round(scale(lsperm)[,1],4),
         weight = round(scale(weight)[,1],4),
         vsperm = round(scale(vsperm)[,1],4),
         tsperm = round(scale(tsperm)[,1],4)
         )

# sum(is.na(preblup$weight))
# sum(is.na(preblup$lsperm))
# 
# var(preblup$lsperm)
# var(preblup$weight)

#output for pheno
  blup = preblup %>% 
    select(iid, locid, lsperm, weight, vsperm, tsperm, gc_id, loc)
  
  write.table(blup, "blup/pheno.txt", 
              col.names = F, row.names = F, quote = F)
              
#read grm
  
  G = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-filter.rel", 
  sep = "", header = F) %>% as.matrix()
  
  Gid = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-filter.rel.id", 
  sep = "", header = F)
  
  colnames(G) = rownames(G) = Gid[,1]
  
  #sum(blup$id %in% colnames(G))
  
  #output relationship matrix
  
  final.mat = G[preblup$gc_id, preblup$gc_id]
  
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
