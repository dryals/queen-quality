library(dplyr)
library(readxl)
library(lubridate)
library(Matrix)

select=dplyr::select


#read phenotypic data
pheno = read_excel("data/phenotypes.xlsx") %>% 
  mutate(pheno_id = gsub(" ", "", QC)) 

#repair dates
  pheno$date.1 = as_date(NA)
  pheno$date.1[grepl("^[0-9]*$", pheno$Received)] = 
    as_date(as.numeric(pheno$Received[grepl("^[0-9]*$", pheno$Received)]), origin = "1899-12-30")
  
  pheno$date.2 = as_date(NA)
  pheno$date.2[grepl("Z$", pheno$Received)] = 
    as_date(pheno$Received[grepl("Z$", pheno$Received)])
  
  
  pheno$date = as_date(NA)
  for(i in 1:nrow(pheno)){
    if(is.na(pheno$date.1[i])){
      pheno$date[i] = pheno$date.2[i]
    } else{
      pheno$date[i] = pheno$date.1[i]
    }
  }
  
  pheno = pheno %>% select(-date.1, -date.2)
  
  pheno$year = year(pheno$date)
  


#standardize locations
loc.trans = data.frame(Location = c("Hawaii", "Georgia", "Southern California", "Minnesota",
                                    "Northern California", "Washington", "West Virginia",
                                    "Michigan", "Unknown", "NCA", "NC", "USA", "MN", 
                                    "GA", "HI", "PA", "CA", "OH", "NY", "WA", "SCA",
                                    "VA", "AL", "FL", "OR", "Oregon"), 
                       
                       loc.fix = c("HI", "GA", "SCA", "MN",
                                   "NCA", "WA", "WV",
                                   "MI", "USA", "NCA", "NC", "USA", "MN", 
                                   "GA", "HI", "PA", "CA", "OH", "NY", "WA", "SCA",
                                   "VA", "AL", "FL", "OR", "OR"))
pheno = pheno %>% left_join(loc.trans, by = "Location")

  #remove duplicate entries
  pheno = pheno[-(which(pheno$pheno_id == "QC2573")[2]),]
  pheno = pheno[-(which(pheno$pheno_id == "QC2422")[2]),]
  
  #print(pheno[,c('Location', 'loc.fix')], n = 500)

#data cleaning and prep

    #pull all names from sequencer
    
    allnames = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-filter.fam", header = F, sep = "") %>% 
    #allnames = read.delim("data/samples-filter.fam", header = F, sep = "") %>% 
      select(gc_id = V1) %>% 
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
      pheno.fix = pheno %>% 
        left_join(allnames %>% 
                    select(pheno_id = new_id, gc_id), by = 'pheno_id')
                    
    #manually drop duplicated pheno id
    pheno.fix = pheno.fix[-which(pheno.fix$pheno_id == "QC2422")[2],]
    
      
      #sum(!is.na(pheno$gc_id))
      #nrow(allnames)
      
#     pheno.fix %>% group_by(gc_id) %>%
#       summarise(n = n()) %>%
#       filter(n > 1) 
#     
#     pheno.fix %>% group_by(pheno_id) %>%
#       summarise(n = n()) %>%
#       filter(n > 1) 
#     
  
  #check for phenotype outliers
      pheno.num = pheno.fix %>% 
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
      
    #collapse some small levels
      #manually fix northern california 2019
      pheno.num$loc.fix[pheno.num$loc.fix == "NCA" &
                          pheno.num$year == 2019] = "CA"
      
      #loc.years with small count become USA.year
      small = pheno.num %>% group_by(loc.fix, year) %>% 
        summarise(n = n()) %>% 
        filter(n<5) %>% 
        ungroup() 
      
      
      
      for(i in 1:nrow(small)){
        pheno.num$loc.fix[pheno.num$loc.fix == small$loc.fix[i] &
                          pheno.num$year == small$year[i]] = "USA"
      }
    #create loc.year effect
    pheno.num$loc.year = paste0(pheno.num$loc.fix, pheno.num$year)
    
    #collapse remaining small USA levels
    small.usa = pheno.num %>% filter(loc.fix == "USA") %>% 
      group_by(loc.year) %>% 
      summarise(n = n()) %>% 
      filter(n<5) %>% 
      ungroup
    
    pheno.num$loc.year[pheno.num$loc.year %in% small.usa$loc.year] = "USA20XX"
    
    #table(pheno.num$loc.year)

    
    #write cleaned phenotypes
    write.csv(pheno.num, "data/cleaned_pheno.csv",
              row.names = F, quote= F)
  

#read plink PCA
  pca.geno = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-gwas.eigenvec",
  #pca.geno = read.delim("data/samples-pca.eigenvec",
                       header = F, sep = "")[,c(1, 3:5)]
    colnames(pca.geno) = c("gc_id", "PC1", "PC2", "PC3")

#prepare gwas  
gwas = pheno.num %>% 
  filter(gc_id %in% pca.geno$gc_id) %>% 
  left_join(pca.geno %>% select(gc_id, PC1, PC2, PC3), by = 'gc_id')

#sapply(gwas, function(x){sum(is.na(x))})

#summary(lm(m.Body ~ loc.year + PC1 + PC2, data = gwas))
#summary(lm(v.Sperm ~ loc.year + PC1 + PC2, data = gwas))


gwas$adj.m.Body = lm(m.Body ~ loc.year + PC1 + PC2, data = gwas)$residuals %>% 
  round(4)
gwas$adj.v.Sperm = lm(v.Sperm ~ loc.year + PC1 + PC2, data = gwas)$residuals %>% 
  round(4)

# gwas$adj.l.Sperm = lm(l.Sperm ~ loc.year + PC1 + PC2, data = gwas)$residuals %>% 
#   round(4)

#write out
gwas.out = data.frame(fid = gwas$gc_id, 
                      iid = gwas$gc_id, 
                      weight = gwas$adj.m.Body, 
                      vsperm = gwas$adj.v.Sperm)#,
                      #lsperm = gwas$adj.l.Sperm)
                      
write.table(file = "data/qq_weight.pheno",
            gwas.out %>% select(fid, iid, weight),
            col.names = F, row.names = F, quote = F,
            sep = "\t")
write.table(file = "data/qq_vsperm.pheno",
            gwas.out %>% select(fid, iid, vsperm),
            col.names = F, row.names = F, quote = F,
            sep = "\t")
# write.table(file = "data/qq_lsperm.pheno",
#             gwas.out %>% select(fid, iid, lsperm),
#             col.names = F, row.names = F, quote = F,
#             sep = "\t")
  

#TODO: try using plink's GRM again but with maf 0.05 ...
#TODO: why doesnt gcta's grm work?!

  #read grm
  
#   G = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-filter.rel", 
#   sep = "", header = F) %>% as.matrix()
#   
#   Gid = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-filter.rel.id", 
#   sep = "", header = T)
#   
#   colnames(G) = rownames(G) = Gid[,1]

        ReadGRMBin=function(prefix, AllN=F, size=4){
            sum_i=function(i){
                return(sum(1:i))
            }
                BinFileName=paste(prefix,".grm.bin",sep="")
                NFileName=paste(prefix,".grm.N.bin",sep="")
                IDFileName=paste(prefix,".grm.id",sep="")
                id = read.table(IDFileName)
                n=dim(id)[1]
                BinFile=file(BinFileName, "rb");
                grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
                NFile=file(NFileName, "rb");
                if(AllN==T){
                    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
                }
                else N=readBin(NFile, n=1, what=numeric(0), size=size)
                i=sapply(1:n, sum_i)
                return(list(diag=grm[i], off=grm[-i], id=id, N=N))
            }
        
        G = ReadGRMBin("/scratch/negishi/dryals/queen-quality/plink/samples-gs")
    
        remove = G$id[G$diag > 1.8,1]
        
        G.mat = matrix(nrow = nrow(G$id), ncol = nrow(G$id))
          colnames(G.mat) = rownames(G.mat) = G$id[,1]
          
        G.mat[lower.tri(G.mat, diag = F)] = round(G$off,4)
        diag(G.mat) = round(G$diag,4)
        
        G.mat = forceSymmetric(G.mat, uplo = "L") %>% as.matrix()


  
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
  filter(gc_id %in%  colnames(G.mat)) 

preblup = preblup %>% 
  select(gc_id, pheno_id, loc = loc.year, 
  lsperm = l.Sperm, weight = m.Body, vsperm = v.Sperm,
  tsperm = t.Sperm) %>% 
  mutate(
         locid = blup_rename(loc),
         lsperm = round(scale(lsperm)[,1],4),
         weight = round(scale(weight)[,1],4),
         vsperm = round(scale(vsperm)[,1],4),
         tsperm = round(scale(tsperm)[,1],4)
         )

#error checking
#   sapply(preblup, function(x){sum(is.na(x))})
# 
#   length(unique(preblup$gc_id))
#   length(unique(preblup$pheno_id))
#   nrow(preblup)
# 
#   preblup %>% group_by(pheno_id) %>%
#     summarise(n = n()) %>%
#     filter(n > 1) 
#     
         
# sum(is.na(preblup$weight))
# sum(is.na(preblup$lsperm))
# 
# var(preblup$lsperm)
# var(preblup$weight)
            
  
  #sum(blup$id %in% colnames(G))
  
  #remove outliers by diag values
    print("remove")
    print(remove)
    
  #output for pheno
  blup = preblup %>% 
    filter(!gc_id %in% remove)
    
    blup = blup %>%
    mutate(iid = 1:nrow(blup)) %>%
    select(iid, locid, lsperm, weight, vsperm, tsperm, loc, gc_id, pheno_id)
  
  write.table(blup, "blup/pheno.txt", 
              col.names = F, row.names = F, quote = F)
    
    
  
  #output relationship matrix
  
  final.mat = G.mat[blup$gc_id, blup$gc_id]
  
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
