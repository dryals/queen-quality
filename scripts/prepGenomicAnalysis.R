library(dplyr)
library(readxl)
library(lubridate)
library(Matrix)

select=dplyr::select

#read phenotypic data
pheno = read.csv("data/cleaned_pheno.csv")
                
#read plink PCA
  pca.geno = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-gwas.eigenvec",
  #pca.geno = read.delim("data/samples-pca.eigenvec",
                       header = F, sep = "")[,c(1, 3:5)]
    colnames(pca.geno) = c("gc_id", "PC1", "PC2", "PC3")

#prepare gwas  
gwas = pheno %>% 
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
  
  G.p = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-gs2.rel", 
  sep = "", header = F) %>% as.matrix()
  
  G.p.id = read.delim("/scratch/negishi/dryals/queen-quality/plink/samples-gs2.rel.id", 
  sep = "", header = T)
  
  colnames(G.p) = rownames(G.p) = G.p.id[,1]
# 
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
#         G = ReadGRMBin("/scratch/negishi/dryals/queen-quality/plink/samples-gs")
#         
#         #TODO: what is the order of G$off???
#         
#         G.mat = matrix(nrow = nrow(G$id), ncol = nrow(G$id))
#           colnames(G.mat) = rownames(G.mat) = G$id[,1]
#           
#         G.mat[lower.tri(G.mat, diag = F)] = round(G$off,4)
#         
#         G.mat[1:10,1:10]
#         head(G$off)
#         
#         G.mat = forceSymmetric(G.mat, uplo = "L") %>% as.matrix()
#         
#         diag(G.mat) = round(G$diag,4)
        
    #remove high diags
    #remove = G$id[G$diag > 1.8,1]
    remove = colnames(G.p)[diag(G.p) > 1.7]

#         #compare the grms
#         G.test = G.mat[lower.tri(G.mat, diag = F)] - 
#                  G.p[lower.tri(G.p, diag = F)]
#         
#         G.testdiag = diag(G.mat) - diag(G.p)
#         
#         range(diag(G.mat))
#         range(diag(G.p))
#                  
#         pdf(file = "grmcompare.pdf")
#           hist(G.test)
#         dev.off()
#         pdf(file = "grmcomparediag.pdf")
#           hist(G.testdiag)
#         dev.off()
#         
#         sum(colnames(G.p) == colnames(G.mat))

  
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
  left_join(pheno)

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
  
#   #remove outliers by diag values
#     print("remove")
#     print(remove)
#     
#   #output for pheno
  blup = preblup %>% 
    filter(!gc_id %in% remove)

    blup = blup %>%
    mutate(iid = 1:nrow(blup)) %>%
    select(iid, locid, lsperm, weight, vsperm, tsperm, loc, gc_id, pheno_id)
  
  write.table(blup, "blup/pheno.txt", 
              col.names = F, row.names = F, quote = F)
    
#     
#   
#   #output relationship matrix
  
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
#   
#   covmat.old = covmat
#   #another version for gcta grm
#     N = dim(G$id)[1]
#     covmat = matrix(ncol = 3, nrow = (N*N-N)/2 + N)
#     cmr = 0
#     k = 0
#     for(i in 1:N){
#       for(j in 1:i){
#         cmr = cmr +1
#         if(i == j){
#             covmat[cmr,] = c(i,j, G$diag[i])
#         } else {
#             k = k+1
#             covmat[cmr,] = c(i,j, G$off[k])
#         }
#         
#       }
#     }
#   
#   
  write.table(covmat, "blup/covmat.txt", row.names = F, col.names = F, sep = " ")
  
max(covmat[,2])
  
sapply(blup, max)
