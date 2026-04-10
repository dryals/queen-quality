#take a command line argument
args = commandArgs(trailingOnly=TRUE)
#args = "wv"

library(dplyr)
library(tidyr)

#split into phenotypes under consideration
args.split = unlist(strsplit(args, ""))

ntrait = length(args.split)

#target script for CV runs
targetParam = paste0(args, "-cv.par1")

# print(targetParam)
# print(getwd())
#wd is queen quality!

#TODO: edit for multiple fixed effs and new pheno file!


# fixed.key = data.frame(effect = c(2:4) %>% as.character(),
#                        name = c("apnumid", "start", "requeen"))
# 
#key for all ids
trait.key = data.frame(tn = c("l", "w", "v", "t"),
                       trait = c(1:4) %>%  as.character())

#key for traits in this run
trait.key.used = data.frame( tn = args.split,
                            trait = c(1:ntrait) %>%  as.character())

#pull fixed effects from full run
fixed.eff = read.delim(paste0("data/sol-", args[1], ".txt"), sep = "", header = F)[-1,]
colnames(fixed.eff) = c("trait", "effect", "level" ,"solution", "se")
fixed.eff = fixed.eff %>% mutate(se = as.numeric(se),
                                 solution = as.numeric(solution)) %>% 
  filter(effect != 1)
  
#load sample info from pheno file
pheno = read.delim("blup/pheno.txt", header = F, sep = " ")
  colnames(pheno) = c("iid", "locid", "PC1","PC2","PC3", trait.key$tn, "loc", "gc_id", "pheno_id" )

#create list for CV
  #only mask individuals with full pheno data
  masked = pheno
  
  set.seed(2026)
  masked$CV = sample(c(1:5), nrow(pheno), replace = T)
  
  
  #table(masked$CV)

#create output object
CVout = list()  

for(CVnum in 1:5){
  
  #CVnum=1
  
  #mask
  pheno.cv = masked
    #mask phenotype columns
    N = dim(pheno.cv)[2]
    pheno.cv[pheno.cv$CV == CVnum, trait.key$tn] = -999

  #output for pheno
  pheno.out = pheno.cv %>% 
    select(-CV)
  
  write.table(pheno.out, "blup/pheno-cv.txt", 
              col.names = F, row.names = F, quote = F)
  
  #RUN BLUP script
  cmd = paste("scripts/cv-multi.sh", targetParam, CVnum)
  system(cmd)
  
  
  #pull solutions: multi-trait
  pull = paste0("data/sol-cv", CVnum, ".txt")
  sol = read.delim(pull, sep = "", header = F)[-1,]
  colnames(sol) = c("trait", "effect", "level" ,"solution", "se")
  sol = sol %>% mutate(se = as.numeric(se),
                       solution = as.numeric(solution))  
  
  #CV correlation
    #pull masked predictions

    sol.tmp = sol %>% filter(effect == 1) %>% 
      left_join(masked %>% 
                  select(level = 1, gc_id, CV, locid, loc) %>% 
                  mutate(level = as.character(level)),
                by = 'level') %>% 
      left_join(trait.key.used, by = "trait") %>%
      filter(CV == CVnum) %>% 
      select(gc_id, pheno.est = solution, tn)
      
#       %>% 
#       #add fixed effect of location
#       left_join( fixed.eff %>% 
#         select(locid = level, loceff = solution, trait) %>%
#         mutate(locid = as.numeric(locid)), by = c('locid', 'trait')) %>%
#       mutate(pheno.est = solution)
    
    #compare against real phenos
    
    #pull real phenos and adjust by all fixed effects
    
    realpheno = masked[masked$CV == CVnum, c("gc_id", args.split)] %>%
      pivot_longer(cols = all_of(args.split), names_to = "tn", values_to = "pheno.raw") %>%
      left_join(masked %>% select(gc_id, PC1, PC2, PC3, locid))
      
      realpheno$pheno.adj = NA
      
      #TODO: hardcoded, fix later
      realpheno = realpheno %>% left_join(trait.key.used)
      
      for(i in 1:nrow(realpheno)){
        #location
        X = fixed.eff$solution[fixed.eff$trait == realpheno$trait[i] &
                               fixed.eff$effect == 2 &
                               fixed.eff$level == realpheno$locid[i]]
        #PC1
        X = X + (fixed.eff$solution[fixed.eff$effect == 3 &
                                   fixed.eff$level == 1 &
                                   fixed.eff$trait == realpheno$trait[i] ] *
                 realpheno$PC1[i])
        #PC2
        X = X + (fixed.eff$solution[fixed.eff$effect == 4 &
                                   fixed.eff$level == 1 &
                                   fixed.eff$trait == realpheno$trait[i] ] *
                 realpheno$PC2[i])
        #PC3
        X = X + (fixed.eff$solution[fixed.eff$effect == 5 &
                                   fixed.eff$level == 1 &
                                   fixed.eff$trait == realpheno$trait[i] ] *
                 realpheno$PC3[i])
                 
      realpheno$pheno.adj[i] = realpheno$pheno.raw[i] - X

      }

    sol.tmp = sol.tmp %>% 
      left_join(realpheno %>% select(gc_id, tn, pheno.adj), by = c('gc_id', 'tn'))
    
    #correlation
    cv = sol.tmp %>% 
      group_by(tn) %>%
      summarise(cor = cor(pheno.est, pheno.adj),
                slope = cor * sd(pheno.adj) / sd(pheno.est))
    
    CVout[[CVnum]] = cv
  
}

#save(CVout, file = "CVresult.Rdat")
#("CVresult15JAN26.Rdat")

#reformat as df
CVdf = CVout[[1]]
CVdf$CV = 1
for(i in 2:5){
  x = CVout[[i]]
  x$CV = i
  CVdf = rbind(CVdf, x)
}

CVdf

write.csv(CVdf, file = "data/CV_summary.csv", row.names = F)

print(CVdf %>% group_by(tn) %>% 
  summarise(meancor = mean(cor),
            meanslope = mean(slope)))
