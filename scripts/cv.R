#take a command line argument
args = commandArgs(trailingOnly=TRUE)
#args = "w"

#target script for CV runs
targetParam = paste0(args[1], "-cv.par1")

print(targetParam)
print(getwd())

library(dplyr)

#wd is queen quality!


# fixed.key = data.frame(effect = c(2:4) %>% as.character(),
#                        name = c("apnumid", "start", "requeen"))
# 
#key for trait ids
trait.key = data.frame(tn = c("l", "w", "v", "t"),
                       trait = c(1:4) %>%  as.character())

#pull fixed effects from full run
fixed.eff = read.delim(paste0("data/sol-", args[1], ".txt"), sep = "", header = F)[-1,]
colnames(fixed.eff) = c("trait", "effect", "level" ,"solution", "se")
fixed.eff = fixed.eff %>% mutate(se = as.numeric(se),
                                 solution = as.numeric(solution)) %>% 
  filter(effect == 2)
  
#load sample info from pheno file
pheno = read.delim("blup/pheno.txt", header = F, sep = " ")
  #grab colnames of sampleids and locations
  samp.n = which(grepl("QC", pheno))
  loc.n = which(grepl("HI", pheno))
  

#create list for CV
  #only mask individuals with full pheno data
  masked = pheno
    colnames(masked)[samp.n] = "id"
    colnames(masked)[loc.n] = "loc"
    colnames(masked)[3:(2 + nrow(trait.key))] = trait.key$tn
  
  
  set.seed(2026)
  masked$CV = sample(c(1:5), nrow(pheno), replace = T)
  table(masked$CV)

#create output object
CVout = list()  

for(CVnum in 1:5){
  
  #CVnum=1
  
  #mask
  pheno.cv = masked
    #mask phenotype columns
    N = dim(pheno.cv)[2]
    pheno.cv[pheno.cv$CV == CVnum, c(3:(N-3))] = -999

  #output for pheno
  pheno.out = pheno.cv %>% 
    select(-CV)
  
  write.table(pheno.out, "blup/pheno-cv.txt", 
              col.names = F, row.names = F, quote = F)
  
  #RUN BLUP script
  cmd = paste("scripts/cv.sh", targetParam, CVnum)
  system(cmd)
  
  
  #pull solutions: 1-trait
  pull = paste0("data/sol-cv", CVnum, ".txt")
  sol = read.delim(pull, sep = "", header = F)[-1,]
  colnames(sol) = c("trait", "effect", "level" ,"solution", "se")
  sol = sol %>% mutate(se = as.numeric(se),
                       solution = as.numeric(solution))  
  
  #CV correlation
    #pull masked predictions
    sol.tmp = sol %>% filter(effect == 1) %>% 
      left_join(masked %>% 
                  select(level = 1, id, CV, locnum = 2, loc) %>% 
                  mutate(level = as.character(level)),
                by = 'level') %>% 
      filter(CV == CVnum) %>% 
      select(-effect, -level) %>% 
      #add fixed effect of location
      left_join( fixed.eff %>% 
        select(locnum = level, loceff = solution) %>%
        mutate(locnum = as.numeric(locnum)), by = 'locnum') %>%
      mutate(pheno.est = solution + loceff)
    
    #compare against real phenos
    realpheno = masked[masked$CV == CVnum, c("id", args[1])]
      colnames(realpheno)[2] = "pheno.real"

    
    sol.tmp = sol.tmp %>% 
      left_join(realpheno, by = 'id')
    
    #correlation
    cv = sol.tmp %>% 
      select(id, pheno.est, pheno.real) %>% 
      summarise(cor = cor(pheno.est, pheno.real),
                slope = cor * sd(pheno.real) / sd(pheno.est))
    
    CVout[[CVnum]] = cv
  
#   #plot
#   pdf(file = "plot.pdf")
#   plot(sol.tmp$pheno.est, sol.tmp$pheno.real)
#   dev.off()
  #TODO: report cv to outdf
  
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

print(CVout)

# #plot cor
# ggplot(CVdf, aes(x = tn, y = cor)) + 
#   geom_boxplot() + 
#   geom_point() + 
#   lims(y = c(0,1))
# 
# 
# # #accuracy
# #   h = data.frame(tn = trait.key$tn,
# #                  g = c( 0.26220, 0.10237, 0.12331),
# #                  p = sapply(pheno.complete[,trait.key$tn], var))
# #   h$h2 = h$g/h$p
# #   h$h = sqrt(h$h2)
# #   
# #   CVdf$acc = apply(CVdf, 1, function(x){
# #     as.numeric(x["cor"]) / h$h[h$tn == x["tn"]]
# #   })
# #   
# #plot all
# CVdf.long = CVdf %>%
#   pivot_longer(cols = c("cor", "slope"), names_to = "parameter")
# 
# CVdf.cols = CVdf.long %>% group_by(tn, parameter) %>% 
#   summarise(mean = mean(value), se = se(value))
# 
# CVdf.long %>% 
#   #filter(parameter %in% c("acc", "cor")) %>% 
#   ggplot(aes(x = tn, y = value, color = parameter)) +
#   geom_col(data = CVdf.cols, 
#            aes(y = mean, fill = parameter), 
#            position = "dodge", alpha = 0.5) +
#   geom_errorbar(data = CVdf.cols,
#                 aes(y = mean, ymin = mean-se, ymax = mean+se),
#                 position = position_dodge(width = 0.85),
#                 width = 0.5) +
#   geom_point(position = position_dodge(width = 0.85)) + 
#   labs(x = "trait")
# 
# 
# CVdf.long %>% group_by(tn, parameter) %>% 
#   summarise(mean = mean(value), se = se(value)) %>% 
#   arrange(parameter)
