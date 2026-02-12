#take a command line argument


#target script for CV runs
targetScript = "qH-HMP-cv.par1"

fixed.key = data.frame(effect = c(2:4) %>% as.character(),
                       name = c("apnumid", "start", "requeen"))

#key for trait ids
trait.key = data.frame(tn = c("honey", "mites.adj", "pattern"),
                       trait = c(1,2,3) %>%  as.character())

#pull fixed effects from full run
fixed.eff = read.delim("blup/sol-quant-2JAN26.txt", sep = "", header = F)[-1,]
colnames(fixed.eff) = c("trait", "effect", "level" ,"solution", "se")
fixed.eff = fixed.eff %>% mutate(se = as.numeric(se),
                                 solution = as.numeric(solution)) %>% 
  filter(effect>1) %>% 
  left_join(fixed.key, by = 'effect')

#create list for CV
#try 10-fold
set.seed(2026)
#only mask individuals with full pheno data
masked = pheno.f
masked$CV = sample(c(1:5), nrow(pheno.f), replace = T)
table(masked$CV)

#create output object
CVout = list()  

for(CVnum in 1:5){
  
  #mask
  pheno.cv = masked
  pheno.cv[pheno.cv$CV == CVnum, trait.key$tn] = -999
  
  #check order for rows
  sum(colnames(H.f) != pheno.cv$colony_id)
  
  #output for pheno
  pheno.out = pheno.cv %>% 
    select(qnumid, apnumid, ynumid, synumid, start, requeen, pattern, honey,
           mites.adj, dwv.bin, bd.bin)
  
  write.table(pheno.out, "blup/cv-H-pheno.txt", 
              col.names = F, row.names = F, quote = F)
  
  #RUN BLUP script
  cmd = paste("blup/cv-H-blup.sh", targetScript, CVnum)
  system(cmd)
  
  
  #pull solutions: 3-trait
  pull = paste0("blup/cv-H-", CVnum, "-sol.txt")
  sol = read.delim(pull, sep = "", header = F)[-1,]
  colnames(sol) = c("trait", "effect", "level" ,"solution", "se")
  sol = sol %>% mutate(se = as.numeric(se),
                       solution = as.numeric(solution))  
  
  #CV correlation
  #tmp mask
  sol.tmp = sol %>% filter(effect == 1) %>% 
    left_join(masked %>% 
                select(level = qnumid, colony_id, CV, fixed.key$name) %>% 
                mutate(level = as.character(level)),
              by = 'level') %>% 
    filter(CV == CVnum) %>% 
    select(-effect, -level)
  
  #adjust solutions for fixed effects: Y = X + G
  sol.tmp$pheno.est = NA
  for(i in 1:nrow(sol.tmp)){
    #loop through fixed effects: maybe a matrix algebra method?
    X = 0
    for(j in 1:nrow(fixed.key)){
      X = X + fixed.eff$solution[
        fixed.eff$name == fixed.key$name[j] &
          fixed.eff$level == sol.tmp[i,fixed.key$name[j]] &
          fixed.eff$trait == sol.tmp$trait[i]
      ]
    }
    sol.tmp$pheno.est[i] = sol.tmp$solution[i] + X
    
  }
  
  #compare against real phenos
  realpheno = masked %>% 
    filter(CV == CVnum) %>% 
    select(colony_id, trait.key$tn) %>% 
    pivot_longer(cols = trait.key$tn, names_to = "tn", values_to = "pheno.real") %>% 
    left_join(trait.key, by = 'tn')
  
  sol.tmp = sol.tmp %>% 
    left_join(realpheno, by = c('trait', 'colony_id'))
  
  #correlation
  cv = sol.tmp %>% 
    select(colony_id, tn, pheno.est, pheno.real) %>% 
    left_join(pheno %>% select(colony_id, year), by = 'colony_id') %>% 
    group_by(tn) %>% 
    summarise(cor = cor(pheno.est, pheno.real),
              slope = cor * sd(pheno.real) / sd(pheno.est))
  
  CVout[[CVnum]] = cv
  
  #plot
  ggplot(sol.tmp, aes(x = pheno.est, y = pheno.real)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = 'lm', se = F) +
    facet_wrap(facets = vars(tn), nrow = 3) +
    labs(y = "Real Phenotype", x = "Estimated Phenotype")
  
  #TODO: report cv to outdf
  
}

save(CVout, file = "CVresult15JAN26.Rdat")
#("CVresult15JAN26.Rdat")

#reformat as df
CVdf = CVout[[1]]
CVdf$CV = 1
for(i in 2:5){
  x = CVout[[i]]
  x$CV = i
  CVdf = rbind(CVdf, x)
}

#plot cor
ggplot(CVdf, aes(x = tn, y = cor)) + 
  geom_boxplot() + 
  geom_point() + 
  lims(y = c(0,1))


# #accuracy
#   h = data.frame(tn = trait.key$tn,
#                  g = c( 0.26220, 0.10237, 0.12331),
#                  p = sapply(pheno.complete[,trait.key$tn], var))
#   h$h2 = h$g/h$p
#   h$h = sqrt(h$h2)
#   
#   CVdf$acc = apply(CVdf, 1, function(x){
#     as.numeric(x["cor"]) / h$h[h$tn == x["tn"]]
#   })
#   
#plot all
CVdf.long = CVdf %>%
  pivot_longer(cols = c("cor", "slope"), names_to = "parameter")

CVdf.cols = CVdf.long %>% group_by(tn, parameter) %>% 
  summarise(mean = mean(value), se = se(value))

CVdf.long %>% 
  #filter(parameter %in% c("acc", "cor")) %>% 
  ggplot(aes(x = tn, y = value, color = parameter)) +
  geom_col(data = CVdf.cols, 
           aes(y = mean, fill = parameter), 
           position = "dodge", alpha = 0.5) +
  geom_errorbar(data = CVdf.cols,
                aes(y = mean, ymin = mean-se, ymax = mean+se),
                position = position_dodge(width = 0.85),
                width = 0.5) +
  geom_point(position = position_dodge(width = 0.85)) + 
  labs(x = "trait")


CVdf.long %>% group_by(tn, parameter) %>% 
  summarise(mean = mean(value), se = se(value)) %>% 
  arrange(parameter)
