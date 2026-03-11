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
    
    #manually remove this problematic individual
    pheno.fix = pheno.fix[-which(pheno.fix$gc_id == "QC0758"),]
    
#TODO: get a count for how many phenotypes actually exist...
      
      
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
              
    #write intersection of pheno and genoids 
    write.table(file = "data/phenotyped.gcnames",
                pheno.num$gc_id[!is.na(pheno.num$gc_id)],
                col.names = F, row.names = F, quote = F)
                
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
