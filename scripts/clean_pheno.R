library(dplyr)
library(readxl)
library(lubridate)
library(Matrix)

select=dplyr::select


#read phenotypic data
pheno = read_excel("pheno/phenotypes.xlsx") %>% 
  mutate(pheno_id = gsub(" ", "", QC)) 
  
#data cleaning and prep

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
    
    #pull vernal day
    pheno$vernal = format(pheno$date, "%m-%d")
    pheno$vernal2 = as.Date(pheno$vernal, format = "%m-%d") 
    
    #assign season
    szn.cut = as.Date(c("2026-03-01", "2026-06-01", "2026-08-15", "2026-10-15"))
    pheno$season = NA
    pheno$season[pheno$vernal2 >= szn.cut[1] & pheno$vernal2 < szn.cut[2]] = "spring"
    pheno$season[pheno$vernal2 >= szn.cut[2] & pheno$vernal2 < szn.cut[3]] = "summer"
    pheno$season[pheno$vernal2 >= szn.cut[3] & pheno$vernal2 < szn.cut[4]] = "fall"
    pheno$season[pheno$vernal2 >= szn.cut[4] | pheno$vernal2 < szn.cut[1]] = "winter"
    
#     #show distribution
#     png(file = "data/season-hist.png")
#     
#       ggplot(pheno, aes(x = vernal2, fill = season)) + geom_histogram(bins=20)
#     
#     dev.off()
    

  #standardize locations
  loc.trans = data.frame(Location = c("Hawaii", "Georgia", "Southern California", "Minnesota",
                                      "Northern California", "Washington", "West Virginia",
                                      "Michigan", "Unknown", "NCA", "NC", "USA", "MN", 
                                      "GA", "HI", "PA", "CA", "OH", "NY", "WA", "SCA",
                                      "VA", "AL", "FL", "OR", "Oregon"), 
                        
                        loc.fix = c("HI", "GA", "CA", "MN",
                                    "CA", "WA", "WV",
                                    "MI", "USA", "CA", "NC", "USA", "MN", 
                                    "GA", "HI", "PA", "CA", "OH", "NY", "WA", "CA",
                                    "VA", "AL", "FL", "OR", "OR"))
  pheno = pheno %>% left_join(loc.trans, by = "Location")
  
  
    
#     regions = data.frame(loc.fix = c( "HI", "GA", "SCA", "MN", "NCA", "WA", "WV", "MI",
#                                       "USA", "NC", "PA", "CA", "OH", "NY", "VA", "AL", "FL", "OR"),
#                                       
#                          region = c( "HI", "SE", "CA", "MW", "CA", "NW", "EC", "MW",
#                                       "USA", "EC", "NE", "CA", "MW", "NE", "EC", "SE", "SE", "NW")
#                         )
#                         
#     pheno = pheno %>% left_join(regions)
    
    

    #remove duplicate entries
    pheno = pheno[-(which(pheno$pheno_id == "QC2573")[2]),]
    pheno = pheno[-(which(pheno$pheno_id == "QC2422")[2]),]
    
    #pull all names from sequencer
    
    allnames = read.delim("data/sample.fnames", header = F, sep = "") %>% 
      select(gc_id = V2) %>% 
      filter(grepl("QC", gc_id))
    
    #fix some incorrect sample ID's from the sequencer to match phenotype ID's
    bradley = read.csv("pheno/bradley_edits.csv")
    
    allnames = allnames %>% left_join(bradley %>% select(gc_id = manifest_id, new_id), by = 'gc_id')
      allnames$new_id[is.na(allnames$new_id)] = allnames$gc_id[is.na(allnames$new_id)] 
    
    # #how many fail to match?
    allnames$gc_id[!allnames$gc_id %in% pheno$pheno_id]
    allnames$new_id[!allnames$new_id %in% pheno$pheno_id]
      
      
    #manually drop duplicated genomic samples
      allnames = allnames[-which(allnames$new_id == "QC2573")[1],]
      allnames = allnames[-which(allnames$new_id == "QC3371")[1],]
  
      
    #amend pheno with genetic ids
      pheno.fix = pheno %>% 
        left_join(allnames %>% 
                    select(pheno_id = new_id, gc_id), by = 'pheno_id')
                    
    #manually drop duplicated pheno id
    pheno.fix = pheno.fix[-which(pheno.fix$pheno_id == "QC2422")[2],]
    
    #numeric class for numeric data
    pheno.num = pheno.fix %>% as.data.frame()
      
      for(COL in c("m.Body", "v.Sperm", "l.Sperm", "t.Sperm", "w.Head",
                   "w.Thorax", "d.Spermatheca", "f.Spermatheca")){
        pheno.num[,COL] = as.numeric(pheno.num[,COL])
      }
    
      
    
    #write all phenotypes
    write.csv(pheno.num, "data/all_pheno.csv",
              row.names = F, quote= F)
    
    #remove outlier phenotypes
    pheno.num = pheno.num[-(which.min(pheno.num$m.Body)),]
    
    #TODO
    #is this still needed?
    
#     #manually remove this problematic individual
#     pheno.fix = pheno.fix[-which(pheno.fix$gc_id == "QC0758"),]
#     
    #write cleaned phenotypes
    write.csv(pheno.num, "data/cleaned_pheno.csv",
              row.names = F, quote= F)
              
    #write intersection of pheno and genoids 
    write.table(file = "data/phenotyped.gcnames",
                pheno.num$gc_id[!is.na(pheno.num$gc_id)],
                col.names = F, row.names = F, quote = F)
                
      # #for each target loc
      #     for (LOC in c("HI","CA", "GA")){
      #       towrite = pheno.num$gc_id[!is.na(pheno.num$gc_id) & pheno.num$loc.fix == LOC]
      #       #format for plink
      #       towrite = data.frame(V1 = towrite, V2 = towrite)
      #     write.table(file = paste0("data/phenotyped_", LOC, ".gcnames"),
      #           towrite,
      #           col.names = F, row.names = F, quote = F)
      #           }
                
  
  
