################################################################################
#
# Template Matching run 10
#
# Date: 2019-06-21
# Author: Daniel Molling <daniel.molling@va.gov>
#
################################################################################

library(survival)
library(tictoc)
library(haven)
library(tidyr)
library(ggplot2)
library(optmatch)
library(plyr) 
library(rcbalance)
library(MatchIt)
library(dplyr)
library(knitr)
library(kableExtra)
library(asbio) #for mahalanobis distance

year = "2017"
grid_path  = "/cifs2/vhacdwfpcfs02/vhaannmollid/Projects/"
setwd(paste(grid_path,
            "ORD_Prescott_Comparing/5. Identifiable Data (sensitive data)", sep = ""))
df = as.data.frame(readRDS(paste("IPEC data/","ipec",year,"_working.rds",sep="")))

###### OPTIMAL MATCHING with RCBalance #####

comorbidities = c("chf", "renal", "liver", "cancer_met", "depression", "paralysis", "pulm")
demographics = c("female", "age") #excluding ethnicity 
admission_source = c("admit_nursing_home","admit_emergency_dept")

match_var_list = c("pred_mort", comorbidities, admission_source, demographics)
#ICU Level corresponds to "Tiers" in text (level 1 = tier 1, ect)
icu_level_list = c("Level 1","Level 2","Level 3","Level 4","No ICU")

#define variable that stores the count of operative patients at that hospital
df = df %>%
  group_by(facility_name) %>%
  mutate(operative_count =  sum(operative))

#Set fine balance groups
l1 = c("pred_mort_cat5")
for (operative_status in c("non_op","op")) {
  
  for (level in icu_level_list){
    
    if ( (operative_status=="op" & level == "Level 4") | (operative_status=="op" & level=="No ICU") ) {
      next
    }
    
    df_icu = df %>% filter(icu_level == level)
    template_id = readRDS(paste("IPEC data/template_icu_mod_",operative_status,"_",level,".rds",sep = ""))
    
    template = filter(df_icu, id_hospitalization %in% template_id)
    template$case = 1
    template_size = dim(template)[1]
    
    if (operative_status == "non_op" ){
      df_icu = df_icu %>% 
        filter(operative == 0)
    } else {
      ### Removing non-operative hospitalizations and hospitals that do not
      ### have enough operative hospitalizations to ensure a 3:1 match ratio
      df_icu = df_icu %>% 
        filter(operative == 1) %>%
        filter(operative_count >= template_size*3)
    }
    
    num_facilities = length(unique(df_icu$facility_id))
    #initialize matrix to store IDs for matched patient hospitalizations
    matched_ids = matrix(c(rep(0,template_size*num_facilities), rep(1:template_size,num_facilities) ), 
                         nrow = num_facilities*template_size, ncol = 2, byrow = F, 
                         dimnames = list(c(1:(num_facilities*template_size)),
                                         c("id_hospitalization", "id_match")
                         ) )
    
    tic("RCBalance matching")
    print(paste("Beginning RCBalance Matching for ICU ",level," ",operative_status, sep = ""))
    
    counter = 0
    for (site in unique(df_icu$facility_name) ) {
      counter = counter+1
      if (counter %% 10 == 0){
        print(paste(counter, " sites completed..."))
      }
      print(paste("Matching: ",site) )
      #keep a single facility to match to
      facility = filter(df_icu, facility_name == site)
      facility$case = 0
      test_data = rbind(facility,template)
      
      #generate mahalanobis distance structure
      my.dist <- build.dist.struct(z=test_data$case, X=test_data[match_var_list], calip.option = "none")
      
      #compute match
      match.out <- rcbalance(my.dist,
                             fb.list = list(l1),
                             near.exact = c("unit_dx_grouped"),
                             treated.info = template,
                             control.info = facility,
                             exclude.treated = FALSE,
                             tol = 1e-3)
      
      matrix_start = (counter-1)*template_size + 1
      result = test_data[match.out$matches,"id_hospitalization"]
      matched_ids[ matrix_start:(matrix_start+template_size-1) ,1] = unlist(result, use.names = F)
      #break
      #if(counter ==2) {break}
    } 
    saveRDS(matched_ids,  
            file = paste("Matching results/matched_ids_icu_mod_",
                         operative_status,"_",level,".rds", sep = ""))
    
    toc()
  }
  
}