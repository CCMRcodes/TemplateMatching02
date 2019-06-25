################################################################################
#
# Template Matching run 7
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
library(asbio) 

year = "2017"
grid_path  = "/cifs/vhacdwfpcfs02/vhaannmollid/Projects/"
setwd(paste(grid_path,
            "ORD_Prescott_Comparing/5. Identifiable Data (sensitive data)", sep = ""))
df = as.data.frame(readRDS(paste("IPEC data/","ipec",year,"_working.rds",sep="")))

#defining variable groups
lab_scores = c("wbc_unit_score", "albumin_unit_score", "bilirubin_unit_score", "bun_unit_score",
               "glucose_unit_score", "hct_unit_score", "sodium_unit_score", "pao2_unit_score", 
               "pco2_unit_score", "ph_unit_score", "gfr_unit_score")
comorbidities = c("chf", "renal", "liver", "cancer_met", "depression", "paralysis", "pulm")
demographics = c("female", "black", "white", "hispanic", "age")

template_size = 300

###### OPTIMAL MATCHING with RCBalance #####

match_var_list = c("pred_mort", lab_scores, comorbidities, "admission_source", demographics)
icu_level_list = c("Level 1","Level 2","Level 3","Level 4","No ICU")


l1 = c("pred_mort_cat5")

for (level in icu_level_list){
  df_icu = df %>% filter(icu_level == level)
  num_facilities = length(unique(df_icu$facility_id))
  template_id = readRDS(paste("IPEC data/template_icu_",level,".rds",sep = ""))
  template = filter(df_icu, id_hospitalization %in% template_id)
  template$case = 1
  
  #initialize matrix to store IDs for matched patient hospitalizations
  matched_ids = matrix(c(rep(0,template_size*num_facilities), rep(1:template_size,num_facilities) ), 
                       nrow = num_facilities*template_size, ncol = 2, byrow = F, 
                       dimnames = list(c(1:(num_facilities*template_size)),c("id_hospitalization", "id_match") ) )
  
  tic("RCBalance matching")
  print(paste("Beginning RCBalance Matching for ICU ",level, sep = ""))
  
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
    
    #make distance structure
    my.dist <- build.dist.struct(z=test_data$case, X=test_data[match_var_list], calip.option = "none")
    
    #compute match
    match.out <- rcbalance(my.dist, 
                           fb.list = list(l1),
                           near.exact = c("unit_dx_grouped","operative"), 
                           treated.info = template, 
                           control.info = facility, 
                           exclude.treated = FALSE,
                           tol = 1e-3)
    
    matrix_start = (counter-1)*template_size + 1
    matched_ids[ matrix_start:(matrix_start+299) ,1] = test_data[match.out$matches,"id_hospitalization"]
    saveRDS(matched_ids, 
            file = paste("Matching results/matched_ids_fb_nearexact_dx_op","_icu_",level,".rds", sep = ""))
  } 
  toc()
}