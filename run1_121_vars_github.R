################################################################################
#
# Template Matching run 1
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
grid_path  = "/cifs/vhacdwfpcfs02/vhaannmollid/Projects/"
setwd(paste(grid_path,"ORD_Prescott_Comparing/5. Identifiable Data (sensitive data)" ,sep = ""))
df = as.data.frame(readRDS(paste("IPEC data/","ipec",year,"_working.rds",sep="")))

#defining variable groups
lab_scores = c("wbc_unit_score", "albumin_unit_score", "bilirubin_unit_score", "bun_unit_score",
               "glucose_unit_score", "hct_unit_score", "sodium_unit_score", "pao2_unit_score", 
               "pco2_unit_score", "ph_unit_score", "gfr_unit_score")
#comorbidities = c("chf", "renal", "liver", "cancer_met", "depression", "paralysis", "pulm")
comorbidities = c("chf", "dm_comp", "dm_uncomp", "hypothyroid", "renal", "liver", "pud", 
                  "immunedef", "lymphoma", "cancer_met", "cancer_nonmet", "ra", "coag", "obesity", "wtloss", "fen", 
                  "anemia_cbl", "anemia_def", "etoh", "drug", "psychoses", "valvular_d2", "depression", "pulm_circ",
                  "pvd", "htn", "paralysis", "neuro", "pulm")

demographics = c("female", "black", "white", "hispanic", "age")

template_size = 300

###### OPTIMAL MATCHING with RCBalance #####
num_facilities = length(unique(df$facility_id))
template_id = readRDS(paste("IPEC data/template_id_",year,".rds",sep = ""))
template = filter(df, id_hospitalization %in% template_id)
template$case = 1
match_var_list = c("pred_mort", lab_scores, comorbidities, "admission_source", demographics, "unit_dx_ccs")
#initialize matrix to store IDs for matched patient hospitalizations
matched_ids = matrix(c(rep(0,template_size*num_facilities), rep(1:template_size,num_facilities) ), 
                     nrow = num_facilities*template_size, ncol = 2, byrow = F, 
                     dimnames = list(c(1:(num_facilities*template_size)),c("id_hospitalization", "id_match") ) )

tic("RCBalance matching")
print("Beginning RCBalance Matching")

counter = 0
for (site in unique(df$facility_name) ) {
  counter = counter+1
  if (counter %% 10 == 0){
    print(paste(counter, " sites completed..."))
  }
  print(paste("Matching: ",site) )
  #keep a single facility to match to
  facility = filter(df, facility_name == site)
  facility$case = 0
  test_data = rbind(facility,template)
  
  #make distance structure
  my.dist <- build.dist.struct(z=test_data$case, X=test_data[match_var_list], calip.option = "none")
  
  #compute match
  match.out <- rcbalance(my.dist,  
                         treated.info = template, 
                         control.info = facility,						
                         exclude.treated = FALSE)
  matrix_start = (counter-1)*template_size + 1
  matched_ids[ matrix_start:(matrix_start+299) ,1] = test_data[match.out$matches,"id_hospitalization"]
  #if(counter ==2) {break}
} 
toc()

saveRDS(matched_ids, file = "Matching results/matched_ids_no_fb_70_vars.rds")