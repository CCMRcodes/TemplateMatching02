################################################################################
#
# Template Matching: Selecting templates
#
# Below code is for runs with a seperate template for each
# hospital tier and a fixed proportion of operative hospitalizations (run 8)
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
library(here)


##### Loading Data and defining variable groups #####
setwd(here("ORD_Prescott_Comparing","5. Identifiable Data (sensitive data)"))
year = "2017"
df = as.data.frame(readRDS(paste("IPEC data/","ipec",year,"_working.rds",sep="")))

#defining variable groups
labs_raw = c("wbc_unit", "albumin_unit", "bilirubin_unit", "bun_unit",
             "glucose_unit", "htc_unit", "sodium_unit", "p02_unit", 
             "pco2_unit", "ph_unit", "gfr_unit")
#full comorbidity list:
# comorbidities = c("chf", "dm_comp", "dm_uncomp", "hypothyroid", "renal", "liver", "pud", 
#                   "immunedef", "lymphoma", "cancer_met", "cancer_nonmet", "ra", "coag", 
#                   "obesity", "wtloss", "fen", "anemia_cbl", "anemia_def", "etoh", 
#                   "drug", "psychoses", "valvular_d2", "depression", "pulm_circ", "pvd", 
#                   "htn", "paralysis", "neuro", "pulm")
#reduced comorbidity list:
comorbidities = c("chf", "renal", "liver",  "cancer_met", "depression", "paralysis", "pulm")

demographics = c("female", "black", "white", "hispanic", "age")
#note: leaving out dxgroup_other indicator to avoid computational issues because of singular matrix
dx_groups = c("dxgroup_cardiovascular", "dxgroup_gastro", "dxgroup_infection", 
              "dxgroup_psych", "dxgroup_renal", "dxgroup_resp") 

##### Choosing template #####
#note: We chose a quantity of operative patients such that the percentage in the template is 
#the same as the percentage in the median hospital within each tier/ICU level, rounding down
#Level 1: 9.5% operative (90.5%) = 272 op
#Level 2: 8.7% operative (91.3%) = 274
#Level 3: 7.9% operative (92.1%) = 277
#Level 4: 0.6% operative (99.45) = 0 (adjusted down to 0 to avoid irrelevant exclusions)
#No ICU:  0.0% operative (100%)  = 0

df_modular = df %>%
  group_by(facility_name) %>%
  mutate(operative_count =  sum(operative))

full_template_size = 300

distVarlist = c("facility_id","id_patient", "id_hospitalization", "pred_mort", demographics,comorbidities, dx_groups, labs_raw,"mort30")

numTemplateSamples = 1000
icu_level_list = c("Level 1","Level 2","Level 3","Level 4","No ICU")

##### Selecting non-operative templates for each tier #####
for (level in icu_level_list){
  df_maha = df %>% 
    filter(icu_level == level) %>% 
    filter(operative==0) %>%
    select(distVarlist)
  numVars = dim(df_maha)[2]-3 
  
  if (level == "Level 1"){
    template_size = ceiling(full_template_size*.905)
  } else if (level == "Level 2") {
    template_size = ceiling(full_template_size*.913)
  } else if (level == "Level 3") {
    template_size = ceiling(full_template_size*.921)
  } else if (level == "Level 4") {
    template_size = 300
  } else if (level == "No ICU") {
    template_size = 300
  }
  
  #generate 500 random samples and calculate mahalanobis distance (based on means) from overall dataset
  set.seed(10)
  tic("Selecting template")
  print(paste("Choosing a Template for ICU ", level, sep = ""))
  
  templateDistances = matrix(c(1:numTemplateSamples, rep(0,numTemplateSamples) ), 
                             nrow = numTemplateSamples, ncol = 2, byrow = F, 
                             dimnames = list(c(1:numTemplateSamples),c("Num", "mDist") ) )
  
  for (i in 1:numTemplateSamples) {
    tmp = df_maha[sample(nrow(df_maha), template_size, replace = F), ]
    #D.sq calculated mahalanobis distance
    templateDistances[i,2] = D.sq(df_maha[,4:dim(df_maha)[2]],tmp[,4:dim(df_maha)[2]])$D.sq
    
    assign(paste("template_",i, sep = "" ),tmp)
    if (i %% 50 == 0){
      print(paste(i," samples generated..."))
    }
  }
  rm(tmp)
  
  bestTemplate = templateDistances[which.min(templateDistances[,2]),]
  template = get(paste("template_",as.character(bestTemplate[1]), sep = ""))
  
  #cleaning environment
  for (i in 1:numTemplateSamples) { rm(list = paste("template_",i,sep = "")) }
  toc()
  
  saveRDS(template$id_hospitalization, file = paste("IPEC data/template_icu_mod_non_op_",level,".rds",sep = ""))
  
}
rm(bestTemplate, templateDistances, df_maha, i, numTemplateSamples, numVars)


##### Selecting operative templates for each tier #####
for (level in icu_level_list){
  df_maha = df %>% 
    filter(icu_level == level) %>% 
    filter(operative==1) %>%
    select(distVarlist)
  numVars = dim(df_maha)[2]-3 
  
  if (level == "Level 1"){
    template_size = 300 - ceiling(full_template_size*.905)
  } else if (level == "Level 2") {
    template_size = 300 - ceiling(full_template_size*.913)
  } else if (level == "Level 3") {
    template_size = 300 - ceiling(full_template_size*.921)
  } else if (level == "Level 4") {
    next #template of 0
  } else if (level == "No ICU") {
    break #template of 0
  }
  print(template_size)
  
  set.seed(11)
  tic("Selecting template")
  print(paste("Choosing a Template for ICU ", level, sep = ""))
  
  templateDistances = matrix(c(1:numTemplateSamples, rep(0,numTemplateSamples) ), 
                             nrow = numTemplateSamples, ncol = 2, byrow = F, 
                             dimnames = list(c(1:numTemplateSamples),c("Num", "mDist") ) )
  
  for (i in 1:numTemplateSamples) {
    tmp = df_maha[sample(nrow(df_maha), template_size, replace = F), ]
    #calculate mahalanobis distance
    templateDistances[i,2] = D.sq(df_maha[,4:dim(df_maha)[2]],tmp[,4:dim(df_maha)[2]])$D.sq
    
    assign(paste("template_",i, sep = "" ),tmp)
    if (i %% 50 == 0){
      print(paste(i," samples generated..."))
    }
  }
  rm(tmp)
  
  bestTemplate = templateDistances[which.min(templateDistances[,2]),]
  template = get(paste("template_",as.character(bestTemplate[1]), sep = ""))
  
  #cleaning environment
  for (i in 1:numTemplateSamples) { rm(list = paste("template_",i,sep = "")) }
  toc() #last run: 
  
  saveRDS(template$id_hospitalization, file = paste("IPEC data/template_icu_mod_op_",level,".rds",sep = ""))
  
}

