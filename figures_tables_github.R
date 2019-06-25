################################################################################
#
# Template Matching Figures and Tables
#
# Date: 2019-06-21
# Author: Daniel Molling <daniel.molling@va.gov>
#
################################################################################

###Note: Code displayed below is for runs with seperate templates for each hospital
###      tier and fixed proportion of operative hospitalizations (runs 8,10)

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
library(here)
library(lattice)  
library(Formula)  
library(ggplot2)  
library(Hmisc)
library(naniar)
library(visdat)
library(lme4) 
library(crossmatch) 
library(data.table)
library(formattable)
library(gridExtra)
library(grid)
library(gtable)

year = "2017"
df = as.data.frame(readRDS(paste("IPEC data/","ipec",year,"_working.rds",sep="")))
df$hospital_tier = recode(df$icu_level,
                          "Level 1" = "Tier 1",
                          "Level 2" = "Tier 2",
                          "Level 3" = "Tier 3",
                          "Level 4" = "Tier 4",
                          "No ICU"  = "Tier 5"
)

##### Cross-match test for tiered runs with fixed operative proportions #####

run = "matched_ids_icu_mod_labs" #run 8

comorbidities = c("chf", "renal", "liver",  "cancer_met", "depression", "paralysis", "pulm")
admission_source = c("admit_nursing_home","admit_emergency_dept")
# demographics = c("female", "black", "white", "hispanic", "age")
demographics = c("female", "age")
labs_raw = c("wbc_unit", "albumin_unit", "bilirubin_unit", "bun_unit",
             "glucose_unit", "htc_unit", "sodium_unit", "p02_unit", 
             "pco2_unit", "ph_unit", "gfr_unit")
###alter list for different runs
match_var_list = c("pred_mort", comorbidities, admission_source, demographics, labs_raw) 
#note: excluding operative for fixed proportion runs

facility_list = unique(left_join(x = data.frame(facility_name = unique(df$facility_name)), 
                                 y = df[,c("facility_name","icu_level")]))

#create matrix to store results for all tiers for a single run
crossmatch_pvals = matrix(c(facility_list$facility_name,
                            facility_list$icu_level,
                            rep(NA,dim(facility_list)[1])), 
                          nrow = dim(facility_list)[1], ncol = 3) 

#ICU levels correspond to hospital tier (Level 1 = Tier 1, ect)
levels = c("Level 1","Level 2","Level 3","Level 4","No ICU")
for( level in levels) {
  #Loading matched data, keeping only facilities with sufficient operative patients
  if (level %in% c("Level 4","No ICU") ) {
    match_ids = as.data.frame(readRDS(paste("Matching results/",run,"non_op_",level,".rds", sep="")))
    df_match = filter(df, id_hospitalization %in% match_ids$id_hospitalization)
    
    template_id_non_op = readRDS(paste("IPEC data/template_icu_mod_non_op_",level,".rds",sep = ""))
    template = filter(df, id_hospitalization %in% template_id_non_op )
    
  } else {
    match_ids_op = as.data.frame(readRDS(paste("Matching results/",run,"op_",level,".rds", sep="")))
    match_ids_non_op = as.data.frame(readRDS(paste("Matching results/",run,"non_op_",level,".rds", sep="")))
    match_ids_non_op$id_match = as.character(as.numeric(match_ids_non_op$id_match) + 
                                               max(as.numeric(match_ids_op$id_match)))
    match_ids = bind_rows(match_ids_non_op, match_ids_op)
    
    df_match_op = filter(df, id_hospitalization %in% match_ids_op$id_hospitalization)
    
    df_match_non_op = df %>% 
      filter(id_hospitalization %in% match_ids_non_op$id_hospitalization) %>%
      filter(facility_name %in% unique(df_match_op$facility_name)) 
    df_match = bind_rows(df_match_op, df_match_non_op)
    
    #loading template
    template_id_op = readRDS(paste("IPEC data/template_icu_mod_op_",level,".rds",sep = ""))
    template_id_non_op = readRDS(paste("IPEC data/template_icu_mod_non_op_",level,".rds",sep = ""))
    template = filter(df, id_hospitalization %in% c(template_id_op,template_id_non_op) )
    
  }
  df_match$case = 0
  template$case = 1
  
  for (hospital in facility_list[which(facility_list$icu_level==level),"facility_name"]){
    hosp = filter(df_match, facility_name == hospital)
    X = rbind(template, hosp)
    z = X$case
    X = select(X, match_var_list)
    ## Rank based Mahalanobis distance between each pair
    X = as.matrix(X)
    n = dim(X)[1]
    k = dim(X)[2]
    for (j in 1:k) X[,j] = rank(X[,j])
    cv = cov(X)
    vuntied = var(1:n)
    rat = sqrt(vuntied/diag(cv))
    cv = diag(rat)%*%cv%*%diag(rat)
    out = matrix(NA,n,n)
    icov = ginv(cv)
    for (i in 1:n) out[i,] = mahalanobis(X, X[i,],icov,inverted=T)
    dis = out
    ## Cross-match test
    crossmatch_pvals[which(crossmatch_pvals==hospital),3] = crossmatchtest(z,dis)$approxpval
  }
}
sum(crossmatch_pvals[,3] <.05,na.rm = T)
table(crossmatch_pvals[,3]<.05, useNA = "always")

##### Post-match tables for tiered runs with fixed operative proportion counting mismatched vars  #####

run = "matched_ids_icu_mod_labs" #run 8

#ICU levels correspond to tier (Level 1 = Tier 1, ect)
levels = c("Level 1","Level 2","Level 3","Level 4","No ICU")

for( level in levels) {
  set.seed(10)
  #Loading matched data, keeping only facilities with sufficient operative patients
  if (level %in% c("Level 4","No ICU") ) {
    match_ids = as.data.frame(readRDS(paste("Matching results/",run,"non_op_",level,".rds", sep="")))
    df_match = filter(df, id_hospitalization %in% match_ids$id_hospitalization)
    
    cat_varlist = c("female", "black", "white", "hispanic", "depression", "liver", "cancer_met", "paralysis", "pulm", "renal",  "admit_nursing_home", "admit_emergency_dept", "dxgroup_cardiovascular", "dxgroup_psych", "dxgroup_infection", "dxgroup_gastro", "dxgroup_resp", "dxgroup_renal", "dxgroup_other") #"chf",
    
  } else {
    match_ids_op = as.data.frame(readRDS(paste("Matching results/",run,"op_",level,".rds", sep="")))
    match_ids_non_op = as.data.frame(readRDS(paste("Matching results/",run,"non_op_",level,".rds", sep="")))
    match_ids_non_op$id_match = as.character(as.numeric(match_ids_non_op$id_match) + 
                                               max(as.numeric(match_ids_op$id_match)))
    match_ids = bind_rows(match_ids_non_op, match_ids_op)
    
    df_match_op = filter(df, id_hospitalization %in% match_ids_op$id_hospitalization)
    
    df_match_non_op = df %>% 
      filter(id_hospitalization %in% match_ids_non_op$id_hospitalization) %>%
      filter(facility_name %in% unique(df_match_op$facility_name)) 
    df_match = bind_rows(df_match_op, df_match_non_op)
    cat_varlist = c("operative", "female", "black", "white", "hispanic", "depression", "liver", 
                    "cancer_met", "paralysis", "pulm", "renal",  "admit_nursing_home", "admit_emergency_dept", 
                    "dxgroup_cardiovascular", "dxgroup_psych", "dxgroup_infection", "dxgroup_gastro", "dxgroup_resp", 
                    "dxgroup_renal", "dxgroup_other") 
    
  }
  
  labs_raw = c("wbc_unit", "albumin_unit", "bilirubin_unit", "bun_unit",
               "glucose_unit", "htc_unit", "sodium_unit", "p02_unit", 
               "pco2_unit", "ph_unit", "gfr_unit")
  continuous_varlist = c("pred_mort","age", labs_raw)
  
  non_binary_varlist = c("pred_mort_cat5")
  
  tab1 = matrix(ncol = 9,nrow = 0)
  for (z in continuous_varlist)  {
    tmp = df %>%
      summarise(minimum = quantile(aggregate(df_match[,z],
                                             by = list(df_match$facility_id),median)$x, prob = c(0)),
                lower_eigth = quantile(aggregate(df_match[,z],
                                                 by = list(df_match$facility_id),median)$x, prob = c(.125)),
                lower_quartile = quantile(aggregate(df_match[,z],
                                                    by = list(df_match$facility_id),median)$x, prob = c(.25)),
                median = quantile(aggregate(df_match[,z],
                                            by = list(df_match$facility_id),median)$x, prob = c(.5)),
                upper_quartile = quantile(aggregate(df_match[,z],
                                                    by = list(df_match$facility_id),median)$x, prob = c(.75)),
                upper_eigth = quantile(aggregate(df_match[,z],
                                                 by = list(df_match$facility_id),median)$x, prob = c(.875)),
                maximum = quantile(aggregate(df_match[,z],
                                             by = list(df_match$facility_id),median)$x, prob = c(1)),
                chisq =   kruskal.test(as.formula(paste(z," ~ as.factor(facility_id)")), 
                                       data = df_match)$statistic,
                p_value = kruskal.test(as.formula(paste(z," ~ as.factor(facility_id)")) , 
                                       data = df_match)$p.value
      ) %>%
      as.matrix()
    tab1 = rbind(tab1,tmp)
  }               
  
  for (z in cat_varlist)  {
    tmp = df %>%
      summarise(minimum = quantile(aggregate(df_match[,z],
                                             by = list(df_match$facility_id),mean)$x, prob = c(0)),
                lower_eigth = quantile(aggregate(df_match[,z],
                                                 by = list(df_match$facility_id),mean)$x, prob = c(.125)),
                lower_quartile = quantile(aggregate(df_match[,z],
                                                    by = list(df_match$facility_id),mean)$x, prob = c(.25)),
                median = quantile(aggregate(df_match[,z],
                                            by = list(df_match$facility_id),mean)$x, prob = c(.5)),
                upper_quartile = quantile(aggregate(df_match[,z],
                                                    by = list(df_match$facility_id),mean)$x, prob = c(.75)),
                upper_eigth = quantile(aggregate(df_match[,z],
                                                 by = list(df_match$facility_id),mean)$x, prob = c(.875)),
                maximum = quantile(aggregate(df_match[,z],
                                             by = list(df_match$facility_id),mean)$x, prob = c(1)),
                chisq =   chisq.test(df_match[,z],df_match[,"facility_id"], simulate.p.value = T)$statistic,
                p_value = chisq.test(df_match[,z],df_match[,"facility_id"], simulate.p.value = T)$p.value
      ) %>%
      as.matrix()
    tab1 = rbind(tab1,tmp)
  }     
  
  for (z in non_binary_varlist)  {
    tmp = df %>%
      summarise(minimum = "N/A",
                lower_eigth = "N/A",
                lower_quartile = "N/A",
                median = "N/A",
                upper_quartile = "N/A",
                upper_eigth = "N/A",
                maximum = "N/A",
                chisq =   chisq.test(df_match[,z],df_match[,"facility_id"], simulate.p.value = T)$statistic,
                p_value = chisq.test(df_match[,z],df_match[,"facility_id"], simulate.p.value = T)$p.value
      ) %>%
      as.matrix()
    tab1 = rbind(tab1,tmp)
  }     
  
  rownames(tab1) = c(continuous_varlist, cat_varlist, non_binary_varlist)
  write.table(tab1, file = paste("Matching results/Charts and tables/",run,level,"_table1.csv", sep = ""), sep=",")
  
}





##### Post-match regressions and figures for tiered runs with fixed operative proportions #####
#note: model is the same for non-tiered or unfixed op proportion runs - just remove looping structure
run = "matched_ids_icu_mod_labs" #run 8

levels = c("Level 1","Level 2","Level 3","Level 4","No ICU")

for( level in levels) {
  tier = case_when(
    level == "Level 1" ~ "Tier 1",
    level == "Level 2" ~ "Tier 2",
    level == "Level 3" ~ "Tier 3",
    level == "Level 4" ~ "Tier 4",
    level ==  "No ICU" ~ "Tier 5"
  )
  #Loading matched data, keeping only facilities with sufficient operative patients
  if (level %in% c("Level 4","No ICU") ) {
    match_ids = as.data.frame(readRDS(paste("Matching results/",run,"non_op_",level,".rds", sep="")))
    df_match = filter(df, id_hospitalization %in% match_ids$id_hospitalization)
    #merge in match id to be used as random effect
    df_match = merge(df_match, match_ids, by = "id_hospitalization", all.x = T)
  } else {
    match_ids_op = as.data.frame(readRDS(paste("Matching results/",run,"op_",level,".rds", sep="")))
    match_ids_non_op = as.data.frame(readRDS(paste("Matching results/",run,"non_op_",level,".rds", sep="")))
    match_ids_non_op$id_match = as.character(as.numeric(match_ids_non_op$id_match) + 
                                               max(as.numeric(match_ids_op$id_match)))
    match_ids = bind_rows(match_ids_non_op, match_ids_op)
    
    df_match_op = filter(df, id_hospitalization %in% match_ids_op$id_hospitalization)
    
    df_match_non_op = df %>% 
      filter(id_hospitalization %in% match_ids_non_op$id_hospitalization) %>%
      filter(facility_name %in% unique(df_match_op$facility_name)) 
    df_match = bind_rows(df_match_op, df_match_non_op)
    df_match = merge(df_match, match_ids, by = "id_hospitalization", all.x = T)
  }
  
  #Template matching model taking advantage of matched data
  tm_model =  glmer(mort30 ~ pred_mort +  facility_name + (1 | id_match), family = binomial, 
                    data = df_match, 
                    control = glmerControl(optimizer = "bobyqa"), 
                    nAGQ = 0)
  tm_facility_coefficients = fixef(tm_model)[3:length(fixef(tm_model))] 
  #adding column for facility names
  tm_facility_coefficients = as.data.frame(cbind(tm_facility_coefficients,attr(tm_facility_coefficients,"names")))
  rownames(tm_facility_coefficients) = NULL
  colnames(tm_facility_coefficients) = c("tm_coefficient","facility")
  tm_facility_coefficients$facility = as.character(tm_facility_coefficients$facility)
  tm_facility_coefficients$tm_coefficient = as.numeric(as.character(tm_facility_coefficients$tm_coefficient))
  #adding in reference hospital with coefficient of zero and sorting
  tm_facility_coefficients = rbind(tm_facility_coefficients, c(0,"facility_1") ) 
  tm_facility_coefficients$tm_coefficient = as.numeric(as.character(tm_facility_coefficients$tm_coefficient))

  #### Compare to non-matched standard hierarchical regression
  #first, filtering dataset to exclude hospitals with insufficient operative cases
  df_level = df %>% 
    filter(icu_level==level) %>%
    filter(facility_name %in% unique(df_match$facility_name))
  
  #Standard hierarchical model using entire dataset
  reg_model =  glmer(mort30 ~ pred_mort + (1 | facility_name), family = binomial(link = "logit"), 
                     data= df_level, 
                     control = glmerControl(optimizer = "bobyqa"), 
                     nAGQ = 0)
  
  reg_facility_coefficients = ranef(reg_model)$facility_name
  reg_facility_coefficients$facility = rownames(reg_facility_coefficients)
  rownames(reg_facility_coefficients) = NULL
  colnames(reg_facility_coefficients) = c("reg_coefficient","facility")  
  reg_facility_coefficients$facility = as.character(reg_facility_coefficients$facility)
  reg_facility_coefficients$reg_coefficient = as.numeric(as.character(reg_facility_coefficients$reg_coefficient))
  
  combined_coefficients = full_join(tm_facility_coefficients, reg_facility_coefficients, by = "facility")
  combined_coefficients = mutate(combined_coefficients,
                                 #reg_odds_ratio = exp(reg_coefficient),
                                 #tm_odds_ratio = exp(tm_coefficient),
                                 reg_rank = dense_rank(reg_coefficient),
                                 tm_rank = dense_rank(tm_coefficient),
                                 reg_quintile = as.numeric(cut_number(reg_coefficient, n = 5)),
                                 tm_quintile = as.numeric(cut_number(tm_coefficient, n = 5)),
                                 reg_quintile_categorized = if_else(reg_quintile == 1, "Quintile 1",
                                                                    if_else(reg_quintile == 5, "Quintile 5",
                                                                            "Quintile 2-4")),
                                 reg_quintile_categorized = factor(reg_quintile_categorized,
                                                                   levels = c("Quintile 1", "Quintile 2-4", "Quintile 5")),
                                 tm_quintile_categorized = if_else(tm_quintile == 1, "Quintile 1",
                                                                   if_else(tm_quintile == 5, "Quintile 5",
                                                                           "Quintile 2-4")),
                                 tm_quintile_categorized = factor(tm_quintile_categorized,
                                                                  levels = c("Quintile 1", "Quintile 2-4", "Quintile 5"))
                                 # agreement = if_else(reg_quintile==5 & tm_quintile==5, "Both highest quintile",
                                 #                     if_else(reg_quintile==5, "Regression highest quintile",
                                 #                             if_else(tm_quintile==5,"TM highest quintile",
                                 #                                     "Neither")))
  )   
  comparison = ggplot(combined_coefficients, aes(y=reg_rank, x = tm_rank)) +
    geom_point() +
    ggtitle(paste("30-day Mortality in ",tier, " Hospitals", sep = "")) +
    theme(plot.title = element_text(hjust = 0.5))+
    ylab("Rank from standard model") +
    xlab("Rank from modular template matching")
  
  # ggsave(filename = paste("Matching results/Charts and tables/",run,level," template matching vs standard regression rank v2.png", sep = ""), 
  #        height = 4, width = 4)
  # 
  quintile_tab = as.data.frame.matrix(table(combined_coefficients$reg_quintile_categorized, 
                                            combined_coefficients$tm_quintile_categorized))
  
  quintile_tab2 = tableGrob(quintile_tab,
                            rows = c("Std Quintile 1","Std Quintile 2-4", "Std Quintile 5"),
                            cols = c("TM Quintile 1","TM Quintile 2-4", "TM Quintile 5"),
                            theme = ttheme_default(base_size = 11,
                                                   colhead = list(fg_params=list(fontface=3)),
                                                   rowhead = list(fg_params=list(fontface=3))))  
  quintile_tab_fisher = as.table(rbind(c(sum(quintile_tab[1]),sum(quintile_tab[,1]) ),
                                       c(sum(quintile_tab[2]),sum(quintile_tab[,2]) ),
                                       c(sum(quintile_tab[3]),sum(quintile_tab[,3]) )))
  print(paste("Fisher test results:", tier))
  print(fisher.test(quintile_tab))
  #set layout for combining figure with table
  lay = rbind(c(1,1,1), c(1,1,1),c(1,1,1),c(1,1,1),c(1,1,1),c(1,1,1), c(1,1,1), c(1,1,1),c(1,1,1),c(1,1,1),
              c(1,1,1), c(1,1,1),c(1,1,1),c(1,1,1),
              #c(3,3,3), 
              c(2,2,2), c(2,2,2), c(2,2,2),c(2,2,2),c(2,2,2) )
  combined_plot = grid.arrange(comparison, quintile_tab2, layout_matrix = lay)
  
  ggsave(filename = paste("Matching results/Charts and tables/",run,level,
                          " _reg comparison.png", 
                          sep = ""),
         plot = combined_plot,
         height =5.8, width = 4.7)
}


