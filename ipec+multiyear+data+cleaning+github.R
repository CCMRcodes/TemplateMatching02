################################################################################
#
# Template Matching Data Cleaning
#
# Date: 2019-06-21
# Author: Daniel Molling <daniel.molling@va.gov>
#
################################################################################

#### loading  packages  #####
library(lattice)  
library(survival) 
library(Formula)  
library(ggplot2)  
library(Hmisc)
library(haven)
library(tidyr)
library(dplyr)
library(lubridate)
library(here)


##### Loading raw multiyear Data, save seperate years #####
df = read_sas("IPEC data/ip.sas7bdat") #orignal raw file 

#renaming sodium variable to avoid namespace issues
names(df)[names(df) == 'NA'] <- 'na' 
df$admit_year = year(df$ADMDATE)
 
df_2017 = filter(df, admit_year==2017)

saveRDS(df_2017, file = "IPEC data/ipec2017_raw.rds")

##### 2017 Filtering hospitalizations that don't meet inclusion criteria #####
df_2017 = readRDS("IPEC data/ipec2017_raw.rds")

df_2017 = filter(df_2017, 
                 OCCURRENCE_ALL==1 #keeping only first unit stay for a hospitalization
                 ) 
#After: 671897 hospitalizations at 180 facilities

##### Keeping only hospitals that meet the VAPD definition of hospital for 2013-2017 #####
VAPD_hosp_list = read.csv("IPEC data/unique_sta6a_2017_20190523.csv", header = T)
VAPD_hosp_list = select(VAPD_hosp_list,
                        facility_official_name = Official__Station_Name,
                        sta6a = sta6a)
### Merging in Official Station Names ###
facility_info = read.csv("IPEC data/facility_information.csv") 
facility_info = select(facility_info, 
                       SITE2 = Location.Descriptive.Name..Common.Name.,
                       facility_official_name = Official.Station.Name)
facility_info$SITE2 = as.character(facility_info$SITE2)
facility_info$facility_official_name = as.character(facility_info$facility_official_name)
df_2017 = merge(df_2017, facility_info, by = "SITE2", all.x = T)
#Correcting official facility names
df_2017 = mutate(df_2017,
            facility_official_name = if_else(SITE2=="Columbia-MIssouri",
                                            "Harry S. Truman Memorial Veterens' Hospital" ,facility_official_name),
            facility_official_name = if_else(SITE2=="Denver-Colorado",
                                            "Denver VA Medical Center" ,facility_official_name))

#recoding official station names to match VAPD names
df_2017 = mutate(df_2017,
                 facility_official_name = recode(facility_official_name,
                                                 "Aleda E. Lutz Department of Veterans Affairs Medical Center" = "Aleda E. Lutz VA Medical Center-Saginaw",
                                                 "Alvin C. York Veterans' Administration Medical Center" = "Alvin C. York VA Medical Center-Murfreesboro",
                                                 "Audie L. Murphy Memorial Veterans' Hospital" = "Audie L. Murphy VA Medical Center-San Antonio",
                                                 "Augusta VA Medical Center-Uptown" = "Charlie Norwood VA Medical Center-Augusta",
                                                 "Baltimore VA Medical Center" = "Maryland VA Medical Center-Baltimore",
                                                 "Biloxi VA Medical Center" = "Gulf Coast VA Medical Center-Biloxi",
                                                 "Bob Stump Department of Veterans Affairs Medical Center" = "Northern Arizona VA Medical Center-Prescott",
                                                 "Brockton VA Medical Center" = "Boston VA Medical Center-Brockton",
                                                 "Brooklyn VA Medical Center" = "New York Harbor VA Medical Center-Brooklyn",
                                                 "Bruce W. Carter Department of Veterans Affairs Medical Center" = "Bruce W. Carter VA Medical Center-Miami",
                                                 "Buffalo VA Medical Center" = "Western New York VA Medical Center-Buffalo",
                                                 "C.W. Bill Young Department of Veterans Affairs Medical Center" = "C.W. Bill Young VA Medical Center-Bay Pines",
                                                 "Captain James A. Lovell Federal Health Care Center" = "Captain James A. Lovell VA Medical Center-North Chicago",
                                                 "Carl T. Hayden Veterans' Administration Medical Center" = "Carl T. Hayden VA Medical Center-Phoenix",
                                                 "Carl Vinson Veterans' Administration Medical Center" = "Carl Vinson VA Medical Center-Dublin",
                                                 "Castle Point VA Medical Center" = "Hudson Valley VA Medical Center-Castle Point",
                                                 "Central Alabama VA Medical Center-Tuskegee" = "Central Alabama VA Medical Center-Montgomery",
                                                 "Charles George Department of Veterans Affairs Medical Center" = "Charles George VA Medical Center-Asheville",
                                                 "Charlie Norwood Department of Veterans Affairs Medical Center" = "Charlie Norwood VA Medical Center-Augusta",
                                                 "Clement J. Zablocki Veterans' Administration Medical Center" = "Clement J. Zablocki VA Medical Center-Milwaukee",
                                                 "Colmery-O'Neil Veterans' Administration Medical Center" = "Colmery-O'Neil VA Medical Center-Topeka",
                                                 "Corporal Michael J. Crescenz Department of Veterans Affairs Medical Center" = "Philadelphia VA Medical Center",
                                                 "Dallas VA Medical Center" = "North Texas VA Medical Center-Dallas",
                                                 "Danville VA Medical Center" = "Illiana VA Medical Center-Danville",
                                                 "Denver VA Medical Center" = "Eastern Colorado VA Medical Center-Denver",
                                                 "Des Moines VA Medical Center" = "Central Iowa VA Medical Center-Des Moines",
                                                 "Doris Miller Department of Veterans Affairs Medical Center" = "",
                                                 "Dwight D. Eisenhower Department of Veterans Affairs Medical Center" = "Dwight D. Eisenhower VA Medical Center-Leavenworth",
                                                 "East Orange VA Medical Center" = "New Jersey VA Medical Center-East Orange",
                                                 "Edward Hines Junior Hospital" = "Edward Hines Jr. VA Medical Center-Hines",
                                                 "Eugene J. Towbin Healthcare Center" = "Eugene J. Tobin VA Medical Center-North Little Rock",
                                                 "Fort Harrison VA Medical Center" = "Montana VA Medical Center-Fort Harrison",
                                                 "Fort Meade VA Medical Center" = "Black Hills VA Medical Center-Fort Meade",
                                                 "Fort Wayne VA Medical Center" = "Northern Indiana VA Medical Center-Fort Wayne",
                                                 "Franklin Delano Roosevelt Hospital" = "",
                                                 "Fresno VA Medical Center" = "Central California VA Medical Center-Fresno",
                                                 "G.V. (Sonny) Montgomery Department of Veterans Affairs Medical Center" = "G.V. (Sonny) Montgomery VA Medical Center-Jackson",
                                                 "George E. Wahlen Department of Veterans Affairs Medical Center" = "George E. Wahlen VA Medical Center-Salt Lake City",
                                                 "Harry S. Truman Memorial Veterens' Hospital" = "Harry S. Truman VA Medical Center-Columbia",
                                                 "Hershel 'Woody' Williams VA Medical Center" = "Huntington VA Medical Center",
                                                 "Hot Springs VA Medical Center" = "Black Hills VA Medical Center-Hot Springs",
                                                 "Hunter Holmes McGuire Hospital" = "Hunter Holmes McGuire VA Medical Center-Richmond",
                                                 "Ioannis A. Lougaris Veterans' Administration Medical Center" = "Ioannis A. Lougaris VA Medical Center-Sierra Nevada Reno",
                                                 "Jack C. Montgomery Department of Veterans Affairs Medical Center" = "Jack C. Montgomery VA Medical Center-Muskogee",
                                                 "James A. Haley Veterans' Hospital" = "James A. Haley VA Medical Center-Tampa",
                                                 "James E. Van Zandt Veterans' Administration Medical Center" = "James E. Van Zandt VA Medical Center-Altoona",
                                                 "James H. Quillen Department of Veterans Affairs Medical Center" = "James H. Quillen VA Medical Center-Mountain Home",
                                                 "James J. Peters Department of Veterans Affairs Medical Center" = "James J. Peters VA Medical Center-Bronx",
                                                 "Jerry L. Pettis Memorial Veterans' Hospital" = "Loma Linda VA Medical Center",
                                                 "Jesse Brown Department of Veterans Affairs Medical Center" = "Jesse Brown VA Medical Center-Chicago",
                                                 "John Cochran Veterans Hospital" = "John Cochran VA Medical Center-St. Louis",
                                                 "John D. Dingell Department of Veterans Affairs Medical Center" = "John D. Dingell VA Medical Center-Detroit",
                                                 "John J. Pershing Veterans' Administration Medical Center" = "John J. Pershing VA Medical Center-Poplar Bluff",
                                                 "John L. McClellan Memorial Veterans' Hospital" = "John L. McClellan VA Medical Center-Little Rock",
                                                 "Lake City VA Clinic" = "Lake City VA Medical Center",
                                                 "Louis A. Johnson Veterans' Administration Medical Center" = "Louis A. Johnson VA Medical Center-Clarksburg",
                                                 "Louis Stokes Cleveland Department of Veterans Affairs Medical Center" = "Louis Stokes VA Medical Center-Cleveland",
                                                 "Lyons VA Medical Center" = "",
                                                 "Malcom Randall Department of Veterans Affairs Medical Center" = "Malcom Randall VA Medical Center-Gainesville",
                                                 "Manhattan VA Medical Center" = "New York Harbor VA Medical Center-Manhattan",
                                                 "Mann-Grandstaff Department of Veterans Affairs Medical Center" = "Mann-Grandstaff VA Medical Center-Spokane",
                                                 "Michael E. DeBakey Department of Veterans Affairs Medical Center" = "Michael E. DeBakey VA Medical Center-Houston",
                                                 "New Orleans VA Medical Center" = "Southeast Louisiana VA Medical Center-New Orleans",
                                                 "North Las Vegas VA Medical Center" = "Southern Nevada VA Medical Center-Las Vegas",
                                                 "Olin E. Teague Veterans' Center" = "Olin E. Teague VA Medical Center-Temple",
                                                 "Omaha VA Medical Center" = "Nebraska-Western Iowa VA Medical Center-Omaha",
                                                 "Oscar G. Johnson Department of Veterans Affairs Medical Facility" = "Oscar G. Johnson VA Medical Center-Iron Mountain",
                                                 "Overton Brooks Veterans' Administration Medical Center" = "Overton Brooks VA Medical Center-Shreveport",
                                                 "Perry Point VA Medical Center" = "Maryland VA Medical Center-Perry Point",
                                                 "Philadelphia VA Domiciliary" = "",
                                                 "Ralph H. Johnson Department of Veterans Affairs Medical Center" = "Ralph H. Johnson VA Medical Center-Charleston",
                                                 "Raymond G. Murphy Department of Veterans Affairs Medical Center" = "Raymond G. Murphy VA Medical Center-Albuquerque",
                                                 "Richard L. Roudebush Veterans' Administration Medical Center" = "Richard L. Roudebush VA Medical Center-Indianapolis",
                                                 "Robert J. Dole Department of Veterans Affairs Medical and Regional Office Center" = "Robert J. Dole VA Medical Center-Wichita",
                                                 "Robley Rex Department of Veterans Affairs Medical Center" = "Robley Rex VA Medical Center-Louisville",
                                                 "Royal C. Johnson Veterans' Memorial Hospital" = "Sioux Falls VA Medical Center",
                                                 "Sacramento VA Medical Center" = "Northern California VA Medical Center-Sacramento",
                                                 "Samuel S. Stratton Department of Veterans Affairs Medical Center" = "Samuel S. Stratton VA Medical Center-Albany",
                                                 "Seattle VA Medical Center" = "Puget Sound VA Medical Center-Seattle",
                                                 "St. Cloud VA Medical Center" = "",
                                                 "Thomas E. Creek Department of Veterans Affairs Medical Center" = "Thomas E. Creek VA Medical Center-Amarillo",
                                                 "Togus VA Medical Center" = "Maine VA Medical Center",
                                                 "Troy Bowling Campus" = "Lexington VA Medical Center-Cooper",
                                                 "Tucson VA Medical Center" = "Southern Arizona VA Medical Center-Tucson",
                                                 "W.G. (Bill) Hefner Salisbury Department of Veterans Affairs Medical Center" = "W.G. (Bill) Hefner VA Medical Center-Salisbury",
                                                 "West Haven VA Medical Center" = "Connecticut VA Medical Center-West Haven",
                                                 "West Los Angeles VA Medical Center" = "Greater Los Angeles VA Medical Center",
                                                 "West Roxbury VA Medical Center" = "Boston VA Medical Center-West Roxbury",
                                                 "William S. Middleton Memorial Veterans' Hospital" = "William S. Middleton Memorial Veterans Medical Center-Madison",
                                                 "Wm. Jennings Bryan Dorn Department of Veterans Affairs Medical Center" = "Wm. Jennings Bryan Dorn VA Medical Center-Columbia"
                                                 ))

### Filtering to match VAPD hospital list
df_2017 = filter(df_2017, facility_official_name %in% VAPD_hosp_list$facility_official_name)
#removing 3 additional hospitals with duplicate official facility names that aren't in the VAPD
df_2017 = filter(df_2017, !SITE2 %in% c("Tuskegee","Augusta Uptown","Marion-Indiana"))
#After: 668592 hospitalizations at 134 hospitals


df_2017 = filter(df_2017,ADMISSION_SOURCE != "VAH") #exclude transfers from another VA hospital
#After: 664941 hospitalizations at 134 hospitals

df_2017 = filter(df_2017, UNITDX != "") #dropping 64 hospitalizations with missing diagnosis codes
#After: 664874 hospitalizations at 134 hospitals

##### Merging in HPUC CCS group diagnosis and procedure codes to 2017 data #####
#formatting codes for crosswalk
df_2017$unit_pr = gsub("[^0-9a-zA-Z]","",df_2017$PROC_ICD9)
df_2017$unit_dx1 = gsub("[^0-9a-zA-Z]","",df_2017$UNITDX)

#load and format file containing crosswalk from ICD diagnosis codes to HCUP Clinical Classification Scheme group codes
DX_icd10 = read.csv("IPEC data/dxref_2018.csv") 
DX_icd10$unit_dx1          = gsub("[^0-9a-zA-Z]","",DX_icd10$X.ICD.10.CM.CODE.)
DX_icd10$unit_dx_ccs       = gsub("[^0-9a-zA-Z]","",DX_icd10$X.CCS.CATEGORY.)
DX_icd10$unit_dx_multiccs1 = gsub("[^0-9a-zA-Z]","",DX_icd10$X.MULTI.CCS.LVL.1.)
DX_icd10 = select(DX_icd10, unit_dx1, unit_dx_ccs, unit_dx_ccs_description = X.CCS.CATEGORY.DESCRIPTION., 
                  unit_dx_multiccs1, unit_dx_multiccs1_description = X.MULTI.CCS.LVL.1.LABEL.)

#merge ccs codes and descriptions into dataset
df_2017 = merge(df_2017, DX_icd10, by = "unit_dx1", all.x = T)

#load and format file containing crosswalk from ICD procedure codes to HCUP Clinical Classification Scheme group codes
PR_icd10 = read.csv("IPEC data/prref_2017.csv") 
PR_icd10$unit_pr     = gsub("[^0-9a-zA-Z]","",PR_icd10$X.ICD.10.PCS.CODE.)
PR_icd10$pr_unit_ccs = gsub("[^0-9a-zA-Z]","",PR_icd10$X.CCS.CATEGORY.)
PR_icd10 = select(PR_icd10, unit_pr, pr_unit_ccs)

df_2017 = merge(df_2017, PR_icd10, by = "unit_pr", all.x = T)


##### Exclude transfers and diagnosis and facilities with low patient count,  #####
#exclude 167 hospitalizations for transplant operations
df_2017 = filter(df_2017, !DRG %in% c(1,2,5,6,7,8,10,14,15)) 
#After: 664707 hospitalizations at 134 hospitals

#add variable for diagnosis count for each CCS group
df_2017 = group_by(df_2017, unit_dx_ccs) %>% mutate(unit_dx_ccs_count = n()) 
#dropping patients with rare diagnosis (less than .33% of hospitalizations)
df_2017 = filter(df_2017, unit_dx_ccs_count/nrow(df_2017) > .0033)
#After: 563716 hospitalizations at 134 hospitals

#add variable for hospitalizatin count by site
df_2017 = group_by(df_2017, SITE) %>% mutate(site_hospitalizations = n()) 
#exclude hospitals with fewer than 900 patients
df_2017 = filter(df_2017, site_hospitalizations >= 900) 
#After: 556857 hospitalizations at 121 facilities

##### creating indicator and directional variables for mortality covariates  ##### 
df_2017 = mutate(df_2017, 
                 female   = if_else(GENDER      == "F",1,0),
                 black    = if_else(RACE_RECODE == "African American/Black",1,0),
                 hispanic = if_else(ETHNICITY   == "HISPANIC OR LATINO",1,0),
                 white    = if_else(RACE_RECODE == "Caucasian/White",1,0))

##### Assign standard score (0) to patients whose lab score was missing #####
#    Those with missing labs were flagged in the IPEC data and assigned
#    a standard value for that lab. 
#    Note: GFR is not scored, so it is used later as a continuous variable. We
#    truncate GFR at a maximum of 120
#    Note: PH is not scored by IPEC. We apply scores from "Anzics core APD 
#    programmer's data dictionary" version 5.6 page 94
df_2017 = mutate(df_2017, 
                 WBC_SC     = if_else(WBC_MISS == 1, 0, WBC_SC),
                 ALBVAL_SC  = if_else(ALBVAL_MISS == 1, 0, ALBVAL_SC),
                 BILI_SC    = if_else(BILI_MISS == 1, 0, BILI_SC),
                 BUN_SC     = if_else(BUN_MISS == 1, 0, BUN_SC),
                 GLUCOSE_SC = if_else(GLUCOSE_MISS == 1, 0, GLUCOSE_SC),
                 HCT_SC     = if_else(HCT_MISS == 1, 0, HCT_SC),
                 NA_SC      = if_else(NA_MISS == 1, 0, NA_SC),
                 PAO2_SC    = if_else(PAO2_MISS == 1, 0, PAO2_SC),
                 PCO2_SC    = if_else(PCO2_MISS == 1, 0, PCO2_SC),
                 ph_unit_score      = case_when(PH < 7.2 & PCO2 < 50 ~ 12,
                                        PH < 7.2 & PCO2 >= 50 ~ 4,
                                        PH < 7.35 & PH >= 7.20 & PCO2 < 30 ~ 9,
                                        PH < 7.30 & PH >= 7.20 & PCO2 < 40 & PCO2 >= 30 ~ 6,
                                        PH < 7.30 & PH >= 7.20 & PCO2 < 50 & PCO2 >= 40 ~ 3,
                                        PH < 7.30 & PH >= 7.20 & PCO2 >= 50 ~ 2,
                                        PH < 7.50 & PH >= 7.35 & PCO2 < 30 ~ 5,
                                        PH < 7.45 & PH >= 7.30 & PCO2 < 45 & PCO2 >= 30 ~ 0,
                                        PH < 7.45 & PH >= 7.30 & PCO2 >= 45 ~ 1,
                                        PH < 7.50 & PH >= 7.45 & PCO2 < 35 & PCO2 >= 30 ~ 0,
                                        PH < 7.50 & PH >= 7.45 & PCO2 < 45 & PCO2 >= 35 ~ 2,
                                        PH < 7.60 & PH >= 7.50 & PCO2 < 40 ~ 3,
                                        PH >= 7.6 & PCO2 < 25 ~ 0,
                                        PH >= 7.6 & PCO2 < 40 & PCO2 >= 25 ~ 3,
                                        PH >= 7.5 & PCO2 >= 40 ~ 12,
                                        PH < 7.50 & PH >= 7.45 & PCO2 >= 45 ~ 12
                                        ),
                 ph_unit_score      = if_else(PH_MISS == 1, 0, ph_unit_score),
                 gfr_unit_score     = as.numeric(if_else(GFR > 120, 120, GFR))
                 )

#####creating indicators for top 20 group diagnosis codes based on single-level CCS #####
df_2017 = mutate(df_2017, 
                 dx_alchohol                = if_else(unit_dx_ccs == 660, 1, 0),
                 dx_mood                    = if_else(unit_dx_ccs == 657, 1, 0),
                 dx_chf_nonhp               = if_else(unit_dx_ccs == 108, 1, 0),
                 dx_chest_pain              = if_else(unit_dx_ccs == 102, 1, 0),
                 dx_dysrhythmia             = if_else(unit_dx_ccs == 106, 1, 0),
                 dx_coron_athero            = if_else(unit_dx_ccs == 101, 1, 0),
                 dx_septicemia              = if_else(unit_dx_ccs ==   2, 1, 0),
                 dx_copd                    = if_else(unit_dx_ccs == 127, 1, 0),
                 dx_htn                     = if_else(unit_dx_ccs ==  99, 1, 0),
                 dx_substance               = if_else(unit_dx_ccs == 661, 1, 0),
                 dx_pneumonia               = if_else(unit_dx_ccs == 122, 1, 0),
                 dx_skin_infection          = if_else(unit_dx_ccs == 197, 1, 0),
                 dx_osteoarthritis          = if_else(unit_dx_ccs == 203, 1, 0),
                 dx_psychotic               = if_else(unit_dx_ccs == 659, 1, 0),
                 dx_resp_fail_adult         = if_else(unit_dx_ccs == 131, 1, 0),
                 dx_device_complication     = if_else(unit_dx_ccs == 237, 1, 0),
                 dx_diabetes_complication   = if_else(unit_dx_ccs ==  50, 1, 0),
                 dx_procedure_complication  = if_else(unit_dx_ccs == 238, 1, 0),
                 dx_acute_renal_fail        = if_else(unit_dx_ccs == 157, 1, 0),
                 dx_uti                     = if_else(unit_dx_ccs == 159, 1, 0))

##### Create variable for unit diagnosis grouped into 7 major categories #####
psychiatric_codes = c("660","657","661","659","651","662","650")
renal_codes = c("157","55","164","160","29","163","161")
gastrointestinal_codes = c("153","149","155","152","143","145","251","151","16","138","154","58")
infection_codes = c(2,"122","159","197","135","201","123","146")

cardiovascular_codes = c("108","102","106","101","99","100","96","109","117","115","103","112","118","110","114")
respiratory_codes = c("127","133","19","129","130","131")
unclassified_codes = c("203","237","50","238","205","95","245","59","653","211","259","257","83","42","47","45","252")

df_2017 = mutate(df_2017, 
                 unit_dx_grouped = case_when(unit_dx_ccs %in% psychiatric_codes
                                                ~ "Psychiatric and Substance Abuse",
                                             unit_dx_ccs %in% renal_codes 
                                                ~ "Renal / Genitorinary",
                                             unit_dx_ccs %in% gastrointestinal_codes
                                                ~ "Gastrointestinal",
                                             unit_dx_ccs %in% infection_codes
                                                ~ "Infection",
                                             unit_dx_ccs %in% cardiovascular_codes
                                                ~ "Cardiovascular",
                                             unit_dx_ccs %in% respiratory_codes
                                                ~ "Respiratory",
                                             unit_dx_ccs %in% unclassified_codes
                                                ~ "Unclassified / Residual")
                   )
#Create indicators for each of these groups
df_2017 = mutate(df_2017, 
                 dxgroup_psych = if_else(unit_dx_grouped == "Psychiatric and Substance Abuse",1,0),
                 dxgroup_renal = if_else(unit_dx_grouped == "Renal / Genitorinary",1,0),
                 dxgroup_gastro = if_else( unit_dx_grouped == "Gastrointestinal",1,0),
                 dxgroup_infection = if_else(unit_dx_grouped == "Infection",1,0),
                 dxgroup_cardiovascular = if_else(unit_dx_grouped == "Cardiovascular",1,0),
                 dxgroup_resp = if_else(unit_dx_grouped == "Respiratory",1,0),
                 dxgroup_other = if_else(unit_dx_grouped == "Unclassified / Residual",1,0))

#Create indicators for admission source either Nursing Home or Emergency Department                   
df_2017 = mutate(df_2017, 
                 admit_nursing_home = if_else(ADMISSION_SOURCE == "VANH" | ADMISSION_SOURCE == "NON-VANH", 1, 0),
                 admit_emergency_dept = if_else(ADMISSION_SOURCE == "VA ED", 1, 0))

#exclude hospitals with greater than 90% psych-related diagnosis admissions
df_2017 = group_by(df_2017, SITE) %>% 
          mutate(psych_percent = sum(dxgroup_psych)/n()) 
df_2017 = filter(df_2017, psych_percent <= .9) 
#After: 550573 hospitalizations at 117 facilities

##### Create 30-day all-cause mortality indicator #####
df_2017 = mutate(df_2017, 
                 mort30 = if_else(DEATHDATE - ADMDATE <= 30, 1,0),
                 mort30 = if_else(is.na(DEATHDATE),0,mort30))

##### Merging in ICU level data #####
### Note: In the paper text we refer to this variable as "Hospital Tier" or "Tier"
icu_levels = read.csv("IPEC data/icu_level.csv")
icu_levels$facility_id = icu_levels$STA5A_parent_number
icu_levels$icu_score = as.factor(icu_levels$icu_score)
icu_levels = select(icu_levels, SITE = facility_id, icu_level, icu_score)
#updating to addinew rows to correct for absences in file
#updates come from email correspondence with VA staff familiar with this data
icu_levels = bind_rows(icu_levels, 
                       c(SITE = "523A4", icu_level = "Level 1", icu_score = "5"),
                       c(SITE = "596A4", icu_level = "Level 1", icu_score = "5"),
                       c(SITE = "612A4", icu_level = "Level 3", icu_score = "3"),
                       c(SITE = "626A4", icu_level = "Level 1", icu_score = "5"),
                       c(SITE = "630A4", icu_level = "Level 1", icu_score = "5"),
                       c(SITE = "523A5", icu_level = "Level 1", icu_score = "5"),
                       c(SITE = "573A4", icu_level = "Level 1", icu_score = "5"),
                       c(SITE = "657A0", icu_level = "Level 1", icu_score = "5"),
                       c(SITE = "674A4", icu_level = "Level 2", icu_score = "4"),
                       c(SITE = "610A4", icu_level = "Level 4", icu_score = "2"),
                       c(SITE = "598A0", icu_level = "Level 1", icu_score = "5"),
                       c(SITE = "589A6", icu_level = "Level 3", icu_score = "3")
)
df_2017 = left_join(df_2017, icu_levels, by = "SITE")
rm(icu_levels)
##### filtering unneeded variables and renaming #####
 df_2017 = select(df_2017,
               admit_date         = ADMDATE,
               admission_source         = ADMISSION_SOURCE,
               admit_time         = ADMTIME,
               age         = AGE,
               albumin_unit         = ALBVAL,
               albumin_unit_miss         = ALBVAL_MISS,
               albumin_unit_score         = ALBVAL_SC,
               bilirubin_unit         = BILI,
               bilirubin_unit_miss         = BILI_MISS,
               bilirubin_unit_score         = BILI_SC,
               bun_unit         = BUN,
               bun_unit_miss         = BUN_MISS,
               bun_unit_score         = BUN_SC,
               cholesterol_unit         = CHOLESTEROL,
               cholesterol_unit_miss         = CHOLESTEROL_MISS,
               cholesterol_unit_score         = CHOLESTEROL_SC,
               complexity_level         = COMPLEXITY_LEVEL,
               cpk_unit         = CPK,
               cpk_unit_miss         = CPK_MISS,
               cpk_unit_score         = CPK_SC,
               cr_unit         = CR,
               creatinine_unit_miss         = CR_MISS,
               creatnine_unit_score         = CR_SC,
               chf         = CX1,
               dm_comp         = CX10,
               dm_uncomp         = CX11,
               hypothyroid         = CX12,
               renal         = CX13,
               liver         = CX14,
               pud         = CX15,
               immunedef         = CX16,
               lymphoma         = CX17,
               cancer_met         = CX18,
               cancer_nonmet         = CX19,
               ra         = CX20,
               coag         = CX21,
               obesity         = CX22,
               wtloss         = CX23,
               fen         = CX24,
               anemia_cbl         = CX25,
               anemia_def         = CX26,
               etoh         = CX27,
               drug         = CX28,
               psychoses         = CX29,
               valvular_d2         = CX3,
               depression         = CX30,
               pulm_circ         = CX4,
               pvd         = CX5,
               htn         = CX6,
               paralysis         = CX7,
               neuro         = CX8,
               pulm         = CX9,
               inhosp_mort         = DEADDIS,
               death_date         = DEATHDATE,
               death_time         = DEATHTIME,
               discharge_date         = DISDATE,
               disposition_place         = DISPLACE,
               district         = DISTRICT,
               date_of_birth         = DOB,
               drg         = DRG,
               dx_unit_description         = DX_ICD9_DESCRIPTION,
               fio2_unit         = FIO2,
               gfr_unit         = GFR,
               glucose_unit         = GLUCOSE,
               glucose_unit_miss         = GLUCOSE_MISS,
               glucose_unit_score         = GLUCOSE_SC,
               htc_unit         = HCT,
               hct_unit_miss         = HCT_MISS,
               hct_unit_score         = HCT_SC,
               id_patient         = ID,
               id_hospitalization         = ID2,
               INPATIENTSID         = INPATIENTSID,
               facility_level         = LEVEL,
               sodium_unit         = na,
               sodium_unit_miss         = NA_MISS,
               sodium_unit_score         = NA_SC,
               unit_stay_order         = OCCURRENCE_ALL,
               operative         = OPERATIVE,
               or_date         = ORDATE,
               p02_unit         = PAO2,
               p02_unit_miss         = PAO2_MISS,
               pao2_unit_score         = PAO2_SC,
               PATIENTSID         = PATIENTSID,
               pco2_unit         = PCO2,
               pco2_unit_miss         = PCO2_MISS,
               pco2_unit_score         = PCO2_SC,
               ph_unit         = PH,
               ph_unit_miss         = PH_MISS,
               platelets_unit         = PLATELETS,
               platelets_unit_miss         = PLATELETS_MISS,
               platelets_unit_score         = PLATELETS_SC,
               pr_unit_description         = PR_ICD9_DESCRIPTION,
               pr_unit         = PROC_ICD9,
               surgery         = PROCEDURE,
               race_cat         = RACE_RECODE,
               ssn_scrambled         = SCRSSN,
               facility_id         = SITE,
               facility_name         = SITE2,
               sta3n_unit         = STA3N,
               unit_dx1         = UNITDX,
               unit_dx10         = UNITDX10,
               unit_dx11         = UNITDX11,
               unit_dx12         = UNITDX12,
               unit_dx13         = UNITDX13,
               unit_dx14         = UNITDX14,
               unit_dx15         = UNITDX15,
               unit_dx16         = UNITDX16,
               unit_dx17         = UNITDX17,
               unit_dx18         = UNITDX18,
               unit_dx19         = UNITDX19,
               unit_dx2         = UNITDX2,
               unit_dx20         = UNITDX20,
               unit_dx21         = UNITDX21,
               unit_dx22         = UNITDX22,
               unit_dx23         = UNITDX23,
               unit_dx24         = UNITDX24,
               unit_dx25         = UNITDX25,
               unit_dx3         = UNITDX3,
               unit_dx4         = UNITDX4,
               unit_dx5         = UNITDX5,
               unit_dx6         = UNITDX6,
               unit_dx7         = UNITDX7,
               unit_dx8         = UNITDX8,
               unit_dx9         = UNITDX9,
               unit_specialty         = UNITSPEC,
               unitspec_rpt     = UNITSPEC_RPT,
               visn         = VISN,
               wbc_unit         = WBC,
               wbc_unit_miss         = WBC_MISS,
               wbc_unit_score         = WBC_SC,
               ph_unit_score,
               white,
               black,
               dx_alchohol,              
               dx_mood,                  
               dx_chf_nonhp,             
               dx_chest_pain,            
               dx_dysrhythmia,           
               dx_coron_athero,          
               dx_septicemia,            
               dx_copd,                  
               dx_htn,                   
               dx_substance,             
               dx_pneumonia,             
               dx_skin_infection,        
               dx_osteoarthritis,        
               dx_psychotic,             
               dx_resp_fail_adult,       
               dx_device_complication,   
               dx_diabetes_complication, 
               dx_procedure_complication,
               dx_acute_renal_fail,      
               dx_uti,
               female,
               gfr_unit_score,
               hispanic,
               ph_unit_score,
               mort30,
               unit_dx_ccs,
               unit_dx_ccs_description,
               unit_dx_multiccs1,
               unit_dx_multiccs1_description,
               pr_unit_ccs,
               unit_dx_grouped,
               icu_level,
               icu_score,
               admit_nursing_home,
               admit_emergency_dept,
               dxgroup_psych,
               dxgroup_renal,
               dxgroup_gastro,
               dxgroup_infection,
               dxgroup_cardiovascular,
               dxgroup_resp,
               dxgroup_other
               )  

##### Estimate predicted probability of 30-day mortality upon admission for each hospitalization #####
lab_scores = c("wbc_unit_score", "albumin_unit_score", "bilirubin_unit_score", "bun_unit_score", "glucose_unit_score", 
               "hct_unit_score", "sodium_unit_score", "pao2_unit_score", "pco2_unit_score", "ph_unit_score") 

#converting lab scores to character variables 
df_2017 = df_2017 %>% 
  as_tibble() %>% 
  mutate_at(lab_scores, as.character) 
  
demographics = c("female", "black", "white", "hispanic", "age")
top20_diagnosis = c("dx_alchohol", "dx_mood", "dx_chf_nonhp", "dx_chest_pain", "dx_dysrhythmia", "dx_coron_athero", "dx_septicemia", "dx_copd", "dx_htn", "dx_substance", "dx_pneumonia", "dx_skin_infection", "dx_osteoarthritis", "dx_psychotic", "dx_resp_fail_adult", "dx_device_complication", "dx_diabetes_complication", "dx_procedure_complication", "dx_acute_renal_fail", "dx_uti")

comorbidities = c("chf", "dm_comp", "dm_uncomp", "hypothyroid", "renal", "liver", "pud", "immunedef", "lymphoma", "cancer_met", "cancer_nonmet", "ra", "coag", "obesity", "wtloss", "fen", "anemia_cbl", "anemia_def", "etoh", "drug", "psychoses", "valvular_d2", "depression", "pulm_circ", "pvd", "htn", "paralysis", "neuro", "pulm")

admission_info = c("admit_nursing_home","admit_emergency_dept", "unitspec_rpt") 

mort30_formula = c(demographics, "operative", comorbidities, top20_diagnosis, admission_info) %>%
  paste(collapse = " + ") %>%
  paste("mort30 ~", . ,collapse = "") 
mort30_formula = as.formula(mort30_formula)

#logistic model
mort_model = glm(mort30_formula, data = df_2017, family = "binomial")
df_2017$pred_mort = predict(mort_model, type = "response") 

# Create categorical predicted mortality variables - 5 and 10 categories
df_2017$pred_mort_cat5 = cut2(df_2017$pred_mort,g = 5)
df_2017$pred_mort_cat10 = cut2(df_2017$pred_mort,g = 10)


##### Saving Cleaned Data #####
saveRDS(df_2017, file = "IPEC data/ipec2017_working.rds")
##### mortality model diagnostics #####
# mort_predictions = predict(mort_model, type = c("response"))
# library(pROC)
# roccurve = roc(df_2017$mort30 ~ mort_predictions)
# plot(roccurve)











