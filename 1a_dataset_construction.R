# Dataset construction
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/12/2024

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("here", "haven", "tidyverse", "janitor")

options(scipen = 9999)

# ----loading data----
source(here("scripts", "0_paths.R"))

aa_raw_data <-
  read_sas(paste0(path_to_raw_data,
                  "aa_adrd_cardiometabolic_tte.sas7bdat")) %>% clean_names()

bp_data <- read_sas(paste0(path_to_bp_data,
                           "c10049_bp.sas7bdat")) %>%
  clean_names()

# ----cleaning wide cardiometabolic data----

# n = 184,929 obs

wide_data <-
  
  aa_raw_data %>% filter(main_dem_v1_sample == 1) %>%
  
  # after filtering out prevalent dementia, n = 180,394
  
  filter(ethnicity_rev %in% c(1, 2, 3, 5, 9)) %>%
  
  # after restricting to South Asian, Chinese, Japanese, Filipino, White, n = 159,477
  
  select(
    subjid,
    
    # exposure
    htn_dx5yr_flag,
    htn_dx7yr_flag,
    htn_dx9yr_flag,
    
    # time vars
    survey_age,
    main_dem_v1_end_age,
    main_dem_v1_dem_age,
    main_dem_v1_fu_time,
    
    # outcome
    main_dem_v1_end_dem_flag,
    main_dem_v1_end_type,
    
    # covariates
    ethnicity_rev,
    asian,
    usaborn_rev,
    female,
    education_rev,
    edu_ge_college,
    
    # income related
    sizeofhh,
    income,
    income_pp,
    
    # additional covariates
    maritalstatus,
    smoking_status,
    usabornfather_rev,
    usabornmother_rev,
    
    # additional comorbidities to include in table 1
    diab_dx5yr_flag,
    diab_dx7yr_flag,
    diab_dx9yr_flag,
    first_prevstroke_flag,
    first_incstroke_flag,
    sr_depress,
    sr_bmi,
    
    # membership criteria
    presurv5yr_sample,
    presurv7yr_sample,
    presurv9yr_sample,
    
    # height variables
    sr_total_height_in,
    ehr_ht_median
  )

wide_data <- wide_data %>% 
  
  # deriving end age and end type
  mutate(
    end_age = case_when(
      main_dem_v1_end_type %in% c("ADMIN CENSORED", "END OF MEMBERSHIP") ~
        main_dem_v1_end_age,
      main_dem_v1_end_type %in% c("DEMENTIA", "DEATH", "CENSORED 90+") &
        main_dem_v1_end_age > 90 ~ 90,
      TRUE ~ main_dem_v1_end_age
    ),
    end_type = case_when(
      main_dem_v1_end_type %in% c("ADMIN CENSORED", "END OF MEMBERSHIP") ~
        str_to_title(main_dem_v1_end_type),
      main_dem_v1_end_type %in% c("DEMENTIA", "DEATH", "CENSORED 90+") &
        main_dem_v1_end_age > 90 ~ "Admin Censored",
      TRUE ~ str_to_title(main_dem_v1_end_type)
    ),
    
    # deriving flags for end type dementia, death, and end of membership
    
    event_endmem = ifelse(end_type == "End of membership", 1, 0),
    event_death = ifelse(end_type == "Death", 1, 0),
    event_dem = ifelse(end_type == "Dementia", 1, 0),
    
    # deriving 3-level education variable
    
    edu_3 = case_when(
      education_rev %in% c(1, 2) ~ "< HS",
      education_rev %in% c(3, 4) ~ "HS + Some college",
      education_rev %in% c(5, 6) ~ "College and above"
    ),
    edu_3 = factor(edu_3,
                   levels = c(
                     "< HS", "HS + Some college", "College and above"
                   ))
  )


# ----cleaning bp data----

length(unique(bp_data$subjid)) # n = 192,534

# filter for subjects that are in the analytic dataset

bp_data <- bp_data %>% filter(subjid %in% wide_data$subjid)

length(unique(bp_data$subjid)) # n = 157,607

# there are 1870 subjects with no bp data from EHR

nrow(wide_data) - length(unique(bp_data$subjid))

bp_data_clean <- bp_data %>%
  
  mutate(
    subjid = as.numeric(subjid),
    # merging categorical bp's with continuous
    continuous_dbp = ifelse(is.na(continuous_dbp), categorical_dbp, continuous_dbp),
    continuous_sbp = ifelse(is.na(continuous_sbp), categorical_sbp, continuous_sbp),
    # calculating age at event
    age = ifelse(
      as.Date(measure_date) == "1600-01-01",
      90,
      interval(as.Date("1900-01-01"), as.Date(measure_date)) / years(1)
    ),
    age_yr = round(age),
    flag_90plus = as.Date(measure_date) == "1600-01-01"
  ) %>%
  
  rename(dbp = continuous_dbp, sbp = continuous_sbp) %>%
  select(-categorical_dbp, -categorical_sbp, -measure_date, -enctype) %>%
  arrange(subjid, age)


# ----creating baseline age dataset----
age_data <-
  wide_data %>% mutate(survey_age_r = round(survey_age)) %>%
  select(subjid, survey_age_r)

# ----calculating median bp for each year----
bp_median <- bp_data_clean %>%
  filter(!flag_90plus) %>% # discarding obs that are 90+ censored
  group_by(subjid, age_yr) %>% # summarizing median bp's
  summarise(dbp_median = median(dbp),
            sbp_median = median(sbp)) %>%
  ungroup() %>%
  rename(dbp = dbp_median, sbp = sbp_median)

# ----merging in survey age rounded by yr----
bp_median <- age_data %>% left_join(bp_median, by = "subjid")

bp_med_vars <- bp_median %>%
  group_by(subjid) %>%
  filter(age_yr >= survey_age_r) %>%
  mutate(survey_age_sbp_diff = min(abs(survey_age_r - age_yr)),
         sbp_closest_to_baseline = if_else(survey_age_sbp_diff == min(survey_age_sbp_diff),
                                           sbp, NA_real_)) %>%
  slice(1) %>%
  select(subjid, survey_age_sbp_diff, sbp_closest_to_baseline)

wide_data_pre_mi <- left_join(wide_data, bp_med_vars, by = "subjid")
# ----saving pre-mi dataset----
saveRDS(wide_data_pre_mi,
        file = paste0(path_to_analytic_data, "wide_data_pre_mi.rds")) # 43 variables
