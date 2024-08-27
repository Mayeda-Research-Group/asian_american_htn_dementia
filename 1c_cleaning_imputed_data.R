# Multiple imputation of wide dataset
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/12/2024

# ----setup and loading packages----

if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("here", "tidyverse", "mice", "mitools")

options(scipen = 9999)

# ----loading data----
source(here("scripts", "0_paths.R"))

pre_mi <- readRDS(paste0(path_to_analytic_data, "wide_data_pre_mi.rds"))

load(paste0(path_to_analytic_data, "imputed_raw.RData"))

# check imputation diagnostics ----
imputed$loggedEvents #NULL

plot(imputed)

densityplot(imputed, ~ income)
densityplot(imputed, ~ sizeofhh)
densityplot(imputed, ~ income_pp)
densityplot(imputed, ~ education_rev)
densityplot(imputed, ~ usaborn_rev)
densityplot(imputed, ~ generalhealth)


impute_var_list <- c(
  "survey_age",
  "female",
  "ethnicity_rev",
  "usaborn_rev",
  "usabornfather_rev",
  "usabornmother_rev",
  "education_rev",
  "sizeofhh",
  "income",
  "income_pp",
  "maritalstatus",
  "generalhealth", 
  "smoking_status",
  "sr_total_height_in",
  "ehr_ht_median",
  "sbp_closest_to_baseline"
)

other_vars <- pre_mi %>% select(-all_of(impute_var_list))

# create stacked dataset for table 1 ----

imp_temp <- lapply(
  1:imputed$m, function(i) {
    complete(imputed, action = i) %>% 
      cbind("subjid" = pre_mi$subjid) %>% 
      mutate(imp = i) %>% 
      relocate(subjid, imp, .before = everything())
  }
)

imp_data_stacked <- do.call(rbind, imp_temp)
colnames(imp_data_stacked)

tbl1_imputed_data <- imp_data_stacked %>%
  left_join(other_vars, by = "subjid") %>%
  filter(presurv5yr_sample == 1) %>% # can change this out for 7 or 9 year membership criteria
  mutate(
    edu_ge_college = ifelse(education_rev %in% c(1, 2, 3, 4), 0, 1),
    edu_3 = case_when(
      education_rev %in% c(1, 2) ~ "< HS",
      education_rev %in% c(3, 4) ~ "HS + Some college",
      education_rev %in% c(5, 6) ~ "College and above"
    ) %>% factor(levels = c("< HS", "HS + Some college", "College and above")),
    married = case_when(
      maritalstatus %in% c(2) ~ "Married",
      maritalstatus %in% c(1,3,4) ~ "Not married"
    ) %>% factor(levels = c("Married", "Not married")),
    eversmoke = case_when(
      smoking_status == 1 ~ "No", 
      smoking_status %in% c(2, 3) ~ "Yes"
    ) %>% factor(levels = c("Yes", "No"))
  )

saveRDS(
  tbl1_imputed_data,
  file = paste0(path_to_analytic_data, "aa_htn_imputed_tbl1_data.rds")
)


# create mild object for analysis ----
aa_htn_wide_mild <- complete(imputed, action = 'all')

# adding variables back into the mild object
for (i in 1:imputed$m) {
  aa_htn_wide_mild[[i]] <- aa_htn_wide_mild[[i]] %>%
    mutate(subjid = pre_mi$subjid) %>% 
    left_join(other_vars, by = "subjid") %>%
    mutate(
      ethnicity_rev = factor(
        ethnicity_rev, levels = c(2, 5, 3, 1, 9),
        labels = c("Chinese", "Filipino", "Japanese", "South Asian", "White")
      ),
      edu_ge_college = ifelse(education_rev %in% c(1, 2, 3, 4), 0, 1),
      edu_3 = case_when(
        education_rev %in% c(1, 2) ~ "< HS",
        education_rev %in% c(3, 4) ~ "HS + Some college",
        education_rev %in% c(5, 6) ~ "College and above"
      ) %>% factor(levels = c("< HS", "HS + Some college", "College and above")),
      eversmoke = case_when(
        smoking_status == 1 ~ "No", 
        smoking_status %in% c(2, 3) ~ "Yes"
      ) %>% factor(levels = c("Yes", "No"))
    ) %>% 
    relocate(subjid, .before = everything())
}

# ----applying different membership criteria and saving imputed datasets----
imputed_tte_data_5 <- lapply(aa_htn_wide_mild, function(df) {
  filter(df, presurv5yr_sample == 1)
})

save(imputed_tte_data_5,
     file = paste0(path_to_analytic_data, "imputed_tte_data_5.rds"))


# did not end up using these datasets 

# imputed_tte_data_7 <- lapply(aa_htn_wide_mild, function(df) {
#   filter(df, presurv7yr_sample == 1)
# })
# 
# imputed_tte_data_9 <- lapply(aa_htn_wide_mild, function(df) {
#   filter(df, presurv9yr_sample == 1)
# })

# save(imputed_tte_data_7,
#      file = paste0(path_to_analytic_data, "imputed_tte_data_7.rds"))
# save(imputed_tte_data_9,
#      file = paste0(path_to_analytic_data, "imputed_tte_data_9.rds"))
