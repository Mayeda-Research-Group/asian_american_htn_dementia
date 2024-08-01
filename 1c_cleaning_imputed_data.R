# Cleaning imputed dataset
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

pre_mi <-
  readRDS(paste0(path_to_analytic_data, "wide_data_pre_mi.rds"))

load(paste0(path_to_analytic_data, "imputed_raw.RData"))

# ----looking at mi object to see if anything went wrong----
imputed$loggedEvents #NULL

plot(imputed)

densityplot(imputed, ~ income)
densityplot(imputed, ~ sizeofhh)
densityplot(imputed, ~ income_pp)
densityplot(imputed, ~ education_rev)
densityplot(imputed, ~ usaborn_rev)


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
  "smoking_status",
  "sr_total_height_in",
  "ehr_ht_median",
  "sbp_closest_to_baseline"
)

imp_temp <- lapply(1:imputed$m, function(i) {
  temp <- complete(imputed, action = i)
  temp <- cbind(pre_mi$subjid, temp[, c(impute_var_list)])
  temp$imp <- i
  return(temp)
})

imp_data_stacked <- do.call(rbind, imp_temp)
colnames(imp_data_stacked)
names(imp_data_stacked)[names(imp_data_stacked) == "pre_mi$subjid"] <-
  "subjid"

# ----creating mild object for analysis----
aa_htn_wide_mild <- complete(imputed, action = 'all')


# adding variables back in
other_vars <- pre_mi %>% select(-all_of(impute_var_list))

# ----merge into stacked imputed dataset for table 1----
tbl1_imputed_data <- imp_data_stacked %>%
  right_join(other_vars, by = "subjid") %>%
  filter(presurv5yr_sample == 1) %>% # can change this out for 7 or 9 year membership criteria
  mutate(
    edu_ge_college = ifelse(education_rev %in% c(1, 2, 3, 4), 0, 1),
    edu_3 = case_when(
      education_rev %in% c(1, 2) ~ "< HS",
      education_rev %in% c(3, 4) ~ "HS + Some college",
      education_rev %in% c(5, 6) ~ "College and above"
    ),
    edu_3 = factor(edu_3,
                   levels = c(
                     "< HS", "HS + Some college", "College and above"
                   )),
    married = as.factor(ifelse(maritalstatus == 2, "Married",
                               ifelse(maritalstatus %in% c(1,3,4), "Not married", NA)))
  )

saveRDS(tbl1_imputed_data,
        file = paste0(path_to_analytic_data,
                      "aa_htn_imputed_tbl1_data.rds"))


# merging into mild object for analysis
for (i in 1:imputed$m) {
  aa_htn_wide_mild[[i]] <- aa_htn_wide_mild[[i]] %>%
    cbind("subjid" = pre_mi$subjid) %>%
    merge(x = ., y = other_vars, by = "subjid") %>%
    # filter(presurv5yr_sample == 1) %>% # can change this out for 7 or 9 year membership criteria
    mutate(
      edu_ge_college = ifelse(education_rev %in% c(1, 2, 3, 4), 0, 1),
      edu_3 = case_when(
        .$education_rev %in% c(1, 2) ~ "< HS",
        .$education_rev %in% c(3, 4) ~ "HS + Some college",
        .$education_rev %in% c(5, 6) ~ "College and above"
      ),
      edu_3 = factor(edu_3,
                     levels = c(
                       "< HS", "HS + Some college", "College and above"
                     ))
    )
}

# ----applying different membership criteria and saving imputed datasets----
imputed_tte_data_5 <- lapply(aa_htn_wide_mild, function(df) {
  filter(df, presurv5yr_sample == 1)
})

imputed_tte_data_7 <- lapply(aa_htn_wide_mild, function(df) {
  filter(df, presurv7yr_sample == 1)
})

imputed_tte_data_9 <- lapply(aa_htn_wide_mild, function(df) {
  filter(df, presurv9yr_sample == 1)
})

save(imputed_tte_data_5,
     file = paste0(path_to_analytic_data, "imputed_tte_data_5.rds"))
save(imputed_tte_data_7,
     file = paste0(path_to_analytic_data, "imputed_tte_data_7.rds"))
save(imputed_tte_data_9,
     file = paste0(path_to_analytic_data, "imputed_tte_data_9.rds"))
