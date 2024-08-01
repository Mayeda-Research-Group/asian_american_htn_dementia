# Table 1 of imputed data
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/12/2024

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("here",
       "tidyverse",
       "labelled",
       "magrittr",
       "table1",
       "openxlsx")

options(scipen = 9999)

# ----loading data and functions----
source(here("scripts", "0_paths.R"))
source(here("scripts", "tbl1_functions.R"))

load(paste0(path_to_analytic_data, "mi_table1_data.RData"))

# ----adding membership criteria and variable labels to post-imputation dataset----
wide_data_imp_5 <- imputed_table1_data %>%
  filter(presurv5yr_sample == 1) %>%
  tbl1_label()
wide_data_imp_7 <- imputed_table1_data %>%
  filter(presurv7yr_sample == 1) %>%
  tbl1_label()
wide_data_imp_9 <- imputed_table1_data %>%
  filter(presurv9yr_sample == 1) %>%
  tbl1_label()

# ----saving datasets for different exposure definitions----
saveRDS(wide_data_imp_5,
        file = paste0(path_to_analytic_data, "wide_data_imp_5.rds"))
saveRDS(wide_data_imp_7,
        file = paste0(path_to_analytic_data, "wide_data_imp_7.rds"))
saveRDS(wide_data_imp_9,
        file = paste0(path_to_analytic_data, "wide_data_imp_9.rds"))

# check dimensions
nrow(wide_data_imp_5) / 20 # n = 136,180

imp_vars <-           "~ female +
                        usaborn_rev +
                        edu_3 +
                        maritalstatus +
                        smoking_status +
                        rx_1yr_pre_surv +
                        diab_dx5yr_flag +
                        first_prevstroke_flag +
                        main_dem_v1_end_type | ethnicity_rev + " # rx is same as complete case

# ----creating table 1s: categorical vars----
tbl1_imp_5 <- tbl1_render(
  wide_data_imp_5,
  imp_vars,
  "htn_dx5yr_flag",
  "Hypertension status (assessed 5+ years before baseline)",
  my.render.cat_imp,
  paste0(
    "Summary statistics stratified on race/ethnicity and hypertension ",
    "diagnosis 5+yrs pre-survey using the pre-imputation dataset based ",
    "on 5+1 yr membership criteria."
  )
)

tbl1_imp_7 <- tbl1_render(
  wide_data_imp_7,
  imp_vars,
  "htn_dx7yr_flag",
  "Hypertension status (assessed 7+ years before baseline)",
  my.render.cat_imp,
  paste0(
    "Summary statistics stratified on race/ethnicity and hypertension ",
    "diagnosis 7+yrs pre-survey using the pre-imputation dataset based ",
    "on 7+1 yr membership criteria."
  )
)

tbl1_imp_9 <- tbl1_render(
  wide_data_imp_9,
  imp_vars,
  "htn_dx9yr_flag",
  "Hypertension status (assessed 9+ years before baseline)",
  my.render.cat_imp,
  paste0(
    "Summary statistics stratified on race/ethnicity and hypertension ",
    "diagnosis 9+yrs pre-survey using the pre-imputation dataset based ",
    "on 9+1 yr membership criteria."
  )
)

cat_post_mi_tbl1s <- list("5+yr" = tbl1_imp_5,
                          "7+yr" = tbl1_imp_7,
                          "9+yr" = tbl1_imp_9)

write.xlsx(
  cat_post_mi_tbl1s,
  rowNames = TRUE,
  file = here(
    "scripts",
    "output",
    "tables",
    "cat_post_mi_tbl1s_subgrp.xlsx"
  )
)

# ----using function to apply Rubin's rules to cont vars----
# ehr_ht_median, income_pp, and sbp_closest_to_baseline are imputed

cont_imp_vars <-
  list(
    "income_pp" = 0,
    "ehr_ht_median" = 1,
    "sbp_closest_to_baseline" = 1
  )

tbl1_cont_post_mi_5 <-
  imp_cont_var_tbl1(wide_data_imp_5,
                    cont_imp_vars,
                    c("ethnicity_rev", "htn_dx5yr_flag"))

tbl1_cont_post_mi_7 <-
  imp_cont_var_tbl1(wide_data_imp_7,
                    cont_imp_vars,
                    c("ethnicity_rev", "htn_dx7yr_flag"))
tbl1_cont_post_mi_9 <-
  imp_cont_var_tbl1(wide_data_imp_9,
                    cont_imp_vars,
                    c("ethnicity_rev", "htn_dx9yr_flag"))

cont_post_mi_tbl1s <- list("5+yr" = tbl1_cont_post_mi_5,
                           "7+yr" = tbl1_cont_post_mi_7,
                           "9+yr" = tbl1_cont_post_mi_9)

write.xlsx(
  cont_post_mi_tbl1s,
  rowNames = TRUE,
  file = here(
    "scripts",
    "output",
    "tables",
    "cont_post_mi_tbl1s_subgrp.xlsx"
  )
)
