# Table 1 of complete case data
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/12/2024

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("here",
       "janitor",
       "tidyverse",
       "labelled",
       "magrittr",
       "table1",
       "openxlsx")

options(scipen = 9999)

# ----loading pre-mi data and functions----
source(here("scripts", "0_paths.R"))
source(here("scripts", "tbl1_functions.R"))

load(paste0(path_to_analytic_data, "cc_table1_data.RData"))

# ----adding membership criteria and variable labels to pre-imputation dataset----
wide_data_pre_mi_5 <- cc_table1_data %>%
  filter(presurv5yr_sample == 1) %>%
  tbl1_label()
wide_data_pre_mi_7 <- cc_table1_data %>%
  filter(presurv7yr_sample == 1) %>%
  tbl1_label()
wide_data_pre_mi_9 <- cc_table1_data %>%
  filter(presurv9yr_sample == 1) %>%
  tbl1_label()

# check dimensions again
dim(wide_data_pre_mi_5) # n = 136,180
dim(wide_data_pre_mi_7) # n = 122,032
dim(wide_data_pre_mi_9) # n = 114,832

# ----saving datasets for different exposure definitions----
saveRDS(wide_data_pre_mi_5,
        file = paste0(path_to_analytic_data, "wide_data_pre_mi_5.rds"))
saveRDS(wide_data_pre_mi_7,
        file = paste0(path_to_analytic_data, "wide_data_pre_mi_7.rds"))
saveRDS(wide_data_pre_mi_9,
        file = paste0(path_to_analytic_data, "wide_data_pre_mi_9.rds"))

cc_vars <-           "~ survey_age +
                        female +
                        usaborn_rev +
                        edu_3 +
                        maritalstatus +
                        smoking_status +
                        income_pp +
                        rx_1yr_pre_surv +
                        sbp_closest_to_baseline +
                        ehr_ht_median +
                        diab_dx5yr_flag +
                        first_prevstroke_flag +
                        main_dem_v1_end_type +
                        main_dem_v1_fu_time | ethnicity_rev + "


            
# ----creating table 1s----
tbl1_pre_mi_5 <- tbl1_render(
  wide_data_pre_mi_5,
  cc_vars,
  "htn_dx5yr_flag",
  "Hypertension status (assessed 5+ years before baseline)",
  my.render.cat_cc,
  paste0(
    "Summary statistics stratified on race/ethnicity and hypertension ",
    "diagnosis 5+yrs pre-survey using the pre-imputation dataset based ",
    "on 5+1 yr membership criteria."
  )
)

tbl1_pre_mi_7 <- tbl1_render(
  wide_data_pre_mi_7,
  cc_vars,
  "htn_dx7yr_flag",
  "Hypertension status (assessed 7+ years before baseline)",
  my.render.cat_cc,
  paste0(
    "Summary statistics stratified on race/ethnicity and hypertension ",
    "diagnosis 7+yrs pre-survey using the pre-imputation dataset based ",
    "on 7+1 yr membership criteria."
  )
)

tbl1_pre_mi_9 <- tbl1_render(
  wide_data_pre_mi_9,
  cc_vars,
  "htn_dx9yr_flag",
  "Hypertension status (assessed 9+ years before baseline)",
  my.render.cat_cc,
  paste0(
    "Summary statistics stratified on race/ethnicity and hypertension ",
    "diagnosis 9+yrs pre-survey using the pre-imputation dataset based ",
    "on 9+1 yr membership criteria."
  )
)

# ----ethnicity by hypertension header n's----
ethn_n <- cbind(
  with(wide_data_pre_mi_5, table(ethnicity_rev)),
  with(wide_data_pre_mi_7, table(ethnicity_rev)),
  with(wide_data_pre_mi_9, table(ethnicity_rev))
)

# ----ethnicity by hypertension header percentages----
ethn_x_htn <- cbind(
  with(wide_data_pre_mi_5, table(ethnicity_rev, htn_dx5yr_flag)) %>%
    proportions(margin = 1),
  with(wide_data_pre_mi_7, table(ethnicity_rev, htn_dx7yr_flag)) %>%
    proportions(margin = 1),
  with(wide_data_pre_mi_9, table(ethnicity_rev, htn_dx9yr_flag)) %>%
    proportions(margin = 1)
) %>% round_pad(digits = 3)

# ----list object with all table 1 info----
pre_mi_tbl1s <- list(
  "5+yr" = tbl1_pre_mi_5,
  "7+yr" = tbl1_pre_mi_7,
  "9+yr" = tbl1_pre_mi_9,
  "ethn_n" = ethn_n,
  "ethn_x_htn" = ethn_x_htn
)

# ----saving all table 1s----
write.xlsx(
  pre_mi_tbl1s,
  rowNames = TRUE,
  file = here("scripts", "output", "tables", "table_1_grouped.xlsx")
)
