# Cleaning Rx data and merging it into the pre-mi dataset for table 1 complete case
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/12/2024

if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("here", "janitor", "haven", "tidyverse")

options(scipen = 9999)

#---- loading data ----
source(here("scripts", "0_paths.R"))

pre_mi <- readRDS(paste0(path_to_analytic_data, "wide_data_pre_mi.rds"))
mi_data <- readRDS(paste0(path_to_analytic_data, "aa_htn_imputed_tbl1_data.rds"))

rx_data <- read_sas(
  paste0(
    path_to_box, 
    "Asian_Americans_dementia_data/diabetes_adrd/Rx_data/",
    "rx_data.sas7bdat"
  )
) %>% clean_names()

# derive htn medication variable ----
# htn medication use within 1-yr period prior to survey 
ages <- pre_mi %>% select(subjid, survey_age)

htn_rx_data <- rx_data %>%
  # keeping only htn treatments
  filter(str_detect(ucla_ahfs, paste(c("12:", "24:", "40:"), collapse = "|"))) %>% 
  right_join(., ages, by = "subjid")

# restrict 1 year prior to baseline
htn_rx_data <- htn_rx_data %>%
  mutate(
    rx_1yr_pre_surv = case_when(
      start_age > survey_age | end_age < survey_age - 1 ~ 0,
      is.na(start_age) | is.na(end_age) ~ 0,
      TRUE ~ 1
    )
  )

htn_rx_wide <- htn_rx_data %>% 
  select(subjid, rx_1yr_pre_surv) %>%
  arrange(subjid, rx_1yr_pre_surv) %>%
  group_by(subjid) %>%
  filter(row_number() == n())

# join variable to table 1 datasets and save ----

cc_table1_data <- left_join(pre_mi, htn_rx_wide, by = "subjid")
imputed_table1_data <- left_join(mi_data, htn_rx_wide, by = "subjid")

save(cc_table1_data, file = paste0(path_to_analytic_data, "cc_table1_data.RData"))
save(imputed_table1_data, file = paste0(path_to_analytic_data, "mi_table1_data.RData"))
