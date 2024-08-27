# Table 1 of complete case data
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/12/2024

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("here", "haven", "janitor", "tidyverse", 
       "labelled", "table1", "openxlsx", "flextable")

options(scipen = 9999)

# ----loading pre-mi data and functions----
source(here("scripts", "0_paths.R"))
source(here("scripts", "tbl1_functions.R"))

load(paste0(path_to_analytic_data, "cc_table1_data.RData"))

# ----adding membership criteria and variable labels to pre-imputation dataset----
wide_data_pre_mi_5 <- cc_table1_data %>%
  filter(presurv5yr_sample == 1) %>%
  tbl1_label()

# wide_data_pre_mi_7 <- cc_table1_data %>%
#   filter(presurv7yr_sample == 1) %>%
#   tbl1_label()
# wide_data_pre_mi_9 <- cc_table1_data %>%
#   filter(presurv9yr_sample == 1) %>%
#   tbl1_label()

# check dimensions again
dim(wide_data_pre_mi_5) # n = 136,180

# dim(wide_data_pre_mi_7) # n = 122,032
# dim(wide_data_pre_mi_9) # n = 114,832

# ----saving datasets for different exposure definitions----
# saveRDS(wide_data_pre_mi_5,
#         file = paste0(path_to_analytic_data, "wide_data_pre_mi_5.rds"))
# saveRDS(wide_data_pre_mi_7,
#         file = paste0(path_to_analytic_data, "wide_data_pre_mi_7.rds"))
# saveRDS(wide_data_pre_mi_9,
#         file = paste0(path_to_analytic_data, "wide_data_pre_mi_9.rds"))

cc_vars <-           "~ survey_age +
                        female +
                        usaborn_rev +
                        edu_3 +
                        maritalstatus +
                        income_pp +
                        generalhealth_3 +
                        smoking_status + eversmoke + 
                        generalhealth_3 + 
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

# tbl1_pre_mi_7 <- tbl1_render(
#   wide_data_pre_mi_7,
#   cc_vars,
#   "htn_dx7yr_flag",
#   "Hypertension status (assessed 7+ years before baseline)",
#   my.render.cat_cc,
#   paste0(
#     "Summary statistics stratified on race/ethnicity and hypertension ",
#     "diagnosis 7+yrs pre-survey using the pre-imputation dataset based ",
#     "on 7+1 yr membership criteria."
#   )
# )
# 
# tbl1_pre_mi_9 <- tbl1_render(
#   wide_data_pre_mi_9,
#   cc_vars,
#   "htn_dx9yr_flag",
#   "Hypertension status (assessed 9+ years before baseline)",
#   my.render.cat_cc,
#   paste0(
#     "Summary statistics stratified on race/ethnicity and hypertension ",
#     "diagnosis 9+yrs pre-survey using the pre-imputation dataset based ",
#     "on 9+1 yr membership criteria."
#   )
# )

# ethnicity by hypertension header perc's----

ethn_htn_perc <- wide_data_pre_mi_5 %>% 
  count(ethnicity_rev, htn_dx5yr_flag) %>% 
  group_by(ethnicity_rev) %>% 
  mutate(ethn_sum = sum(n)) %>% 
  ungroup() %>% 
  mutate(
    perc = n / ethn_sum * 100,
    perc = format(round(perc, 1), nsmall = 1) %>% paste0("%"),
    n = paste0("N=", n)
  ) %>% 
  select(-ethn_sum)

# ----saving all table 1s----

# save_as_docx(
#   tbl1_pre_mi_5 %>% t1flex() %>% autofit(),
#   ethn_htn_perc %>% flextable(),
#   path = here("scripts", "output", "tables", "table_1_pre_MI.docx")
# )

write.xlsx(
  list(
    "pre-MI table 1" = tbl1_pre_mi_5,
    "header percentages" = ethn_htn_perc
  ),
  rowNames = TRUE,
  file = here("scripts", "output", "tables", "table_1_pre_MI.xlsx")
)

# look at additional descriptives ----

## survey age distribution by ethn and htn ----

# histogram is harder to look at
# facet_grid: can't allow free_y scale
# facet_wrap: free_y scale looks messy and hard to compare
wide_data_pre_mi_5 %>%
  mutate(htn_dx5yr_flag = ifelse(htn_dx5yr_flag == "Yes", "Hypertension", "No hypertension")) %>% 
  ggplot(aes(x = survey_age, fill = htn_dx5yr_flag, group = ethnicity_rev)) + 
  geom_histogram(alpha = 0.7, show.legend = FALSE) +
  facet_grid(row = vars(htn_dx5yr_flag), col = vars(ethnicity_rev)) + # 1
  # facet_wrap(vars(htn_dx5yr_flag, ethnicity_rev), nrow = 2, scales = "free_y") + # 2
  theme_bw() + 
  labs(x = "Survey age (years)", y = "Count")

ggsave(
  here("scripts", "output", "figures", "supp_surveyage_dist_hist_1.png"),
  height = 4, width = 8
)

# density
wide_data_pre_mi_5 %>%
  mutate(htn_dx5yr_flag = ifelse(htn_dx5yr_flag == "Yes", "Hypertension", "No hypertension")) %>% 
  ggplot(aes(x = survey_age, color = htn_dx5yr_flag, group = ethnicity_rev)) + 
  geom_density(show.legend = FALSE) +
  facet_grid(row = vars(htn_dx5yr_flag), col = vars(ethnicity_rev)) +
  # facet_wrap(vars(htn_dx5yr_flag, ethnicity_rev), nrow = 2, scales = "free_y") + 
  theme_bw() + 
  labs(x = "Survey age (years)", y = "Density")

ggsave(
  here("scripts", "output", "figures", "supp_surveyage_dist_dens.png"),
  height = 4, width = 8
)
 
## htn at baseline and anti-htn use ----

aa_raw_data <- read_sas(
  paste0(path_to_raw_data, "aa_adrd_cardiometabolic_tte.sas7bdat")
) %>% clean_names()

wide_data_pre_mi_5 <- aa_raw_data %>% 
  select(subjid, first_htn_dx_age) %>% 
  right_join(wide_data_pre_mi_5, by = "subjid")

summary(wide_data_pre_mi_5$first_htn_dx_age)

wide_data_pre_mi_5 <- wide_data_pre_mi_5 %>% 
  mutate(
    htn_dx_by_survey = case_when(
      first_htn_dx_age <= survey_age ~ "Yes",
      TRUE ~ "No"
      # some subjects have first htn dx after baseline, 
      # some don't have htn dx at all
    ) %>% 
      factor(levels = c("Yes", "No"))
  ) 

wide_data_pre_mi_5 %>% 
  filter(htn_dx5yr_flag == "No") %>% 
  table1(
    # ~ rx_1yr_pre_surv | ethnicity_rev + htn_dx_by_survey, 
    ~ htn_dx_by_survey | ethnicity_rev + rx_1yr_pre_surv, 
    data = .
  )

