# setup and loading packages ----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load(
  "tidyverse", "here", "broom", 
  "survival", "splines"
  # "mice", "mitools", 
  # "purrr", "openxlsx",
  # "RColorBrewer",
  # "flextable"
)

options(scipen = 9999)

# set up paths and prep data ----

source(here("scripts", "0_paths.R"))
load(paste0(path_to_analytic_data, "imputed_tte_data_5.rds"))

# transform data
data_cut <- function(data) {
  cuts <- c(min(data$survey_age), # the minimum start age
            seq(60, 102, 1),     # all the intervals in the middle
            max(data$main_dem_v1_end_age))   # the maximum end of follow-up age
  
  longdata <- survSplit(
    data    = data,
    formula = Surv(survey_age, main_dem_v1_end_age, main_dem_v1_end_dem_flag) ~ .,
    cut     = cuts,
    # renaming the variables that indicate the start and end of each interval
    start   = "age_start",
    end     = "age_end"
  ) %>%
    mutate(
      current_age_c75 = age_start - 75,
      htn_age_int = htn_dx5yr_flag * current_age_c75
    )
  
  return(longdata)
}

longdata <- imputed_tte_data_5[[1]] %>% data_cut() 

longformula <- paste0(
  "Surv(age_start, age_end, main_dem_v1_end_dem_flag) ~ ",
  "htn_dx5yr_flag + htn_age_int + female + edu_3 + ",
  "usaborn_rev + ehr_ht_median + eversmoke + diab_dx5yr_flag"
)


fit_cox_ns <- function(d) {
  coxmod <- coxph(
    Surv(age_start, age_end, main_dem_v1_end_dem_flag) ~ 
      htn_dx5yr_flag + htn_dx5yr_flag:ns(current_age_c75, df = d) + 
      female + edu_3 + usaborn_rev + ehr_ht_median + eversmoke + diab_dx5yr_flag,
    subset = (ethnicity_rev == "White"), # subset to one ethnicity
    data = longdata
  )
  
  age_ns <- ns(seq(60, 90, by = 1) - 75, df = d) 
  cmat <- cbind(
    rep(1, nrow(age_ns)), 
    age_ns
  )
  
  index <- names(coxmod$coefficients) %>% str_detect("htn_dx5yr_flag")
  
  coef <- coef(coxmod)[index]
  vcov <- vcov(coxmod)[index, index]
  
  hr_age_df <- tibble(
    age = seq(60, 90, by = 1),
    est = cmat %*% coef %>% as.numeric(),
    se = diag(cmat %*% vcov %*% t(cmat)) %>% sqrt() %>% as.numeric()
  ) %>% 
    mutate(
      lower = est - 1.96 * se, 
      upper = est + 1.96 * se,
      across(c(est, lower, upper), exp),
      ns_df = d
    ) 
  
  return(hr_age_df)
}

# models with splines
test <- lapply(c(3:7), fit_cox_ns) %>% bind_rows() 

# fit the linear age-interaction model
coxmod <- coxph(
  Surv(age_start, age_end, main_dem_v1_end_dem_flag) ~ 
    htn_dx5yr_flag + htn_dx5yr_flag:current_age_c75 + 
    female + edu_3 + usaborn_rev + ehr_ht_median + eversmoke + diab_dx5yr_flag,
  subset = (ethnicity_rev == "White"), # subset to one ethnicity
  data = longdata
)

age_c75 <- seq(60, 90, by = 1) - 75
cmat <- cbind(rep(1, length(age_c75)), age_c75)
index <- names(coxmod$coefficients) %>% str_detect("htn_dx5yr_flag")
coef <- coef(coxmod)[index]
vcov <- vcov(coxmod)[index, index]

hr_age_df <- tibble(
  age = seq(60, 90, by = 1),
  est = cmat %*% coef %>% as.numeric(),
  se = diag(cmat %*% vcov %*% t(cmat)) %>% sqrt() %>% as.numeric()
) %>% 
  mutate(
    lower = est - 1.96 * se, 
    upper = est + 1.96 * se,
    across(c(est, lower, upper), exp)
  ) 

# plot together
test %>% 
  mutate(ns_df = paste0("df = ", ns_df)) %>% 
  ggplot(aes(x = age, y = est)) + 
  geom_line(aes(linetype = "splines")) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) + 
  geom_line(data = hr_age_df, aes(x = age, y = est, linetype = "linear")) +
  facet_wrap(~ ns_df, nrow = 1) +
  scale_x_continuous(breaks = seq(60, 90, 5), limits = c(65, 90)) + 
  scale_y_continuous(limits = c(0.4, 2.1), breaks = seq(0.5, 3, 0.25)) +
  scale_linetype_manual(values = c("splines" = 1, "linear" = 2)) +
  guides(linetype = "none") + 
  theme_bw() + 
  labs(x = "Age (years)", y = "Hazard ratio")

ggsave(
  width = 8, height = 3,
  file = here("scripts", "output", "figures", "supp_nonlinear_age_intxn.png")
)


