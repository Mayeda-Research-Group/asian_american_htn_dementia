
# load packages ----

if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("tidyverse", "here")

options(scipen = 9999)

# read results ----

boot_files <- list.files(
  here("scripts", "output", "bootstrap", "results"),
  pattern = "bootstrap_imp", full.names = TRUE
)

bootstraps <- lapply(boot_files, read_csv) %>% bind_rows()

# brief check
bootstraps %>% View()
bootstraps %>% count(ethnicity, imp) %>% View()

# calculate age-spec HR associated hypertension ----

bootstraps <- bootstraps %>% 
  mutate(
    htn_65 = htn_dx5yr_flag + htn_age_int * (65 - 75), 
    htn_70 = htn_dx5yr_flag + htn_age_int * (70 - 75), 
    htn_75 = htn_dx5yr_flag + htn_age_int * (75 - 75), 
    htn_80 = htn_dx5yr_flag + htn_age_int * (80 - 75), 
    htn_85 = htn_dx5yr_flag + htn_age_int * (85 - 75), 
    htn_90 = htn_dx5yr_flag + htn_age_int * (90 - 75)
  )

# prepare point estimates and confidence intervals ----

pe <- bootstraps %>%
  filter(t == 0) %>% 
  group_by(ethnicity) %>% 
  summarise(across(-c(t, imp), function(x) mean(x), .names = "{.col}.pe")) 


# RR is missing for certain samples since inc in the no htn arm is 0
# I think this is inherent in resampling
# I kind of want to treat these missings as extreme values at either end
# and include them when taking the percentiles
# since otherwise the CI would be narrower

# all the missingness is coming from RR at 65 and 70
bootstraps %>% filter(t != 0) %>% 
  group_by(ethnicity) %>% 
  summarise(
    nan_65 = sum(is.nan(RR_65)), 
    nan_70 = sum(is.nan(RR_70)),
    nan_75 = sum(is.nan(RR_75))
  )

# when the RR's are missing (NaN due to dividing by 0),
# I randomly assign -999 or 999 so that they lie on the extreme ends of the 
# bootstrap distribution
bootstraps_edit <- bootstraps %>% 
  filter(t != 0) %>%
  rowwise() %>% 
  mutate(
    RR_65 = ifelse(is.nan(RR_65), sample(c(-999, 999), size = 1, replace = TRUE), RR_65),
    RR_70 = ifelse(is.nan(RR_70), sample(c(-999, 999), size = 1, replace = TRUE), RR_70)
  ) 

# check missingness again to make sure
bootstraps_edit %>% na.omit() %>% nrow()

confint <- bootstraps_edit %>%
  filter(t <= 1000) %>% # take 1000 bootstraps for each ethn and each imp
  group_by(ethnicity) %>% 
  summarise(
    across(
      -c(t, imp), 
      list("lower" = function(x) quantile(x, 0.025), 
           "upper" = function(x) quantile(x, 0.975)),
      .names = "{.col}.{.fn}"
    )
  ) 

confint %>% 
  select(ethnicity, starts_with("RR"))
# it seems that RR at 65 has unprecise CI bounds (both ends) 
# for Chinese, Japanese, and South Asian,
# and RR at 70 has unprecise lower bound for South Asian. 

# might not want to include in the main results
# compare with the restricted bootstraps -- seems overly precise
# we can probably match with the age-spec incidence rates table
# in terms of what estimates to show and what not to
bootstraps %>% 
  group_by(ethnicity) %>% 
  summarise(
    across(
      -c(t, imp), 
      list("lower" = function(x) quantile(x, 0.025, na.rm = TRUE), 
           "upper" = function(x) quantile(x, 0.975, na.rm = TRUE)),
      .names = "{.col}.{.fn}"
    )
  ) %>% 
  select(ethnicity, starts_with("RR"))

# join and save the point estimates and the CIs
allest <- pe %>% left_join(confint, by = "ethnicity")

save(allest, file = here("scripts", "output", "bootstrap_results.RData"))
