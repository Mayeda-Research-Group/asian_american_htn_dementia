
# some notes about runtime
# within each imputation, the Asian ethnicities generally take 1 hour in total 
# to finish running 1100 bootstraps, while for the White group, it takes 7.5 hrs


# Hoffman setup ----
# Read in the arguments listed in the:
# R CMD BATCH --no-save --no-restore "--args scenario_num=$SGE_TASK_ID"  
## expression:
args = (commandArgs(TRUE))

# Check to see if arguments are passed and set default values if not.
# If so, parse the arguments. (I only have one argument here.)
if (length(args) == 0) {
  print("No arguments supplied.")
  ##supply default values
  scenario_num <- 1
} else {
  for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
  }
}
# print values just to make sure:
print(scenario_num)

# also make sure that scenario_num is not greater than 20
if (scenario_num > 20) {
  print("out of bound senario number")
  # supply default value
  scenario_num <- 1
}


# setup and loading packages----

## locally ----

# if (!require("pacman"))
#   install.packages("pacman", repos = 'http://cran.us.r-project.org')
# 
# p_load(
#   "tidyverse", "here", 
#   "boot", "broom", "survival"
# )
# 
# options(scipen = 9999)

## on Hoffman ----
library(tidyverse)
library(boot)
library(broom)
library(survival)

# load and prep data----

# locally 

# source(here("scripts", "0_paths.R"))
# load(paste0(path_to_analytic_data, "imputed_tte_data_5.rds"))

## on Hoffman 
load("/u/home/y/yixzhou/AA_hypertension/data/imputed_tte_data_5.rds")

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


# bootstrap function ----

# fit a Cox model, extract point estimates of HRs (age-spec if so), 
# set up a hypothetical subject and calculate inc curves with/without htn,
# and calculate RR and RD


cox_boot <- function(data, formula, subj, indices) {
  # input data should be filtered to one imputation and a certain ethnicity
  # the input data should also be wide so that we re-sample by subject
  
  # # testing arguments
  # boot_data <- imputed_tte_data_long[[1]] %>% filter(ethnicity_rev == "Japanese")
  # # long format model formula
  # formula <- paste0(
  #   "Surv(age_start, age_end, main_dem_v1_end_dem_flag) ~ ",
  #   "htn_dx5yr_flag + htn_age_int + female + edu_3 + ",
  #   "usaborn_rev + ehr_ht_median + eversmoke + diab_dx5yr_flag"
  # )
  # subj <- test_subj_long
  
  boot_data <- data[indices,] %>% data_cut() # make long dataset
  
  # fit a Cox model
  coxmod <- coxph(
    as.formula(formula),
    # subset = (ethnicity_rev == ethn), # subset to one ethnicity
    data = boot_data
  )
  
  # collect the point estimates of the log(HR)s
  hr_pe <- coef(coxmod) # we can calculate age-spec HRs after bootstrapping
  
  # start prediction 
  pred_htn_0 <- survfit(
    coxmod, 
    newdata = subj %>% 
      mutate(htn_dx5yr_flag = 0, htn_age_int = htn_dx5yr_flag * current_age_c75),
    id = newid
  ) %>% 
    tidy() %>% 
    mutate(htn = 0, time = time + 60)
  
  pred_htn_1 <- survfit(
    coxmod, 
    newdata = subj %>% 
      mutate(htn_dx5yr_flag = 1, htn_age_int = htn_dx5yr_flag * current_age_c75),
    id = newid
  ) %>% 
    tidy() %>% 
    mutate(htn = 1, time = time + 60)
  
  # clean up the predicted incidence curves
  pred_inc <- bind_rows(pred_htn_0, pred_htn_1) %>% 
    mutate(time_cut = cut(time, seq(59, 90, 1))) %>%
    group_by(htn, time_cut) %>%
    slice_tail() %>%
    ungroup() %>%
    select(time, estimate, htn) %>% 
    mutate(
      estimate = 1 - estimate,
      colname = paste0("time_", time, "_htn_", htn)
    ) 
  
  # transform into a vector for output
  inc <- pred_inc$estimate
  names(inc) <- pred_inc$colname
  
  # calculate RR and RD at 5-year intervals
  RR_RD_5 <- pred_inc %>%
    pivot_wider(
      id_cols = "time", 
      names_from = "htn",
      names_prefix = "htn_", 
      values_from = "estimate"
    ) %>% 
    mutate(
      RD = htn_1 - htn_0, 
      RR = htn_1 / htn_0,
      .keep = "unused"
    ) %>% 
    filter(time %in% seq(60, 90, 5)) 
  
  # transform into a vector for output
  RR <- RR_RD_5$RR
  RD <- RR_RD_5$RD
  names(RR) <- paste0("RR_", RR_RD_5$time)
  names(RD) <- paste0("RD_", RR_RD_5$time)
  
  out <- c(hr_pe, inc, RR, RD)
  
  return(out)
}

# set up subject to make predictions on ----

test_subj <- tibble(
  female = factor(1), 
  edu_3 = factor(
    "HS + Some college", 
    levels = c("< HS", "HS + Some college", "College and above")
  ), 
  usaborn_rev = 1, 
  ehr_ht_median = 66,
  eversmoke = factor("No", levels = c("Yes", "No")), 
  diab_dx5yr_flag = 0
)

test_subj_long <- test_subj %>% 
  mutate(newid = 1) %>%
  right_join(expand_grid(age_start = 60:89, newid = 1), by = "newid") %>%
  mutate(
    current_age_c75 = age_start - 75,
    age_end = age_start + 1,
    main_dem_v1_end_dem_flag = 0
  )


# bootstrap starts here ----

ethnicities <- levels(imputed_tte_data_5[[1]]$ethnicity_rev)

# set up the index of imputations using scenario_num
# scenario_num <- 2
# imp_start <- (scenario_num - 1) * 2 + 1
# imp_end <- scenario_num * 2
imp_start <- scenario_num
imp_end <- scenario_num
seed <- scenario_num * 1234
set.seed(seed)

# testing locally
# set.seed(1234)
# imp_start <- 1
# imp_end <- 1

longformula <- paste0(
  "Surv(age_start, age_end, main_dem_v1_end_dem_flag) ~ ",
  "htn_dx5yr_flag + htn_age_int + female + edu_3 + ",
  "usaborn_rev + ehr_ht_median + eversmoke + diab_dx5yr_flag"
)

for (i in imp_start:imp_end) {
  # i <- 1
  print(i)
  
  for (ethn in ethnicities) {
    # ethn <- ethnicities[2]
    
    print(ethn)
    start <- Sys.time() # record starting time
    
    # subset to the specific ethnicity 
    input_data <- imputed_tte_data_5[[i]] %>% filter(ethnicity_rev == ethn) 
    
    # run bootstrapping
    res <- boot(
      data = input_data,
      statistic = cox_boot,
      R = 1100,
      formula = longformula, 
      subj = test_subj_long
    )
    
    end <- Sys.time() # record ending time
    
    # save bootstrapping results to output
    # res$t0 is results using full dataset, no resampling
    # res$t is results with resampling
    out_df <- rbind(res$t0, res$t)
    # append with sampling number and imputation number (i)
    out_df <- cbind(out_df, t = c(0:res$R), imp = i, ethnicity = ethn) %>% as_tibble()
    
    write_csv(
      out_df, 
      file = paste0(
        # Hoffman path
        "/u/home/y/yixzhou/AA_hypertension/results/",
        "bootstrap_imp_", i, "_", ethn, ".csv"
      )
    )
    
    # save runtime
    time <- tibble(imp = i, ethnicity = ethn, runtime = end - start)
    
    write_csv(
      time,
      # create new when it's the first imp and the first ethn; append when not 
      append = !(i == 1 & ethn == ethnicities[[1]]), 
      file = paste0(
        # Hoffman path
        "/u/home/y/yixzhou/AA_hypertension/results/",
        "bootstrap_runtime.csv"
      )
    )
    
  }
}

