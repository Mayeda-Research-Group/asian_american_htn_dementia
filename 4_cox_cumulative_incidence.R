# Purpose : Using cumulative incidence to calculate RR and RD
# Created by: Juliet Zhou
# Creation date: 7/16/2024

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load(
  "tidyverse",
  "here",
  "broom",
  "survival",
  "mice", 
  "mitools", 
  "msm" # for deltamethod()
)

options(scipen = 9999)

# ----loading data and function----

source(here("scripts", "0_paths.R"))
load(paste0(path_to_analytic_data, "imputed_tte_data_5.rds"))

# prepare data ----

data_cut <- function(data) {
  cuts <- c(min(data$survey_age), # the minimum start age
            seq(60, 102, 1),     # all the intervals in the middle
            max(data$end_age))   # the maximum end of follow-up age
  
  longdata <- survSplit(
    data    = data,
    formula = Surv(survey_age, end_age, event_dem) ~ .,
    cut     = cuts,
    # renaming the variables that indicate the start and end of each interval
    start   = "age_start",
    end     = "age_end"
  ) %>%
    mutate(current_age_c75 = age_start - 75,
           htn_age_int = htn_dx5yr_flag * current_age_c75)
  
  return(longdata)
}

imputed_tte_data_long <- list()
for (i in seq_along(imputed_tte_data_5)) {
  imputed_tte_data_long[[i]] <- data_cut(imputed_tte_data_5[[i]])
}

# set up test subjects to predict on ----

# prection function ----
predict_cumulative_inc <- function(data, ethn, test_subj, coxformula, longdata) {
  
  # # testing arguments
  # data <- imputed_tte_data_5
  # longdata <- FALSE
  # coxformula <- "Surv(survey_age, end_age, event_dem) ~ htn_dx5yr_flag + female + edu_3 + usaborn_rev + ehr_ht_median"
  # test_subj <- imputed_tte_data_5[[1]] %>% slice_head()
  # # or
  # data <- imputed_tte_data_long
  # longdata <- TRUE
  # coxformula <- "Surv(age_start, age_end, event_dem) ~ htn_dx5yr_flag + htn_age_int + female + edu_3 + usaborn_rev + ehr_ht_median"
  # test_subj <- imputed_tte_data_5[[1]] %>% 
  #   slice_head() %>% 
  #   select(female, edu_3, usaborn_rev, ehr_ht_median) %>% 
  #   mutate(newid = 1) %>% 
  #   right_join(expand_grid(age_start = 60:89, newid = 1), by = "newid") %>% 
  #   mutate(
  #     current_age_c75 = age_start - 75,
  #     age_end = age_start + 1, 
  #     event_dem = 0
  #   )
  
  # ethn <- 3
  
  m <- length(data)
  sep_models <- vector(mode = "list", length = m)
  sep_pred <- vector(mode = "list", length = m)
  
  
  for (i in 1:m) {
    # i <- 1
    
    data_i <- data[[i]]
    
    # those with event_dem = 1 developed dementia at end_age
    # those with event_dem = 0 ended followup for non-dementia reasons
    sep_models[[i]] <-
      coxph(
        as.formula(coxformula),
        subset = (ethnicity_rev == ethn), # subset to one ethnicity
        data = data_i
      )
    
    # survfit.coxph generates survival/cumulative incidence based on the model fit
    # similar to a predict function
    # this prediction is performed for the same subject 
    # once assuming they do not have htn, and once assuming they do
    
    # code varies slightly based on whether we use a long dataset 
    # to accommodate age-varying htn effect
    if (longdata) {
      
      pred_htn_0 <- survfit(
        sep_models[[i]], 
        newdata = test_subj %>% 
          mutate(htn_dx5yr_flag = 0, htn_age_int = htn_dx5yr_flag * current_age_c75),
        se.fit = TRUE,
        id = newid
      ) %>% 
        tidy() %>% 
        mutate(htn = 0, time = time + 60)
      
      pred_htn_1 <- survfit(
        sep_models[[i]], 
        newdata = test_subj %>% 
          mutate(htn_dx5yr_flag = 1, htn_age_int = htn_dx5yr_flag * current_age_c75),
        se.fit = TRUE,
        id = newid
      ) %>% 
        tidy() %>% 
        mutate(htn = 1, time = time + 60)
      
    } else {
      
      pred_htn_0 <- survfit(
        sep_models[[i]], 
        newdata = test_subj %>% mutate(htn_dx5yr_flag = 0),
        se.fit = TRUE
      ) %>% 
        tidy() %>% 
        mutate(htn = 0)
      
      pred_htn_1 <- survfit(
        sep_models[[i]], 
        newdata = test_subj %>% mutate(htn_dx5yr_flag = 1),
        se.fit = TRUE
      ) %>% 
        tidy() %>% 
        mutate(htn = 1)
    }
    
    # collect this result into a list 
    sep_pred[[i]] <- bind_rows(pred_htn_0, pred_htn_1) %>% 
      # slice at select time points
      mutate(time_cut = cut(time, seq(59, 90, 1))) %>%
      group_by(htn, time_cut) %>%
      slice_tail() %>%
      ungroup() %>%
      select(-time_cut) 
    
  }
  
  # use MIcombine to pool these predictions
  ests <- lapply(sep_pred, function (x) x$estimate)
  vars <- lapply(sep_pred, function (x) diag(x$std.error ^ 2)) # need to be a matrix
  MIcombine_res <- MIcombine(ests, vars) 
  
  RD_RR_results <- bind_cols(
    # add in time and htn status 
    sep_pred[[1]] %>% select(time, htn),
    est = 1 - MIcombine_res$coefficients, # want incidence instead of survival
    var = MIcombine_res$variance %>% diag() # convert matrix back to vector
  ) %>% 
    # # slice at time points of interest
    # mutate(time_cut = cut(time, seq(60, 90, 5))) %>%
    # group_by(htn, time_cut) %>%
    # slice_tail() %>%
    # ungroup() %>%
    # widen the dataset by htn status
    pivot_wider(
      id_cols = time,
      names_from = htn,
      values_from = c(est, var), 
      names_glue = "{.value}_{htn}"
    ) %>% 
    # calculate RD and RR
    mutate(
      # RD is on % scale
      RD_est = (est_1 - est_0) * 100, 
      RD_se = sqrt(var_1 + var_0) * 100, # approximation; gives the same result as deltamethod
      RR_est = est_1 / est_0,
      RR_se = NA # initialize; fill with deltamethod() below
    ) 
  
  for (i in 1:nrow(RD_RR_results)) {
    RD_RR_results[i, "RR_se"] <- deltamethod(
      ~ x1 / x2, mean = c(RD_RR_results$est_1[i], RD_RR_results$est_0[i]), 
      cov = diag(c(RD_RR_results$var_1[i], RD_RR_results$var_0[i]))
    )
  }
  
  RD_RR_results <- RD_RR_results %>% 
    mutate(
      # suppress RD estimates before age of first dementia case
      across(
        -c(time, starts_with("est"), starts_with("var"), RR_est, RR_se), 
        function(x) ifelse(is.nan(RR_est), NaN, x)
      ),
      RD_lower = RD_est - 1.96 * RD_se, 
      RD_upper = RD_est + 1.96 * RD_se, 
      RR_lower = RR_est - 1.96 * RR_se,
      RR_upper = RR_est + 1.96 * RR_se
    ) 
  
  cum_inc_results <- RD_RR_results %>% 
    select(time, starts_with("est_"), starts_with("var_"))
  
  RD_RR_results <- RD_RR_results %>% 
    select(time, RD_est, RD_lower, RD_upper, RR_est, RR_lower, RR_upper) %>% 
    pivot_longer(
      cols = -time, 
      names_to = c("type", ".value"),
      names_sep = "_"
    )

  return(list("cum_inc" = cum_inc_results, "RD_RR" = RD_RR_results))
}


# test function ----

# test <- predict_cumulative_inc(
#   data = imputed_tte_data_5,
#   ethn = 2, 
#   test_subj = imputed_tte_data_5[[1]] %>% slice_head(), 
#   coxformula = "Surv(survey_age, end_age, event_dem) ~ htn_dx5yr_flag + female + edu_3 + usaborn_rev + ehr_ht_median", 
#   longdata = FALSE
# )
# 
# test$cum_inc %>% 
#   ggplot(aes(x = time)) + 
#   geom_line(aes(y = est_0, color = "no hypertension")) + 
#   geom_line(aes(y = est_1, color = "hypertension")) 
# 
# test$RD_RR %>% 
#   mutate(hliney = ifelse(type == "RD", 0, 1)) %>% 
#   ggplot(aes(x = time, y = est)) + 
#   geom_line() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + 
#   geom_hline(aes(yintercept = hliney), linetype = "dashed", color = "grey30") + 
#   labs(x = "Age (years)", y = NULL) + 
#   facet_wrap(
#     ~ type, scales = "free_y",
#     strip.position = "left", 
#     labeller = as_labeller(c(RD = "Risk difference (%)", RR = "Risk ratio"))
#   ) +
#   theme_bw() + 
#   theme(
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     strip.text = element_text(size = 11)
#   ) 
# 
# testlong <- predict_cumulative_inc(
#   data = imputed_tte_data_long,
#   ethn = 2, 
#   test_subj = imputed_tte_data_5[[1]] %>% 
#       slice_head() %>%
#       select(female, edu_3, usaborn_rev, ehr_ht_median) %>%
#       mutate(newid = 1) %>%
#       right_join(expand_grid(age_start = 60:89, newid = 1), by = "newid") %>%
#       mutate(
#         current_age_c75 = age_start - 75,
#         age_end = age_start + 1,
#         event_dem = 0
#       ),
#   coxformula = "Surv(age_start, age_end, event_dem) ~ htn_dx5yr_flag + htn_age_int + female + edu_3 + usaborn_rev + ehr_ht_median",
#   longdata = TRUE
# )
# 
# testlong$cum_inc %>% 
#   ggplot(aes(x = time)) + 
#   geom_line(aes(y = est_0, color = "no hypertension")) + 
#   geom_line(aes(y = est_1, color = "hypertension")) 
# 
# testlong$RD_RR %>% 
#   mutate(hliney = ifelse(type == "RD", 0, 1)) %>% 
#   ggplot(aes(x = time, y = est)) + 
#   geom_line() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + 
#   geom_hline(aes(yintercept = hliney), linetype = "dashed", color = "grey30") + 
#   labs(x = "Age (years)", y = NULL) + 
#   facet_wrap(
#     ~ type, scales = "free_y",
#     strip.position = "left", 
#     labeller = as_labeller(c(RD = "Risk difference (%)", RR = "Risk ratio"))
#   ) +
#   theme_bw() + 
#   theme(
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     strip.text = element_text(size = 11)
#   ) 

# run function on different ethnicities ----

# prediction is based on a hypothetical subject with 
# the median/most common value for each covariate included in the model

imputed_tte_data_5[[1]] %>% 
  select(female, edu_3, usaborn_rev, ehr_ht_median) %>% 
  summary()

test_subj <- tibble(
  female = factor(1), 
  edu_3 = factor(
    "HS + Some college", 
    levels = c("< HS", "HS + Some college", "College and above")
  ), 
  usaborn_rev = 1, 
  ehr_ht_median = 66
)

test_subj_long <- test_subj %>% 
  mutate(newid = 1) %>%
  right_join(expand_grid(age_start = 60:89, newid = 1), by = "newid") %>%
  mutate(
    current_age_c75 = age_start - 75,
    age_end = age_start + 1,
    event_dem = 0
  )

ethnicities <- levels(imputed_tte_data_5[[1]]$ethnicity_rev)

no_int_results <- list()
int_results <- list()

# looping through each ethnicity, running all 4 models using fx
for (i in seq_along(ethnicities)) {
  # i <- 2
  no_int_results[[i]] <- predict_cumulative_inc(
    data = imputed_tte_data_5,
    ethn = ethnicities[i], 
    test_subj = test_subj, 
    coxformula = "Surv(survey_age, end_age, event_dem) ~ htn_dx5yr_flag + female + edu_3 + usaborn_rev + ehr_ht_median", 
    longdata = FALSE
  )
  
  int_results[[i]] <- predict_cumulative_inc(
    data = imputed_tte_data_long,
    ethn = ethnicities[i], 
    test_subj = test_subj_long,
    coxformula = "Surv(age_start, age_end, event_dem) ~ htn_dx5yr_flag + htn_age_int + female + edu_3 + usaborn_rev + ehr_ht_median",
    longdata = TRUE
  )

}

names(no_int_results) <- ethnicities
names(int_results) <- ethnicities

# save(
#   no_int_results, int_results,
#   file = here("scripts", "output", "cumulative_inc_RD_RR.RData")
# )

# plot cumulative incidence and RD ----

load(here("scripts", "output", "cumulative_inc_RD_RR.RData"))

## model without age interaction ----
# this model is not presented in the final manuscript
lapply(no_int_results, function(x) x$cum_inc) %>% 
  bind_rows(.id = "ethnicity") %>% 
  mutate(
    ethnicity = factor(ethnicity, levels = c(2, 3, 5, 1, 9),
                       labels = c("Chinese", "Japanese",
                                  "Filipino", "South Asian", "White"))
  ) %>% 
  ggplot(aes(x = time)) + 
  geom_line(aes(y = est_0, color = "no hypertension")) + 
  geom_line(aes(y = est_1, color = "hypertension")) + 
  labs(x = "Age (years)", y = "Cumulative incidence", color = NULL,
       subtitle = "Fully adjusted model without age interaction") + 
  facet_grid(cols = vars(ethnicity)) + 
  # facet_wrap(~ethnicity) + 
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    legend.position = "bottom",
    # legend.position.inside = c(0.85, 0.2)
  )

# ggsave(
#   here("scripts", "output", "figures", "supp_cumulative_inc_no_ageint.png")
# )


lapply(no_int_results, function(x) x$RD_RR) %>% 
  bind_rows(.id = "ethnicity") %>% 
  mutate(
    ethnicity = factor(ethnicity, levels = c(2, 3, 5, 1, 9),
                       labels = c("Chinese", "Japanese",
                                  "Filipino", "South Asian", "White")),
    hliney = ifelse(type == "RD", 0, 1),
    type = factor(type, levels = c("RD", "RR"), 
                  labels = c("Risk difference (%)", "Risk ratio"))
  ) %>% 
  filter(type == "Risk difference (%)") %>% 
  ggplot(aes(x = time, y = est)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  geom_hline(aes(yintercept = hliney), linetype = "dashed", color = "grey30") + 
  labs(x = "Age (years)", y = "Risk difference (%)"
       # subtitle = "Fully adjusted model without age interaction"
       ) + 
  facet_grid(cols = vars(ethnicity), scales = "free_y", switch = "y") + 
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 11)
  )

# ggsave(
#   here("scripts", "output", "figures", "supp_RD_no_ageint.png")
# )

## model with age interaction ----

lapply(int_results, function(x) x$cum_inc) %>% 
  bind_rows(.id = "ethnicity") %>% 
  mutate(
    ethnicity = factor(ethnicity, levels = c(2, 3, 5, 1, 9),
                       labels = c("Chinese", "Japanese",
                                  "Filipino", "South Asian", "White"))
  ) %>% 
  ggplot(aes(x = time)) + 
  geom_line(aes(y = est_0 * 100, color = "no hypertension")) + 
  geom_line(aes(y = est_1 * 100, color = "hypertension")) + 
  labs(x = "Age (years)", y = "Cumulative incidence (%)", color = NULL,
       # subtitle = "Fully adjusted model without age interaction"
       ) + 
  facet_grid(cols = vars(ethnicity)) + 
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    legend.position = "bottom",
    # legend.position.inside = c(0.85, 0.2)
  )

ggsave(
  here("scripts", "output", "figures", "supp_cumulative_inc_w_ageint.png")
)

lapply(int_results, function(x) x$RD_RR) %>% 
  bind_rows(.id = "ethnicity") %>% 
  mutate(
    ethnicity = factor(ethnicity, levels = c(2, 3, 5, 1, 9),
                       labels = c("Chinese", "Japanese",
                                  "Filipino", "South Asian", "White")),
    hliney = ifelse(type == "RD", 0, 1),
    type = factor(type, levels = c("RD", "RR"), 
                  labels = c("Risk difference (%)", "Risk ratio"))
  ) %>% 
  filter(type == "Risk difference (%)") %>% 
  ggplot(aes(x = time, y = est)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + 
  geom_hline(aes(yintercept = hliney), linetype = "dashed", color = "grey30") + 
  labs(x = "Age (years)", y = "Risk difference (%)"
       # subtitle = "Fully adjusted model without age interaction"
  ) + 
  facet_grid(cols = vars(ethnicity), scales = "free_y", switch = "y") + 
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 11)
  )

ggsave(
  here("scripts", "output", "figures", "supp_RD_w_ageint.png")
)

# clean up tabular results for RR and RD ----

raw_table <- lapply(int_results, function(x) x$RD_RR) %>% 
  bind_rows(.id = "ethnicity") %>% 
  filter(time %in% seq(60, 90, by = 5)) %>% 
  mutate(
    ethnicity = factor(ethnicity, levels = c(2, 3, 5, 1, 9),
                       labels = c("Chinese", "Japanese",
                                  "Filipino", "South Asian", "White")),
    across(c(est, lower, upper), 
           function(x) format(round(x, 2), nsmall = 2) %>% str_trim())
  ) 

RD_table <- raw_table %>% filter(type == "RD") %>% 
  mutate(
    RD = paste0(est, " (", lower, ", ", upper, ")"),
    RD = ifelse(est == "NaN", "--", RD)
  ) %>% 
  rename(age = time) %>% 
  select(ethnicity, age, RD) %>% 
  pivot_wider(
    id_cols = age, 
    names_from = ethnicity, 
    values_from = RD
  ) 

RR_table <- raw_table %>% filter(type == "RR") %>% 
  mutate(
    # RR = paste0(est, " (", lower, ", ", upper, ")"),
    RR = ifelse(est == "NaN", "--", est)
  ) %>% 
  rename(age = time) %>% 
  select(ethnicity, age, RR) %>% 
  pivot_wider(
    id_cols = age, 
    names_from = ethnicity, 
    values_from = RR
  ) 

# save the tables
write.xlsx(
  list(RD = RD_table, RR = RR_table), 
  here("scripts", "output", "tables", "RD_RR_estimates.xlsx")
)

# checking dementia ages
imputed_tte_data_5[[1]] %>% 
  # filter(ethnicity_rev == 1) %>% 
  filter(main_dem_v1_end_dem_flag == 1) %>% 
  group_by(ethnicity_rev, htn_dx5yr_flag) %>% 
  summarise(min = min(main_dem_v1_end_age, na.rm = TRUE)) 
# I think we definitely should not show estimates for age 60
# then for the Japanese group, also don't show age 65






