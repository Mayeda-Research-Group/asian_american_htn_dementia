# Purpose : Cox models looking at association between hypertension and dementia
# using age as timescale and adjusting for time-varying age
# Created by: Natalie Gradwohl
# Creation date: 2/13/2024

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load(
  "tidyverse", "here", "broom", 
  "survival", 
  "mice", "mitools", 
  "purrr", "openxlsx",
  "RColorBrewer",
  "flextable"
)

options(scipen = 9999)

# ----loading data and function (separate script in project)----
source(here("scripts", "0_paths.R"))

# CODE SAMPLE NOTE: this function is attached as a separate script
source(here("scripts", "cox_multiple_models_function.R"))

load(paste0(path_to_analytic_data, "imputed_tte_data_5.rds"))

# ---- create long dataset with current age (time-varying age)----

# look at the tte variables to use
imputed_tte_data_5[[1]] %>% names()
imputed_tte_data_5[[1]] %>% select(subjid, ethnicity_rev, starts_with("main_dem")) %>% View()

# the Cox model will use Surv(survey_age, main_dem_v1_end_age, main_dem_v1_end_dem_flag)
# main_dem_v1_end_age: age at end of followup, either due to dementia or other events
# main_dem_v1_end_dem_flag: whether the end of followup is dementia (1) or not (0)

summary(imputed_tte_data_5[[1]]$main_dem_v1_end_age)

data_cut <- function(data) {
  cuts <- c(
    min(data$survey_age), # the minimum start age
    seq(60, 102, 1),     # all the intervals in the middle
    max(data$main_dem_v1_end_age)  # the maximum end of follow-up age
  )  
  
  longdata <- survSplit(
    data    = data,
    formula = Surv(survey_age, main_dem_v1_end_age, main_dem_v1_end_dem_flag) ~ .,
    cut     = cuts,
    # renaming the variables that indicate the start and end of each interval
    start   = "age_start",
    end     = "age_end"
  ) %>%
    mutate(current_age_c75 = age_start - 75)
  
  return(longdata)
}

imputed_tte_data_long <- lapply(imputed_tte_data_5, data_cut)

# check to make sure it looks right
# imputed_tte_data_long[[1]] %>% View()

# ----setting model specifications for Cox model----
exposure <- "htn_dx5yr_flag"

age_surv_long <- "Surv(age_start, age_end, main_dem_v1_end_dem_flag)" 

ethnicities <- levels(imputed_tte_data_5[[1]]$ethnicity_rev)

covariates <- c(
  "htn_dx5yr_flag:current_age_c75", #interaction term to allow for effect of htn to vary by age
  "female + edu_3",
  "usaborn_rev + ehr_ht_median",
  "eversmoke + diab_dx5yr_flag"
)

nested_models <- c(
  exposure, 
  lapply(
    seq_along(covariates), 
    function(i) {
      paste(c(exposure, covariates[1:i]), collapse = " + ")
    }
  )
)

# load(here("scripts", "output", "3c_end90_cox_results.RData"))
# 
# ja4 <- multiple_models(
#   data = imputed_tte_data_long,
#   predictors = nested_models[4],
#   ethnicity = "Japanese",
#   response = eval(age_surv_long)
# )
# 
# sa4 <- multiple_models(
#   data = imputed_tte_data_long,
#   predictors = nested_models[4],
#   ethnicity = "South Asian",
#   response = eval(age_surv_long)
# )


# creating an empty list to store results
all_model_results <- list()
age_spec_hrs <- list()


# ----running each of the models for each ethnicity----
for (ethn in ethnicities) {
  # ethn <- ethnicities[1]
  
  model_results <- list()
  result <- list()
  
  # applying multiple models function (from separate script!)
  for (j in seq_along(nested_models)) {
    result[[j]] <- multiple_models(
      data = imputed_tte_data_long,
      predictors = nested_models[j],
      ethnicity = ethn,
      response = eval(age_surv_long)
    )
    
    # save hr summary from each model
    model_results[[j]] <- result[[j]]$hr_summary %>% 
      mutate(model_number = paste("model", j, sep = ""))
  }
  
  # collapse the hr summary into one tibble
  all_model_results[[ethn]] <- do.call(rbind, model_results)
  
  # calculate the age-specific hrs for this ethnicity
  # for models that include age-interaction (2:5)
  
  combined_age_hr <- lapply(
    2:length(nested_models), # model numbers to loop over
    function(x) {
      # obtain coef and vcov from each model
      coef_summary <- result[[x]]$coef_summary
      # loop over selected ages to calculate age-spec HRs
      lapply(
        c(60, 65, 70, 75, 80, 85, 90), 
        function(y) {
          age_hr(y, coef_summary$coefficients, coef_summary$variance) %>% 
            mutate(age = y)
        }
      ) %>% 
        bind_rows() %>% 
        mutate(model = x)
    }
  ) %>% bind_rows() %>% 
    mutate(ethnicity = ethn)
  
  # storing age-spec hrs for each ethnicity
  age_spec_hrs[[ethn]] <- combined_age_hr
}

# save results needed for subsequent steps in case 
# re-running analysis is not possible

# save(
#   age_spec_hrs, age_results_long,
#   file = here("scripts", "output", "3c_cox_model_results.RData")
# )

# save(
#   all_model_results, age_spec_hrs, 
#   file = here("scripts", "output", "3c_cox_model_results_addcovariates.RData")
# )

load(
  here("scripts", "output", "3c_cox_model_results_addcovariates.RData")
)

# output table for fully-adjusted HRs for each ethnicity----
# these HRs are for someone age 75 (interaction term is centered at 75)
hr_tables <- bind_rows(all_model_results) %>% 
  filter(term == "htn_dx5yr_flag", model_number != "model3") %>%
  select(ethnicity, model_number, estimate_table) %>%
  pivot_wider(names_from = model_number, values_from = estimate_table) %>% 
  mutate(
    ethnicity = factor(ethnicity, levels = ethnicities, 
                       labels = c(ethnicities[1:4], "Non-Latino White"))
  ) %>% 
  arrange(ethnicity) %>% 
  flextable()
  

hr_tables_w_int <- bind_rows(all_model_results) %>% 
  filter(str_detect(term, "htn_dx5yr_flag"), model_number != "model3") %>% 
  mutate(
    HR = ifelse(term == "htn_dx5yr_flag", "main effect", "age-interaction") 
  ) %>% 
  select(ethnicity, model_number, HR, estimate_table) %>%
  pivot_wider(
    names_from = c(model_number, HR), 
    values_from = estimate_table) %>% 
  mutate(
    ethnicity = factor(ethnicity, levels = ethnicities, 
      labels = c(ethnicities[1:4], "Non-Latino White"))
  ) %>% 
  arrange(ethnicity) %>% 
  flextable() %>% 
  separate_header(split = "_")

# output figure and table for age-spec hrs----

# combining age-spec results
stacked_df_age <- bind_rows(age_spec_hrs) 

# adding data for where to place CI's + arrows for when CI is too large for chart
forest_plot_data <- stacked_df_age %>%
  filter(age != 60) %>% 
  mutate(
    ethnicity = factor(ethnicity, levels = ethnicities, 
                       labels = c(ethnicities[1:4], "Non-Latino White")),
    across(c(hr, lower_ci, upper_ci), 
           function(x) ifelse(age == 65 & ethnicity == "Japanese", NA_real_, x)),
    
    # conf.high.ceiling = ifelse(upper_ci <= 3.6, upper_ci, NA),
    # conf.low.ceiling = ifelse(lower_ci >= 0.25, lower_ci, NA),
    # upper_arrow_pos = ifelse(upper_ci > 3.6, 3.6, NA),
    # lower_arrow_pos = ifelse(lower_ci < 0.25, 0.3, NA),
    
    conf.high.ceiling = ifelse(upper_ci <= 2.5, upper_ci, NA),
    conf.low.ceiling = ifelse(lower_ci >= 0.25, lower_ci, NA),
    upper_arrow_pos = ifelse(upper_ci > 2.5, 2.5, NA),
    lower_arrow_pos = ifelse(lower_ci < 0.25, 0.3, NA),
    
    # markpos = ifelse(age == 65 & ethnicity == "Japanese", 1.25, NA_real_), 
    markpos = ifelse(
      (age == 65 & ethnicity != "Non-Latino White") | (age == 70 & ethnicity == "South Asian"), 
      1.25, NA_real_
    ),
    across(
      -c(age, model, ethnicity, markpos), function(x) {
        ifelse(
          (age == 65 & ethnicity != "Non-Latino White") | (age == 70 & ethnicity == "South Asian"), 
          NA_real_, x
        )
      }
    )
  )


# figure for qge-specific (calculated every 5 years), ethnicity-specific hazard ratios

age_spec_plot_all_mods <- forest_plot_data %>% 
  ggplot(aes(x = age, y = hr, ymin = 0)) +
  geom_point(size = 2) +
  geom_segment(
    aes(x = age, xend = age, y = lower_ci, yend = upper_arrow_pos),
    arrow = arrow(length = unit(0.2, "cm")),
    show.legend = FALSE
  ) +
  # geom_segment(
  #   aes(x = age, xend = age, y = upper_ci, yend = lower_arrow_pos),
  #   arrow = arrow(length = unit(0.2, "cm")),
  #   show.legend = FALSE
  # ) +
  # geom_segment(
  #   aes(x = age, xend = age, y = upper_arrow_pos, yend = lower_arrow_pos),
  #   arrow = arrow(length = unit(0.2, "cm"), ends = "both"),
  #   show.legend = FALSE
  # ) +
  geom_text(
    aes(x = age, y = upper_arrow_pos + 0.1, label = round(upper_ci, 2)),
    size = 3.5
  ) +
  # geom_text(
  #   aes(x = age, y = lower_arrow_pos - 0.1, label = round(lower_ci, 2)),
  #   size = 3.5
  # ) +
  geom_errorbar(
    aes(ymin = conf.low.ceiling, ymax = conf.high.ceiling), width = 1
  ) +
  geom_point(
    aes(y = markpos), shape = 4, show.legend = FALSE
  ) + 
  facet_grid(cols = vars(ethnicity), rows = vars(model)) +
  #scale_color_manual(values = c("#E69F00", "#0072B2")) +
  guides(color = "none") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  scale_y_continuous(
    breaks = seq(0.25, 3.75, by = 0.25),
    limits = c(0.25, 3.8),
    expand = c(0.01, 0.01)
  ) +
  scale_x_continuous(limits = c(63, 92), breaks = seq(65, 90, 5)) + 
  labs(
    # title = "Age-specific hazard ratios for dementia incidence (hypertension vs. no hypertension at ≥5 years prior to baseline)",
    y = "Hazard ratio (95% CI)", x = "Age (years)"
  ) + 
  theme_bw() +
  theme(
    text = element_text(size = 12),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) 

age_spec_plot_final_mod <- forest_plot_data %>% 
  filter(model == 5) %>% 
  ggplot(aes(x = age, y = hr, ymin = 0)) +
  geom_point(size = 2) +
  geom_segment(
    aes(x = age, xend = age, y = lower_ci, yend = upper_arrow_pos),
    arrow = arrow(length = unit(0.2, "cm")),
    show.legend = FALSE
  ) +
  # geom_segment(
  #   aes(x = age, xend = age, y = upper_ci, yend = lower_arrow_pos),
  #   arrow = arrow(length = unit(0.2, "cm")),
  #   show.legend = FALSE
  # ) +
  # geom_segment(
  #   aes(x = age, xend = age, y = upper_arrow_pos, yend = lower_arrow_pos),
  #   arrow = arrow(length = unit(0.2, "cm"), ends = "both"),
  #   show.legend = FALSE
  # ) +
  geom_text(
    aes(x = age, y = upper_arrow_pos + 0.1, label = round(upper_ci, 2)),
    size = 3.2
  ) +
  # geom_text(
  #   aes(x = age, y = lower_arrow_pos - 0.1, label = round(lower_ci, 2)),
  #   size = 3.2
  # ) +
  geom_errorbar(
    aes(ymin = conf.low.ceiling, ymax = conf.high.ceiling), width = 1
  ) +
  geom_point(
    aes(y = markpos), shape = 4, show.legend = FALSE
  ) + 
  facet_grid(cols = vars(ethnicity)) +
  #scale_color_manual(values = c("#E69F00", "#0072B2")) +
  guides(color = "none") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  # scale_y_continuous(
  #   breaks = seq(0.25, 3.75, by = 0.25),
  #   limits = c(0.25, 3.8),
  #   expand = c(0.01, 0.01)
  # ) +
  scale_y_continuous(
    breaks = seq(0.25, 2.5, by = 0.25),
    limits = c(0.4, 2.7),
    expand = c(0.01, 0.01)
  ) +
  scale_x_continuous(limits = c(63, 92), breaks = seq(65, 90, 5)) + 
  labs(
    # title = "Age-specific hazard ratios for dementia incidence (hypertension vs. no hypertension at ≥5 years prior to baseline)",
    y = "Hazard ratio", x = "Age (years)"
  ) + 
  theme_bw() +
  theme(
    text = element_text(size = 11),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) 

# table for supplement (hazard ratios every 5 years for each ethnicity as in figure)
age_spec_tables <- stacked_df_age %>%
  filter(age != 60) %>% 
  mutate(
    across(c(hr, lower_ci, upper_ci), function(x) sprintf("%0.2f", x)), 
    est_ci = paste0(hr, " (", lower_ci, ", " , upper_ci, ")"),
    ethnicity = factor(ethnicity, levels = ethnicities, 
                       labels = c(ethnicities[1:4], "Non-Latino White"))
  ) %>%
  select(model, age, ethnicity, est_ci) %>%
  arrange(ethnicity) %>% 
  mutate(
    est_ci = ifelse(
      (age == 65 & ethnicity != "Non-Latino White") | (age == 70 & ethnicity == "South Asian"),
      "--", est_ci
    )
  ) %>% 
  pivot_wider(names_from = age, values_from = est_ci) %>% 
  filter(model %in% c(4,5)) %>% 
  arrange(model) %>% 
  flextable()


# ----saving outputs----

save_as_docx(
  hr_tables %>% autofit(), 
  hr_tables_w_int %>% autofit(),
  age_spec_tables %>% autofit(), 
  path = here(
    "scripts", "output", "tables",
    "cox_age_current_age_results_5yr_imp.docx"
  )
)

# age-spec, ethnicity-spec HR figure (main figure for manuscript)
ggsave(
  age_spec_plot_all_mods,
  width = 10, height = 12,
  file = here("scripts", "output", "figures", "age_spec_hr_plot_all_mods.png")
)

ggsave(
  age_spec_plot_final_mod,
  # width = 8, height = 4,
  width = 8, height = 3, 
  file = here("scripts", "output", "figures", "age_spec_hr_plot_final_mod.png")
)


# OLD ----
# plot of age-specific hazard ratios

# PLOT: Hazard ratios by ethnicity, across ages 60-90, with confidence bands
plot <- stacked_df_age %>%
  ggplot(aes(x = age, y = hr, group = ethnicity)) +
  geom_line(aes(color = ethnicity)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = ethnicity), alpha = 0.3) +
  facet_wrap( ~ ethnicity, ncol = 5) +
  labs(title = "Hazard Ratio by Ethnicity",
       x = "Age",
       y = "HR") +
  guides(fill = guide_legend(title = "Ethnicity"), color = "none") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
