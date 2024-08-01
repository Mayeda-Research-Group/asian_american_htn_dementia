# Purpose : Cox models looking at association between hypertension and dementia
# using age as timescale and adjusting for time-varying age
# Created by: Natalie Gradwohl
# Creation date: 2/13/2024

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load(
  "tidyverse",
  "here",
  "broom",
  "survival",
  "survminer",
  "mice",
  "mitools",
  "purrr",
  "openxlsx",
  "RColorBrewer"
)

options(scipen = 9999)

# ----loading data and function (separate script in project)----
source(here("scripts", "0_paths.R"))

# CODE SAMPLE NOTE: this function is attached as a separate script
source(here("scripts", "cox_multiple_models_function.R"))

load(paste0(path_to_analytic_data, "imputed_tte_data_5.rds"))

# ----function to create long dataset with current age (time-varying age)----
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
    mutate(current_age_c75 = age_start - 75)
  
  return(longdata)
}

# ----applying long data function to imputed datasets----
imputed_tte_data_long <- list()

for (i in seq_along(imputed_tte_data_5)) {
  data <- data_cut(imputed_tte_data_5[[i]])
  imputed_tte_data_long[[i]] <- data
}

# ----setting model specifications for Cox model----
exposure <- "htn_dx5yr_flag"

age_surv_long <- "Surv(age_start, age_end, event_dem)" # accounting for left-truncation

ethnicities <- levels(data$ethnicity_rev)

ethnicity_names <-
  c("South Asian", "Chinese", "Japanese", "Filipino", "White")

covariates <- c("htn_dx5yr_flag:current_age_c75", #interaction term to allow for effect of htn to vary by age
                "female + edu_3",
                "usaborn_rev + ehr_ht_median")

nested_models <-
  c(exposure, lapply(seq_along(covariates), function(i) {
    paste(c(exposure, covariates[1:i]), collapse = " + ")
  }))

# creating an empty list to store results
age_results_long <- list()
age_spec_hrs <- list()
result <- list()

# ----running each of the models for each ethnicity----
for (i in seq_along(ethnicities)) {
  ethnicity <- ethnicities[i]
  
  model_results <- list()
  
  # applying multiple models function (from separate script!)
  for (j in seq_along(nested_models)) {
    result[[j]] <- multiple_models(
      data = imputed_tte_data_long,
      predictors = nested_models[j],
      ethnicity = ethnicity,
      response = eval(age_surv_long)
    )
    
    # create a dataframe from the hr summary
    model_df <- result[[j]]$hr_summary
    
    model_df$model_number <- paste("model", j, sep = "")
    
    #storing df
    model_results[[j]] <- model_df
  }
  
  # combining dataframes from all models into single df for each ethnicity
  combined_results <- do.call(rbind, model_results)
  
  # adding combined results to main results list
  age_results_long[[ethnicity]] <- combined_results
  
  # storing the age-specific hrs for this ethnicity
  coef_summary <- result[[length(nested_models)]]$coef_summary
  
  age_hr_list <- list()
  
  # looping through all ages
  for (age in c(60, 65, 70, 75, 80, 85, 90)) {
    age_hr_list[[as.character(age)]] <-
      age_hr(age, coef_summary$coefficients, coef_summary$variance)
  }
  
  # combining all age-spec hrs into 1 dataframe per ethnicity
  combined_age_hr <- do.call(rbind, age_hr_list)
  
  combined_age_hr <- combined_age_hr %>%
    mutate(age = as.numeric(row.names(combined_age_hr)),
           ethnicity = ethnicity_names[i])
  
  # storing age-spec hrs for each ethnicity
  age_spec_hrs[[ethnicity]] <- combined_age_hr
}

# ----plot of age-specific hazard ratios----

# combining age-spec results
stacked_df_age <- bind_rows(age_spec_hrs)

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


# ----output tables for fully-adjusted HRs for each ethnicity----
# these HRs are for someone age 75 (interaction term is centered at 75)
hr_tables <- bind_rows(age_results_long) %>%
  mutate(ethnicity = factor(ethnicity, levels = ethnicities, labels = ethnicity_names)) %>%
  filter(term == "htn_dx5yr_flag", model_number != "model3") %>%
  select(ethnicity, model_number, estimate_table) %>%
  pivot_wider(names_from = model_number, values_from = estimate_table)

# ----output table for age-spec hrs----

# adding data for where to place CI's + arrows for when CI is too large for chart
forest_plot_data <- stacked_df_age %>%
  mutate(
    conf.high.ceiling = ifelse(upper_ci <= 3.6, upper_ci, NA),
    conf.low.ceiling = ifelse(lower_ci >= 0.7, lower_ci, NA),
    upper_arrow_pos = ifelse(upper_ci > 3.6, 3.6, NA),
    lower_arrow_pos = ifelse(lower_ci < 0.7, 0.7, NA)
  )

# creating a theme for the ggplot
tree_plot_theme <- theme_bw() + theme(
  legend.position = c(0.1, 0.85),
  legend.background = element_blank(),
  legend.box.background = element_rect(color = "white"),
  axis.title.y = element_text(vjust = 3),
  #axis.title.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.text.x = element_text(size = 11),
  axis.text.y = element_text(size = 9),
  strip.text = element_text(size = 12),
  plot.title = element_text(hjust = 0.5, vjust = 0.5),
)

# PLOT: Age-specific (calculated every 5 years), ethnicity-specific hazard ratios
age_spec_plot <-
  ggplot(forest_plot_data, aes(
    x = age,
    y = hr,
    ymin = 0,
    color = age
  )) +
  geom_point(aes(y = hr),
             position = position_nudge(x = 0.15),
             size = 2) +
  geom_segment(
    aes(
      x = age,
      xend = age,
      y = lower_ci,
      yend = upper_arrow_pos
    ),
    position = position_nudge(x = forest_plot_data$line_nudge),
    arrow = arrow(length = unit(0.2, "cm")),
    show.legend = FALSE
  ) +
  geom_segment(
    aes(
      x = age,
      xend = age,
      y = upper_ci,
      yend = lower_arrow_pos
    ),
    position = position_nudge(x = forest_plot_data$line_nudge),
    arrow = arrow(length = unit(0.2, "cm")),
    show.legend = FALSE
  ) +
  geom_segment(
    aes(
      x = age,
      xend = age,
      y = upper_arrow_pos,
      yend = lower_arrow_pos
    ),
    position = position_nudge(x = forest_plot_data$line_nudge),
    arrow = arrow(length = unit(0.2, "cm"), ends = "both"),
    show.legend = FALSE
  ) +
  geom_text(aes(
    x = age,
    y = upper_arrow_pos + 0.1,
    label = round(upper_ci, 2)
  ),
  size = 3) +
  geom_text(aes(
    x = age,
    y = lower_arrow_pos - 0.1,
    label = round(lower_ci, 2)
  ),
  size = 3) +
  geom_errorbar(aes(ymin = conf.low.ceiling, ymax = conf.high.ceiling),
                #position = position_nudge(x = forest_plot_data$line_nudge),
                width = 2) +
  facet_grid(. ~ ethnicity) +
  #scale_color_manual(values = c("#E69F00", "#0072B2")) +
  guides(color = "none") +
  geom_hline(aes(yintercept = 1),
             linetype = "dashed") +
  scale_y_continuous(
    breaks = seq(0.5, 3.75, by = 0.25),
    limits = c(0.5, 3.75),
    expand = c(0, 0)
  ) +
  labs(y = "Hazard ratio (95% CI)",
       x = "Age (years)") +
  # labs(color = "Current age") +
  ggtitle(
    "Age-specific hazard ratios for dementia incidence (hypertension vs. no hypertension at ≥5 years prior to baseline)"
  ) +
  tree_plot_theme

# table for supplement (hazard ratios every 5 years for each ethnicity as in figure)
age_spec_tables <- stacked_df_age %>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(est_ci = paste0(hr, " (", lower_ci, ", " , upper_ci, ")")) %>%
  select(age, ethnicity, est_ci) %>%
  pivot_wider(names_from = ethnicity, values_from = est_ci)


# ----saving outputs----

# overall HRs by ethnicity
write.xlsx(
  hr_tables,
  file = here(
    "scripts",
    "output",
    "tables",
    "cox_age_current_age_results_5yr_imp.xlsx"
  )
)

# age-specific HRs by ethnicity
write.xlsx(age_spec_tables,
           file = here("scripts",
                       "output",
                       "tables",
                       "age_spec_hrs_5yr_imp.xlsx"))

# age-spec, ethnicity-spec HR figure (main figure for manuscript)
ggsave(
  age_spec_plot,
  width = 15,
  height = 8,
  file = here("scripts",
              "output",
              "figures",
              "age_spec_hr_plot.png")
)
