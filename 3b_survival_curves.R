# Survival curves/cumulative incidence
# Created by: Natalie Gradwohl
# Creation date: 2/14/2024

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load(
  "tidyverse",
  "here",
  "survival",
  "survminer",
  "ggsurvfit", 
  "tidycmprsk",
  "cowplot",
  "ggsci"
)

options(scipen = 9999)

# ----loading data----
source(here("scripts", "0_paths.R"))

data <-
  readRDS(paste0(path_to_analytic_data, "wide_data_pre_mi_5.rds"))

# creating multiple event flag for cumulative inidence
data$event_multi <-
  factor(
    ifelse(
      data$event_dem == 1 & data$event_death == 0,
      1,
      ifelse(data$event_death == 1, 2, 0)
    ),
    labels = c("s(0)", "dementia", "death")
  )

table(data$end_type, data$event_multi)

eths <- levels(data$ethnicity_rev)

# surv objects
km_surv <- "Surv(survey_age, end_age, event_dem)"
cum_inc_surv <- "Surv(survey_age, end_age, event_multi)"

# ----function necessary to save ggsurvplot objects-----
grid.draw.ggsurvplot <- function(x) {
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

# ----loop to generate km curve and cumulative incidence for each ethnicity----

curves <- list()
fit <- list()

for (i in seq_along(eths)) {
  # subset the data for the specific race
  subset_data <- data[data$ethnicity_rev == eths[i],]
  
  # fit the survival model
  curves[[i]] <-
    survfit(formula(paste0(km_surv, " ~ htn_dx5yr_flag")),
            data = subset_data,
            id = subjid)
  
  # plot the survival curve
  km <- ggsurvplot(
    curves[[i]],
    data = subset_data,
    conf.int = TRUE,
    risk.table = FALSE,
    break.time.by = 5,
    xlim = c(60, max(data$end_age)),
    title = paste0("KM curve for dementia-free survival: ", eths[i]),
    #legend.labs = c("Non-hypertensive", "Hypertensive"),
    xlab = "Age (years)"
  )
  
  print(km)
  
  ggsave(
    filename = here(
      "scripts",
      "output",
      "figures",
      paste0(eths[i], "_km_dementia.png")
    ),
    km,
    dpi = 300,
    width = 10,
    height = 7,
    units = "in"
  )
  
  # multiple outcomes
  fit[[i]] <-
    survfit(formula(paste0(cum_inc_surv, " ~ htn_dx5yr_flag")), id = subjid,
            data = subset_data_new_id)
  
  
  # plot cumulative incidence functions
  comp_risk <- ggcompetingrisks(fit[[i]],
                                title = paste("Cumulative Incidence Functions, Ethnicity:", eths[i])) +
    theme_cowplot() +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    scale_fill_jco()
  
  print(comp_risk)
  
  ggsave(
    filename = here(
      "scripts",
      "output",
      "figures",
      "survival_curves",
      paste0(eths[i], "_cuminc_dementia.png")
    ),
    comp_risk,
    dpi = 300,
    width = 10,
    height = 7,
    units = "in"
  )
  
}
