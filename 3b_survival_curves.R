# Survival curves/cumulative incidence
# Created by: Natalie Gradwohl
# Creation date: 2/14/2024

# These figures are included in the supplement

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load(
  "tidyverse", "here", "survival", "broom"
  # "survminer", "ggsurvfit",  "tidycmprsk",
  # "cowplot", "ggsci"
)

options(scipen = 9999)

# ----loading data----
source(here("scripts", "0_paths.R"))

load(paste0(path_to_analytic_data, "imputed_tte_data_5.rds"))
data <- imputed_tte_data_5[[1]]
# we are only using complete tte info, so we just take one dataset


# KM curves and plot ----

KM <- survfit(
  Surv(survey_age, main_dem_v1_end_age, main_dem_v1_end_dem_flag) ~ 
    strata(ethnicity_rev, htn_dx5yr_flag), data = data
)

KM %>% tidy() %>% 
  filter(time <= 90) %>% 
  mutate(
    cif = 100 * (1 - estimate),   
    strata = str_sub(strata, start = 53),
    htn = ifelse(str_detect(strata, "0"), "No hypertension", "Hypertension") %>% 
      factor(levels = c("Hypertension", "No hypertension")),
    ethnicity = str_split_i(strata, ",", 1) %>% 
      factor(levels = c("Chinese", "Filipino", "Japanese", 
                        "South Asian", "White"),
             labels = c("Chinese", "Filipino", "Japanese", 
                        "South Asian", "Non-Latino White"))
  ) %>%
  ggplot(aes(x = time, y = cif, group = htn, color = htn)) + 
  geom_line() + 
  scale_y_continuous(
    breaks = seq(0, 40, 10), minor_breaks = seq(0, 40, 5)
  ) + 
  facet_grid(cols = vars(ethnicity)) +
  theme_bw() + 
  # scale_color_discrete(breaks = c("Hypertension", "No hypertension")) + 
  theme(
    legend.position = "bottom",
    # legend.position.inside = c(0.85, 0.3)
  ) +
  labs(
    x = "Age (years)",
    y = "Cumulative incidence (%)", 
    color = NULL, 
    fill = NULL
  ) 
  
ggsave(
  width = 10, height = 4,
  file = here("scripts", "output", "figures", "KM_curve.png")
)

# cause specific cum inc and plot ----

unique(data$main_dem_v1_end_type)

# create proper event end types for competing risks
data <- data %>% 
  mutate(
    end_type_3 = case_when(
      main_dem_v1_end_type %in% 
        c("ADMIN CENSORED", "END OF MEMBERSHIP", "CENSORED 90+") ~ 0,
      main_dem_v1_end_type == "DEMENTIA" ~ 1,  
      main_dem_v1_end_type == "DEATH" ~ 2
    ) %>% 
      factor(levels = 0:2, labels = c("censored", "dementia", "death"))
  )

cause_spec <- survfit(
  Surv(survey_age, main_dem_v1_end_age, end_type_3) ~ 
    strata(ethnicity_rev, htn_dx5yr_flag), data = data
)

cause_spec %>% tidy() %>%
  filter(state != "(s0)", time <= 90) %>% 
  mutate(
    strata = str_sub(strata, start = 53),
    htn = ifelse(str_detect(strata, "0"), "No hypertension", "Hypertension") %>% 
      factor(levels = c("No hypertension", "Hypertension")),
    ethnicity = str_split_i(strata, ",", 1) %>% 
      factor(levels = c("Chinese", "Filipino", "Japanese", 
                        "South Asian", "White"),
             labels = c("Chinese", "Filipino", "Japanese", 
                        "South Asian", "Non-Latino White"))
  ) %>% 
  group_by(htn, ethnicity, time) %>% 
  mutate(cum_sum = sum(estimate)) %>% 
  ungroup() %>% 
  filter(state != "death") %>% 
  select(ethnicity, htn, time, estimate, cum_sum) %>% 
  rename(dementia = estimate, death = cum_sum) %>% 
  mutate(
    across(c(dementia, death), function(x) x * 100) 
  ) %>% 
  ggplot(aes(x = time)) + 
  geom_ribbon(
    aes(ymin = 0, ymax = dementia, color = "Dementia", fill = "Dementia"),
    alpha = 0.6, linetype = 0
  ) + 
  geom_ribbon(
    aes(
      ymin = dementia, ymax = death, 
      color = "Death without dementia", fill = "Death without dementia"
    ),
    alpha = 0.6, linetype = 0
  ) + 
  geom_ribbon(
    aes(
      ymin = death, ymax = 100, 
      color = "Alive without dementia", fill = "Alive without dementia"
    ),
    alpha = 0.6, linetype = 0
  ) + 
  scale_y_continuous(
    breaks = seq(0, 100, 20), minor_breaks = seq(10, 90, 20)
  ) + 
  scale_fill_viridis_d() + 
  theme_bw() + 
  # either this
  # facet_grid(cols = vars(htn), rows = vars(ethnicity)) +
  # or this
  facet_grid(rows = vars(htn), cols = vars(ethnicity)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Age (years)",
    y = "Cumulative incidence (%)", 
    color = NULL, 
    fill = NULL
  ) 

ggsave(
  width = 6, height = 8,
  file = here("scripts", "output", "figures", "cause_spec_cum_inc.png")
)

ggsave(
  width = 8, height = 4,
  file = here("scripts", "output", "figures", "cause_spec_cum_inc_alt.png")
)

# OLD CODE ----
# creating multiple event flag for cumulative icidence

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
