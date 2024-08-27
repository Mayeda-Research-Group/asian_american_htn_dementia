# IRs and IRRs
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/13/2024

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("tidyverse", "here", "openxlsx")

options(scipen = 9999)

# ----loading data and functions----
source(here("scripts", "0_paths.R"))
source(here("scripts", "ir_irr_functions.R"))

load(paste0(path_to_analytic_data, "imputed_tte_data_5.rds"))

# since incidence rate calculation does not depend on variables with missingness
# we can use any imputed dataset to carry out the calculation
# however, note that the tte variables are after ending FU at age 90

data <- imputed_tte_data_5[[1]]

# age groups
age_cat <- c(59.99, 65, 70, 75, 80, 85, 120)
# standard population from 2000 census that falls in each age group
std_pop <- c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587)

# note that the last age group is 85+ 
# I set 120 to be the upper limit since we don't have event ages older than that
# the std pop count for the last group is also 85+ 

# ----calculating person-years----

result <- calculate_pys(
  data,
  age_cat,
  fu_start = "survey_age",
  fu_end = "main_dem_v1_end_age",
  event_flag = "main_dem_v1_end_dem_flag"
)

age_cat_labels <- c(
  "60-64 years",
  "65-69 years",
  "70-74 years",
  "75-79 years",
  "80-84 years",
  "85+ years"
)

data <- result$data %>% 
  select(contains(c("py_", "case_")), ethnicity_rev, htn_dx5yr_flag, survey_age)

# ---- age-specific IRs for each ethnicity, hypertension status----
eths <- levels(data$ethnicity_rev)
all_asian <- eths[c("Chinese", "Japanese", "Filipino", "South Asian")]

race_age_spec_irs <- list()

# expect warnings from loop below "! NAs introduced by coercion"
# this is happening because we are introducing NAs instead of calculating IRs where
# there are fewer than 5 cases

for (i in seq_along(eths)) {
  race_age_spec_irs[[paste0(eths[i], "_htn_0")]] <-
    calculate_age_specific_ir(data, age_cat_labels, eths[i], "0") %>% mutate(exposure = "No hypertension")
  race_age_spec_irs[[paste0(eths[i], "_htn_1")]] <-
    calculate_age_specific_ir(data, age_cat_labels, eths[i], "1") %>% mutate(exposure = "Hypertension")
}

## using function to pull IRs into cleaner df----

overall_df <- collect_age_spec_ir(race_age_spec_irs, exposure = "0", eths = eths) %>%
  rbind(collect_age_spec_ir(race_age_spec_irs, exposure = "1", eths = eths)) %>% 
  # label the ethnicities for easier formatting
  mutate(
    ethnicity = factor(
      ethnicity, levels = eths, 
      labels = c(eths[1:4], "Non-Latino White")
    )
  ) %>% 
  arrange(exposure, ethnicity)

overall_cases_pys <- bind_rows(
  collect_age_spec_pys_cases(race_age_spec_irs, exposure = "0", eths = eths),
  collect_age_spec_pys_cases(race_age_spec_irs, exposure = "1", eths = eths)
) %>% 
  # additional formatting
  mutate(
    ethnicity = factor(
      ethnicity, levels = eths, 
      labels = c(eths[1:4], "Non-Latino White")
    )
  ) %>% 
  pivot_wider(
    names_from = age_range,
    values_from = c(pys, cases), 
    names_sep = "_", 
    names_glue = "{age_range}_{.value}"
  ) %>% 
  select(
    ethnicity, exposure, 
    starts_with("60"), starts_with("65"), starts_with("70"), 
    starts_with("75"), starts_with("80"), starts_with("85")
  ) %>% 
  arrange(exposure, ethnicity)



# age-adj irs----

htn_age_adj <- calculate_age_adj_ir(
  data = result$data,
  age_cat_labels = result$age_cat_labels,
  htn = 1,
  std_pop = std_pop,
  per_pys = 1000
) %>% 
  bind_rows(.id = "ethnicity") %>% 
  mutate(exposure = 1)

no_htn_age_adj <- calculate_age_adj_ir(
  data = result$data,
  age_cat_labels = result$age_cat_labels,
  htn = 0,
  std_pop = std_pop,
  per_pys = 1000
) %>% 
  bind_rows(.id = "ethnicity") %>% 
  mutate(exposure = 0)

all_age_adj <- rbind(htn_age_adj, no_htn_age_adj) %>% 
  mutate(
    exposure = factor(
      exposure, levels = c(1, 0), 
      labels = c("Hypertension", "No hypertension")
    ),
    ethnicity = factor(
      ethnicity, levels = eths, 
      labels = c(eths[1:4], "Non-Latino White")
    ),
    across(where(is.numeric), ~ round(., 1))
  )  
  
# append age-adjusted IRs to age-specific IRs table ----

overall_df <- all_age_adj %>% 
  mutate(adj_ir_est = paste0(adj_ir, " (", adj_ir_l, ", ", adj_ir_u, ")")) %>% 
  select(exposure, ethnicity, adj_ir_est) %>% 
  right_join(overall_df, by = c("exposure", "ethnicity")) %>% 
  relocate(adj_ir_est, .after = everything()) %>% 
  arrange(exposure, ethnicity)

# save results ---- 

write.xlsx(
  list(
    "overall age-spec irs" = overall_df,
    "overall cases and pys" = overall_cases_pys
  ),
  here(
    "scripts", "output", "tables",
    "overall_age_stratified_irs.xlsx"
  )
)

# IR figure---

figure_data <- race_age_spec_irs %>% 
  bind_rows(.id = "ethnicity") %>% 
  mutate(
    ethnicity = str_split_i(ethnicity, "_", 1) %>% 
      factor(
        levels = eths,
        labels = c(eths[1:4], "Non-Latino White")
      ),
    exposure = factor(exposure, levels = c("No hypertension", "Hypertension")),
    age_range = str_remove(age_range, " years")
  )

figure_data <- figure_data %>%
  mutate(
    ir = case_when(
      ir == "<5 cases" ~ NA_real_,  
      TRUE ~ as.numeric(ir) 
    ),
    mark_le5 = case_when(
      is.na(ir) ~ 2.5,
      TRUE ~ NA
    ), 
    line_nudge = case_when(
      exposure == "No hypertension" ~ -0.25,
      exposure == "Hypertension" ~ 0.25
    ),
    conf.high.ceiling = ifelse(ir_ul < 65, ir_ul, NA),
    arrow_pos = ifelse(ir_ul > 65, 65, NA)
  )


colors <- c("#2E86C1", "#AED6F1")
# colors <- c("#FFDB6D","#00AFBB")

# IR figure with vertical bars (horizontal commented out)----
age_spec_ir_plot <- figure_data %>% 
  ggplot(aes(x = age_range, y = ir, ymin = 0, group = exposure, fill = exposure)) +
  geom_bar(stat = "identity", position = "dodge", aes(alpha = 0.99)) +
  geom_errorbar(
    aes(ymin = ir_ll, ymax = conf.high.ceiling),
    position = position_nudge(x = figure_data$line_nudge),
    width = 0.15
  ) +
  geom_segment(
    aes(x = age_range, xend = age_range, y = ir_ll,yend = arrow_pos),
    position = position_nudge(x = figure_data$line_nudge),
    arrow = arrow(length = unit(0.2, "cm")),
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = age_range, y = arrow_pos + 1, label = round(ir_ul, 2)),
    position = position_nudge(x = figure_data$line_nudge),
    size = 2.75
  ) +
  # maybe we can add something to mark places with <= 5 events
  geom_point(
    aes(x = age_range, y = mark_le5), shape = 4,
    position = position_nudge(x = figure_data$line_nudge),
    show.legend = FALSE
  ) + 
  scale_y_continuous(
    breaks = seq(0, 65, by = 5),
    limits = c(0, 67.5),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = colors) +
  facet_grid(. ~ ethnicity) +
  labs(
    x = "Age group (years)",
    y = "Incidence rate (per 1,000 PY)",
    # title = "Age- and race/ethnicity-specific incidence rates, stratified by hypertension diagnosis at 5+ years prior to baseline",
    fill = "Exposure"
  ) +
  guides(color = "none") +
  scale_alpha(guide = "none") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggsave(
  age_spec_ir_plot,
  file = here("scripts", "output", "figures", "ir_fig.png"),
  # width = 15, height = 8
  width = 12, height = 5
)
