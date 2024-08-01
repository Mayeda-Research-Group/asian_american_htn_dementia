# IRs and IRRs
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/13/2024

# ----setup and loading packages----
if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("tidyverse",
       "here",
       "openxlsx")

options(scipen = 9999)

# ----loading data and functions----
source(here("scripts", "0_paths.R"))
source(here("scripts", "ir_irr_functions.R"))

load(paste0(path_to_analytic_data, "imputed_tte_data_5.rds"))

data <- imputed_tte_data_5[[1]]

# ----calculating person-years----
age_cat <- c(59.99, 65, 70, 75, 80, 85, 90)

result <-
  calculate_pys(data,
                age_cat,
                "survey_age",
                "main_dem_v1_end_age",
                "main_dem_v1_end_dem_flag")

age_cat_labels <-
  c("60-65 years",
    "65-70 years",
    "70-75 years",
    "75-80 years",
    "80-85 years",
    "85-90+ years")

data <-
  result$data %>% select(contains(c("py_", "case_")), ethnicity_rev, htn_dx5yr_flag, survey_age)

# ----calculating age-specific irs for each ethnicity, hypertension status----
eths <- levels(data$ethnicity_rev)
all_asian <-
  eths[c("Chinese", "Japanese", "Filipino", "South Asian")]

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

# ----using function to pull irs into cleaner df----
overall_df <-
  collect_age_spec_ir(race_age_spec_irs, exposure = "0", eths = eths) %>%
  rbind(collect_age_spec_ir(race_age_spec_irs, exposure = "1", eths = eths))

overall_cases_pys <- bind_rows(
  collect_age_spec_pys_cases(race_age_spec_irs, exposure = "0", eths = eths),
  collect_age_spec_pys_cases(race_age_spec_irs, exposure = "1", eths = eths)
)  %>%
  pivot_wider(
    names_from = ethnicity,
    values_from = c(pys, cases),
    names_sep = "_",
    names_glue = "{ethnicity}_{.value}"
  )

# ----adding all ir results to 1 excel file----
irs_wb <- createWorkbook()

# adding each dataframe as a separate sheet
addWorksheet(irs_wb, "overall age-spec irs")
writeData(irs_wb, sheet = "overall age-spec irs", overall_df)

addWorksheet(irs_wb, "overall cases and pys")
writeData(irs_wb, sheet = "overall cases and pys", overall_cases_pys)

saveWorkbook(
  irs_wb,
  here(
    "scripts",
    "output",
    "tables",
    "overall_age_stratified_irs.xlsx"
  ),
  overwrite = TRUE
)

# ----age-adj irs----

# standard population from 2000 census
std_pop <- c(10805447, 9533545, 8857441, 7415813, 4945367, 4239587)

htn_age_adj <-
  calculate_age_adj_ir(
    data = result$data,
    age_cat_labels = result$age_cat_labels,
    htn = 1,
    std_pop = std_pop,
    per_pys = 1000
  ) %>% do.call(rbind, .) %>%
  mutate(
    exposure = 0,
    ethnicity = c("South Asian", "Chinese", "Japanese", "Filipino", "White")
  )

no_htn_age_adj <-
  calculate_age_adj_ir(
    data = result$data,
    age_cat_labels = result$age_cat_labels,
    htn = 0,
    std_pop = std_pop,
    per_pys = 1000
  ) %>% do.call(rbind, .) %>%
  mutate(
    exposure = 1,
    ethnicity = c("South Asian", "Chinese", "Japanese", "Filipino", "White")
  )

all_age_adj <- rbind(htn_age_adj, no_htn_age_adj) %>%
  #  mutate(across(c(crude_ir, crude_ir_ci_l, crude_ir_ci_u), ~ . / 1000)) %>%
  mutate(across(where(is.numeric), ~ round(., 1)))

write.xlsx(all_age_adj,
           file =   here("scripts",
                         "output",
                         "tables",
                         "age_adjusted_irs.xlsx"))

# ---ir figure---

figure_data <- race_age_spec_irs %>% do.call(rbind, .) %>%
  mutate(
    ethnicity = rep(
      c("South Asian", "Chinese", "Japanese", "Filipino", "White"),
      each = 12
    ),
    exposure = relevel(factor(exposure), ref = "No hypertension"),
    age_range = gsub(" years", "", age_range)
  )


figure_data <- figure_data %>%
  mutate(
    ir = case_when(
      ir == "<5 cases" ~ NA_real_,  
      TRUE ~ as.numeric(ir) 
    ),
    line_nudge = case_when(
      exposure == "No hypertension" ~ -0.25,
      exposure == "Hypertension" ~ 0.25
    ),
    conf.high.ceiling = ifelse(ir_ul < 60, ir_ul, NA),
    arrow_pos = ifelse(ir_ul > 60, 60, NA)
  )


colors <- c("#2E86C1", "#AED6F1")
# colors <- c("#FFDB6D","#00AFBB")

# ---ir figure with vertical bars (horizontal commented out)----
age_spec_ir_plot <-
  ggplot(figure_data,
         aes(
           x = age_range,
           y = ir,
           ymin = 0,
           group = exposure,
           fill = exposure
         )) +
  geom_bar(stat = "identity", position = "dodge", aes(alpha = 0.99)) +
  geom_errorbar(
    aes(ymin = ir_ll, ymax = conf.high.ceiling),
    position = position_nudge(x = figure_data$line_nudge),
    width = 0.15
  ) +
  geom_segment(
    aes(
      x = age_range,
      xend = age_range,
      y = ir_ll,
      yend = arrow_pos
    ),
    position = position_nudge(x = figure_data$line_nudge),
    arrow = arrow(length = unit(0.2, "cm")),
    show.legend = FALSE
  ) +
  geom_text(
    aes(
      x = age_range,
      y = arrow_pos + 1,
      label = round(ir_ul, 2)
    ),
    position = position_nudge(x = figure_data$line_nudge),
    size = 2.75
  ) +
  scale_y_continuous(
    breaks = seq(0, 60, by = 5),
    limits = c(0, 62.5),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = colors) +
  facet_grid(. ~ ethnicity) +
  labs(x = "Age range (years)",
       y = "Incidence rate (per 1,000 PY)",
       title = "Age- and race/ethnicity-specific incidence rates, stratified by hypertension diagnosis at 5+ years prior to baseline",
       fill = "Exposure")  +
  guides(color = "none") +
  scale_alpha(guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

ggsave(
  age_spec_ir_plot,
  file = here("scripts", "output", "figures", "ir_fig.png"),
  width = 15,
  height = 8
)
