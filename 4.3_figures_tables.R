# prepares figures and tables for the manuscript

# load packages ----

if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("tidyverse", "here", "flextable")

options(scipen = 9999)

# load results and prep formatting ----

load(here("scripts", "output", "bootstrap_results.RData"))

ethn_order <- c("Chinese", "Filipino", "Japanese", "South Asian", "White")
ethn_label <- c("Chinese", "Filipino", "Japanese", "South Asian", "Non-Latino White")

allest <- allest %>% 
  mutate(ethnicity = factor(ethnicity, levels = ethn_order, labels = ethn_label)) %>% 
  arrange(ethnicity)

# HR ----

# this is only for comparison with results obtained without bootstrapping
# since bootstrapping is run only on the last model 

hr <- allest %>% 
  select(ethnicity, matches("^htn_\\d"), starts_with("htn_age_int")) %>%
  pivot_longer(
    cols = -ethnicity, 
    names_to = c("est", ".value"),
    names_sep = "\\."
  ) %>% 
  mutate(
    across(-c(ethnicity, est), exp),
    table_est = paste0(
      sprintf("%0.2f", pe), 
      " (", sprintf("%0.2f", lower), ", ", sprintf("%0.2f", upper), ")"
    )
  )

hr_table <- list()

# main effect (HR at age 75) only
hr_table[[1]] <- hr %>% 
  filter(est != "htn_age_int") %>% 
  mutate(age = parse_number(est)) %>% 
  select(ethnicity, age, table_est) %>% 
  pivot_wider(names_from = ethnicity, values_from = table_est)

# main effect and interaction term 
hr_table[[2]] <- hr %>% filter(est %in% c("htn_75", "htn_age_int")) %>% 
  mutate(est = ifelse(est == "htn_75", "main effect at age 75", "age-interaction")) %>% 
  select(ethnicity, est, table_est) %>% 
  pivot_wider(
    id_cols = ethnicity,
    names_from = est, 
    values_from = table_est
  ) 

# predicted incidence ----

inc <- allest %>% 
  select(ethnicity, starts_with("time")) %>% 
  pivot_longer(
    cols = -ethnicity, 
    names_to = c("est", ".value"),
    names_sep = "\\."
  ) %>% 
  mutate(
    age = str_sub(est, 6, 7) %>% as.numeric(), 
    htn = str_sub(est, 13), 
    htn = ifelse(htn == 0, "No hypertension", "Hypertension") 
  ) 

inc %>% 
  ggplot(aes(x = age)) + 
  geom_line(aes(y = pe, color = htn)) + 
  # geom_ribbon(aes(ymin = lower, ymax = upper, fill = htn), alpha = 0.3) + 
  labs(
    x = "Age (years)", y = "Cumulative incidence (%)", color = NULL, fill = NULL, 
    # subtitle = "Fully adjusted model without age interaction"
  ) + 
  facet_grid(cols = vars(ethnicity)) + 
  theme_bw() +
  theme(
    # strip.background = element_blank(),
    strip.text = element_text(size = 11),
    legend.position = "bottom",
    # legend.position.inside = c(0.85, 0.2)
  )

ggsave(
  here("scripts", "output", "figures", "supp_cumulative_inc_w_ageint.png"),
  height = 4, width = 10
)

# ggsave(
#   here("scripts", "output", "figures", "supp_cumulative_inc_w_ageint_seband.png"),
#   height = 4, width = 10
# )

# RD ----

RD <- allest %>% 
  select(ethnicity, starts_with("RD")) %>% 
  pivot_longer(
    cols = -ethnicity, 
    names_to = c("est", ".value"),
    names_sep = "\\."
  ) %>% 
  mutate(
    age = parse_number(est),
    across(c(pe, lower, upper), function(x) x * 100)
  )

## table formatting ----
RD_table <- RD %>%
  mutate(
    table_est = paste0(
      sprintf("%0.2f", pe), 
      " (", sprintf("%0.2f", lower), ", ", sprintf("%0.2f", upper), ")"
    ), 
    table_est = ifelse(
      (age == 65 & ethnicity != "Non-Latino White") | (age == 70 & ethnicity == "South Asian"), 
      "--", table_est
    )
  ) %>% 
  pivot_wider(
    id_cols = ethnicity, 
    names_from = age, 
    values_from = table_est
  )

## figure ----
RD_upper_limit <- 10
RD_lower_limit <- -10

RD_plot <- RD %>%
  mutate(
    
    # suppress estimates at age 65 for all Asian ethnicities
    # and at 70 for South Asians since these have problematic CI
    # or we can also suppress all estimates at age 65 and age 70 for Asian ethns
    # ER will decide once she is reviewing the manuscript
    across(
      c(pe, lower, upper), 
      function(x) ifelse(
        (age == 65 & ethnicity != "Non-Latino White") | (age == 70 & ethnicity == "South Asian"),
        # age %in% c(65, 70) & ethnicity != "Non-Latino White", 
        NA_real_, x
      )
    ), 
    markpos = ifelse(
      (age == 65 & ethnicity != "Non-Latino White") | (age == 70 & ethnicity == "South Asian"),
      # age %in% c(65, 70) & ethnicity != "Non-Latino White",
      1, NA_real_
    ), 
    conf.high.ceiling = ifelse(upper <= RD_upper_limit, upper, NA),
    conf.low.ceiling = ifelse(lower >= RD_lower_limit, lower, NA),
    upper_arrow_pos = ifelse(upper > RD_upper_limit, RD_upper_limit, NA),
    lower_arrow_pos = ifelse(lower < RD_lower_limit, RD_lower_limit, NA)
  )

RD_plot %>% 
  ggplot(aes(x = age, y = pe)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey30") + 
  geom_point(size = 2) + 
  geom_errorbar(
    aes(ymin = conf.low.ceiling, ymax = conf.high.ceiling), width = 1
  ) +
  geom_segment(
    aes(x = age, xend = age, y = lower_arrow_pos, yend = upper_arrow_pos),
    arrow = arrow(length = unit(0.2, "cm"), ends = "both")
  ) + 
  geom_segment(
    aes(x = age, xend = age, y = conf.low.ceiling, yend = upper_arrow_pos),
    arrow = arrow(length = unit(0.2, "cm"))
  ) + 
  geom_segment(
    aes(x = age, xend = age, y = min(upper, upper_arrow_pos), yend = lower_arrow_pos),
    arrow = arrow(length = unit(0.2, "cm")),
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = age, y = upper_arrow_pos + 0.8, label = round(upper, 1)),
    size = 3
  ) + 
  geom_text(
    aes(x = age, y = lower_arrow_pos - 0.8, label = round(lower, 1)),
    size = 3
  ) + 
  geom_point(
    aes(y = markpos), shape = 4, show.legend = FALSE
  ) + 
  scale_x_continuous(
    breaks = seq(65, 90, 5), minor_breaks = NULL,
    limits = c(63, 92)
  ) + 
  scale_y_continuous(
    breaks = seq(-10, 10, 2), # minor_breaks = NULL,
    limits = c(-11.5, 11.5), expand = c(0.01, 0.01)
  ) + 
  labs(
    x = "Age (years)", y = "Risk difference (%)"
    # subtitle = "Fully adjusted model without age interaction"
  ) + 
  facet_grid(cols = vars(ethnicity)) + 
  theme_bw() +
  theme(
    # strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 11)
  )

ggsave(
  here("scripts", "output", "figures", "supp_RD_w_ageint.png"),
  height = 3, width = 10
)


# RR ----

RR <- allest %>% 
  select(ethnicity, starts_with("RR")) %>% 
  pivot_longer(
    cols = -ethnicity, 
    names_to = c("est", ".value"),
    names_sep = "\\."
  ) %>% 
  mutate(
    age = parse_number(est)
  )

## table formatting ----
RR_table <- RR %>%
  mutate(
    table_est = paste0(
      sprintf("%0.2f", pe), 
      " (", sprintf("%0.2f", lower), ", ", sprintf("%0.2f", upper), ")"
    ), 
    table_est = ifelse(
      (age == 65 & ethnicity != "Non-Latino White") | (age == 70 & ethnicity == "South Asian"), 
      "--", table_est
    )
  ) %>% 
  pivot_wider(
    id_cols = ethnicity, 
    names_from = age, 
    values_from = table_est
  )

## figure ----

RR_upper_limit <- 3.5
RR_lower_limit <- 0

RR_plot <- RR %>%
  mutate(
    across(
      c(pe, lower, upper), 
      function(x) ifelse(
        (age == 65 & ethnicity != "Non-Latino White") | (age == 70 & ethnicity == "South Asian"), 
        NA_real_, x
      )
    ), 
    markpos = ifelse(
      (age == 65 & ethnicity != "Non-Latino White") | (age == 70 & ethnicity == "South Asian"),
      1.2, NA_real_
    ), 
    conf.high.ceiling = ifelse(upper <= RR_upper_limit, upper, NA),
    conf.low.ceiling = ifelse(lower >= RR_lower_limit, lower, NA),
    upper_arrow_pos = ifelse(upper > RR_upper_limit, RR_upper_limit, NA),
    lower_arrow_pos = ifelse(lower < RR_lower_limit, RR_lower_limit, NA)
  )

RR_plot %>% 
  ggplot(aes(x = age, y = pe)) + 
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "grey30") + 
  geom_point(size = 2) + 
  geom_errorbar(
    aes(ymin = conf.low.ceiling, ymax = conf.high.ceiling), width = 1
  ) +
  geom_segment(
    aes(x = age, xend = age, y = conf.low.ceiling, yend = upper_arrow_pos),
    arrow = arrow(length = unit(0.2, "cm"))
  ) +
  geom_text(
    aes(x = age, y = upper_arrow_pos + 0.15, label = round(upper, 1)),
    size = 3
  ) +
  geom_point(
    aes(y = markpos), shape = 4, show.legend = FALSE
  ) + 
  scale_x_continuous(
    breaks = seq(65, 90, 5), minor_breaks = NULL,
    limits = c(63, 92)
  ) + 
  scale_y_continuous(
    breaks = seq(0, 3.5, 0.5), # minor_breaks = seq(0, 3, 0.5),
    limits = c(0, 3.8), expand = c(0.01, 0.01)
  ) +
  labs(
    x = "Age (years)", y = "Risk ratio"
    # subtitle = "Fully adjusted model without age interaction"
  ) + 
  facet_grid(cols = vars(ethnicity)) + 
  theme_bw() +
  theme(
    # strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 11)
  )

ggsave(
  here("scripts", "output", "figures", "supp_RR_w_ageint.png"),
  height = 3, width = 10
)

# save the tables ----

save_as_docx(
  age_spec_HR = hr_table[[1]] %>% flextable() %>% autofit(), 
  hr_main_int = hr_table[[2]] %>% flextable() %>% autofit(), 
  RD_table = RD_table %>% flextable() %>% autofit(), 
  RR_table = RR_table %>% flextable() %>% autofit(), 
  path = here("scripts", "output", "tables", "bootstrap_HR_RD_RR.docx")
)
