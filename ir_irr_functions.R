# IR/IRR functions
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/13/2024

# ----calculating number of person-years contributed by each person----
calculate_pys <- function(data, age_cat, fu_start, fu_end, event_flag) {
  
  # testing
  # data <- wide_data_pre_mi
  # age_cat <- c(59, 65, 70, 75, 80, 85, 120)
  # fu_start <- "survey_age"
  # fu_end = "main_dem_v1_end_age"
  # event_flag <- "main_dem_v1_end_dem_flag"
  
  age_cat_labels <- c()
  
  for (i in 1:(length(age_cat) - 1)) {
    # testing 
    # i <- 1
    
    # initialize the age interval
    cat_start <- age_cat[i]
    cat_end <- age_cat[i+1]
    # create the variable names for py and case contribution
    age_int <- paste0(cat_start, "_", cat_end)
    vars <- paste0(c("py_", "case_"), age_int)
    py_var <- vars[1]
    case_var <- vars[2]
    
    # collect the age intervals for setting up the IR output table
    age_cat_labels <- c(age_cat_labels, age_int)
    
    # calculate py and case contribution by age category
    data <- data %>% 
      mutate(
        !! py_var := case_when(
          # start of fu after age category or end of fu before age category: 
          # contribute 0 pys
          get(fu_start) > cat_end | get(fu_end) <= cat_start ~ 0, 
          
          # start of fu before age category and end of fu during age category: 
          get(fu_start) <= cat_start & get(fu_end) <= cat_end ~ 
            get(fu_end) - cat_start, 
          
          # start of fu before age category and end of fu after age category: 
          get(fu_start) <= cat_start & get(fu_end) > cat_end ~ 
            cat_end - cat_start, 
          
          # start and end of fu during age category:
          get(fu_start) <= cat_end & get(fu_end) <= cat_end ~ 
            get(fu_end) - get(fu_start), 
          
          # start of fu during age category and end of fu after age category:
          get(fu_start) <= cat_end & get(fu_end) > cat_end ~ 
            cat_end - get(fu_start), 
          
          TRUE ~ NA_real_),  
        
        # if end of fu is during age category, case contribution is the dem flag. 
        # otherwise, case contribution is zero
        !! case_var := ifelse(
          get(fu_end) <= cat_end & get(fu_end) > cat_start, 
          get(event_flag), 0
        )
      )
  }
  
  return(list(data = data, age_cat_labels = age_cat_labels))
}

# ----using person-years output to calculate age-specific irs----
calculate_age_specific_ir <- function(data, age_cat_labels, ethnicity, htn_status) {
  data_filtered <- data %>%
    filter(ethnicity_rev %in% ethnicity & htn_dx5yr_flag == htn_status)
  
  age_spec_ir <- tibble(
    age_range = age_cat_labels,
    pys = colSums(select(data_filtered, starts_with("py_"))),
    cases = colSums(select(data_filtered, starts_with("case_")))
  )
  
  # calculate age-specific irs
  age_spec_ir <- age_spec_ir %>%
    mutate(
      ir = ifelse(cases >= 5, cases / pys * 1000, NA),
      ir_se = ifelse(cases >= 5, sqrt(cases) / pys * 1000, NA),
      ir_ll = ifelse(cases >= 5, round(ir - 1.96 * ir_se, 2), NA),
      ir_ul = ifelse(cases >= 5, round(ir + 1.96 * ir_se, 2), NA)
    ) %>%
    mutate(
      ir = case_when(
        cases < 5 ~ "<5 cases",
        TRUE ~ as.character(round(ir, 1))
      )
    ) %>%
    mutate(across(where(is.numeric), ~ round(., 1)))
  
  return(age_spec_ir)
}


# ----collecting age-specific irs into a cleaner dataframe----
collect_age_spec_ir <- function(ir_data, exposure, eths) {
  age_spec_list <- lapply(eths, function(x) {
    if (!is.null(ir_data[[paste0(x, "_htn_", exposure)]])) {
      ir_data[[paste0(x, "_htn_", exposure)]] %>%
        mutate(age_spec_ir = case_when(
          cases < 5 ~ "< 5 events",
          TRUE ~ paste0(ir, " (", ir_ll, ", ", ir_ul, ")")
        )) %>%
        select(age_range, age_spec_ir) %>%
        pivot_wider(names_from = age_range, values_from = age_spec_ir)
    } else {
      NULL
    }
  })
  
  age_spec_df <- do.call(rbind, Filter(Negate(is.null), age_spec_list)) %>%
    mutate(
      ethnicity = eths,
      exposure = ifelse(exposure == "0", "No hypertension", "Hypertension")
    ) %>%
    relocate(exposure, ethnicity)
  
  return(age_spec_df)
}


# ----using irs to calculate irrs comparing htn to no htn----
# calculate_irr <- function(ir_hyp, ir_no_hyp) {
#   combined_ir <- merge(
#     as.data.frame(ir_hyp),
#     as.data.frame(ir_no_hyp),
#     by = "age_range",
#     suffixes = c("_hyp", "_no_hyp")
#   )
#   
#   combined_ir <- combined_ir %>%
#     mutate(
#       irr = case_when(
#         ir_hyp != "<5 cases" & ir_no_hyp != "<5 cases" ~ round(as.numeric(ir_hyp) / as.numeric(ir_no_hyp), 2),
#         TRUE ~ NA_real_
#       ),
#       se = case_when(
#         !is.na(irr) ~ sqrt(1 / as.numeric(ir_hyp) + 1 / as.numeric(ir_no_hyp)),
#         TRUE ~ NA_real_
#       ),
#       ci_lower = case_when(
#         !is.na(irr) ~ round(exp(log(irr) - 1.96 * se), 2),
#         TRUE ~ NA_real_
#       ),
#       ci_upper = case_when(
#         !is.na(irr) ~ round(exp(log(irr) + 1.96 * se), 2),
#         TRUE ~ NA_real_
#       )
#     )
#   
#   return(combined_ir)
# }

# ----collecting person-years and cases into dataframe----
collect_age_spec_pys_cases <- function(ir_data, exposure, eths) {
  age_spec_list <- lapply(eths, function(x) {
    if (!is.null(ir_data[[paste0(x, "_htn_", exposure)]])) {
      ir_data[[paste0(x, "_htn_", exposure)]] %>%
        select(age_range, starts_with("pys"), starts_with("cases")) %>%
        mutate(ethnicity = x)
    } else {
      NULL
    }
  })
  
  age_spec_df <- do.call(rbind, Filter(Negate(is.null), age_spec_list)) %>% 
    mutate(exposure = ifelse(exposure == "0", "No hypertension", "Hypertension"))
  
  return(age_spec_df)
}

# ----age-adj irs----
calculate_age_adj_ir <- function(data, age_cat_labels, std_pop, htn,
                                 per_pys) {
  
  data <- data %>% 
    filter(htn_dx5yr_flag == htn) %>% 
    split(., list(.$ethnicity_rev)) 
  
  out <- lapply(data, function(x) {
    age_spec_ir <- tibble(
      age_range = age_cat_labels, 
      pop_count = std_pop, 
      pop_total = sum(pop_count), 
      pop_prop = pop_count / pop_total, 
      unw_pys = NA, 
      unw_cases = NA
    )
    
    for (i in 1:length(age_cat_labels)) {
      cat <- age_cat_labels[i]
      age_spec_ir[i, "unw_pys"] <- sum(x[, paste0("py_", cat)]) 
      age_spec_ir[i, "unw_cases"] <- sum(x[, paste0("case_", cat)])
    }
    
    age_spec_ir <- age_spec_ir %>% 
      mutate(
        unw_adj_ir = ifelse(unw_cases >= 5, unw_cases / unw_pys, NA),
        unw_adj_ir_se = ifelse(unw_cases >= 5, sqrt(unw_cases / unw_pys^2), NA),
        unw_adj_ir_ci_l = ifelse(unw_cases >= 5, unw_adj_ir - 1.96 * unw_adj_ir_se, NA),
        unw_adj_ir_ci_u = ifelse(unw_cases >= 5, unw_adj_ir + 1.96 * unw_adj_ir_se, NA)
      ) 
    
    adj_ir_out <- age_spec_ir %>% 
      mutate(
        unw_adj_ir_pop = unw_cases / unw_pys * pop_prop,
        unw_adj_ir_pop_var = unw_cases / unw_pys^2 * pop_prop^2
      ) %>% 
      summarise(
        across(c(unw_cases, unw_pys), sum), 
        adj_ir = sum(unw_adj_ir_pop) * per_pys,
        adj_ir_var = sum(unw_adj_ir_pop_var),
        crude_ir = sum(unw_cases) / sum(unw_pys) * per_pys,
        crude_ir_ci_l = (sum(unw_cases) / sum(unw_pys) - 1.96 * sqrt(sum(unw_cases) / sum(unw_pys)^2)) * per_pys,
        crude_ir_ci_u = (sum(unw_cases) / sum(unw_pys) + 1.96 * sqrt(sum(unw_cases) / sum(unw_pys)^2)) * per_pys
      ) %>% 
      mutate(
        adj_ir_se = sqrt(adj_ir_var),
        adj_ir_l = adj_ir - 1.96 * adj_ir_se * per_pys,
        adj_ir_u = adj_ir + 1.96 * adj_ir_se * per_pys
      ) 
    
    out <- adj_ir_out %>% 
      select(unw_cases, unw_pys, adj_ir, adj_ir_l, adj_ir_u, crude_ir, crude_ir_ci_l, crude_ir_ci_u) %>% 
      mutate(across(c(adj_ir, adj_ir_l, adj_ir_u, crude_ir, crude_ir_ci_l, crude_ir_ci_u), 
                    ~ ifelse(!is.na(.x), round(.x, 2), NA)))
  }
  )
  
  return(out)
}
