# Table 1 Functions
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/13/2024

# ----variable/value labelling function----
tbl1_label <- function (data) {
  tb1_data <- data %>%
    mutate(
      female = factor(
        female,
        levels = c(1, 0),
        labels = c("Yes", "No")
      ),
      
      asian = factor(
        asian,
        levels = c(0, 1),
        labels = c("White", "Asian")
      ),
      
      ethnicity_rev = factor(
        ethnicity_rev,
        levels = c(2, 5, 3, 1, 9),
        labels = c("Chinese", "Filipino", "Japanese",
                   "South Asian", "White")
        
      ),
      
      usaborn_rev = factor(
        usaborn_rev,
        levels = c(1, 0),
        labels = c("Born in the US", "Foreign born")
      ),
      
      edu_3 = case_when(
        education_rev %in% c(1, 2) ~ "< HS",
        education_rev %in% c(3, 4) ~ "HS degree or some college",
        education_rev %in% c(5, 6) ~ "College degree"
      ) %>% factor(levels = c("< HS", "HS degree or some college", "College degree")),
      
      maritalstatus = case_when(
        maritalstatus %in% c(1, 3, 4) ~ "Not married",
        maritalstatus == 2 ~ "Married",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("Married", "Not married")),
      
      generalhealth_3 = case_when(
        generalhealth %in% c(1,2) ~ "Excellent or very good", 
        generalhealth %in% c(3) ~ "Good",
        generalhealth %in% c(4,5) ~ "Fair or poor"
      ) %>% factor(levels = c("Excellent or very good", "Good", "Fair or poor")),
      
      diab_dx5yr_flag = factor(
        diab_dx5yr_flag,
        levels = c(1, 0),
        labels = c("Yes", "No")
      ),
      
      diab_dx7yr_flag = factor(
        diab_dx7yr_flag,
        levels = c(1, 0),
        labels = c("Yes", "No")
      ),
      
      diab_dx9yr_flag = factor(
        diab_dx9yr_flag,
        levels = c(1, 0),
        labels = c("Yes", "No")
      ),
      
      htn_dx5yr_flag = factor(
        htn_dx5yr_flag,
        levels = c(1, 0),
        labels = c("Yes", "No")
      ),
      
      htn_dx7yr_flag = factor(
        htn_dx7yr_flag,
        levels = c(1, 0),
        labels = c("Yes", "No")
      ),
      
      htn_dx9yr_flag = factor(
        htn_dx9yr_flag,
        levels = c(1, 0),
        labels = c("Yes", "No")
      ),
      
      first_prevstroke_flag = factor(
        first_prevstroke_flag,
        levels = c(1, 0),
        labels = c("Yes", "No")
      ),
      
      main_dem_v1_end_type = factor(
        main_dem_v1_end_type,
        levels = c(
          "DEMENTIA",
          "DEATH",
          "ADMIN CENSORED",
          "END OF MEMBERSHIP",
          "CENSORED 90+"
        ),
        labels = c(
          "Dementia",
          "Death",
          "Administratively censored",
          "End of membership",
          "Censored 90+"
        )
      ),
      
      smoking_status = factor(
        smoking_status,
        levels = 1:3,
        labels = c("Never", "Former", "Current")
      ),
      
      eversmoke = case_when(
        smoking_status == 1 ~ "No", 
        smoking_status %in% c(2, 3) ~ "Yes"
      ) %>% factor(levels = c("Yes", "No")),
  
      sr_depress = factor(
        sr_depress,
        levels = c(1, 0),
        labels = c("Yes", "No/Missing")
      ),
      
      rx_1yr_pre_surv = factor(
        rx_1yr_pre_surv,
        levels = c(1, 0),
        labels = c("Yes", "No")
      )
      
    ) %>%
    
    # releveling factors
    
    mutate(
     # htn_dx5yr_flag = relevel(htn_dx5yr_flag, ref = "No"),
      female = relevel(female, ref = "No"),
      edu_3 = relevel(edu_3, ref = "< HS"),
      usaborn_rev = relevel(usaborn_rev, ref = "Born in the US")
    )
  
  
  var_label(tb1_data) <- list(
    survey_age = "Survey age, years [mean (SD)]",
    female = "Female, n (%)",
    education_rev = "Education attainment, n (%)",
    edu_ge_college = "College degree or more, n (%)",
    edu_3 = "Education attainment, n (%)",
    usaborn_rev = "US-born, n (%)",
    maritalstatus = "Marital status, n (%)",
    income_pp = "Household-adjusted income, dollars [mean (SD)]",
    generalhealth_3 = "General health, n (%)", 
    main_dem_v1_end_type = "End of follow up event, n (%)",
    main_dem_v1_fu_time = "Follow up time, years [mean (SD)]",
    diab_dx5yr_flag = "Diabetes exposure 5+ years pre-survey",
    diab_dx7yr_flag = "Diabetes exposure 7+ years pre-survey",
    diab_dx9yr_flag = "Diabetes exposure 9+ years pre-survey",
    htn_dx5yr_flag = "Hypertension diagnosis 5+ years pre-survey, n (%)",
    htn_dx7yr_flag = "Hypertension diagnosis 7+ years pre-survey, n (%)",
    htn_dx9yr_flag = "Hypertension diagnosis 9+ years pre-survey, n (%)",
    first_prevstroke_flag = "History of stroke, n (%)",
    ehr_ht_median = "EHR height, in [mean (SD)]",
    smoking_status = "Smoking status, n (%)",
    eversmoke = "Ever smoked, n (%)", 
    sbp_closest_to_baseline = "SBP in the year closest to baseline [mean (SD)]",
    sr_bmi = "BMI [mean (SD)]",
    sr_depress = "Self-rated depression, n (%)",
    rx_1yr_pre_surv = "Taking anti-hypertensive 1 year prior to survey, n (%)"
  )
  
  return(tb1_data)
  
}

# ----custom function to render continuous variables----
my.render.cont <- function(x) {
  with(
    stats.apply.rounding(
      stats.default(x),
      digits = 1,
      rounding.fn = round_pad
    ),
    c("", "Mean (SD)" = sprintf("%s (%s)", MEAN, SD))
  )
}

# customized function for rendering categorical variables in table 1's post-MI
my.render.cat_cc <- function(x, value.prefix = "") {
  N <- length(x)
  freqs <- table(x)
  pcts <- round(freqs / N * 100, 1)
  freqs_formatted <-
    paste0(round(freqs, 0), " (", as.numeric(pcts), "%)")
  names(freqs_formatted) <- paste0(value.prefix, names(freqs))
  
  out <- c("", freqs_formatted)
  return(out)
}

my.render.cat_imp <- function(x, value.prefix = "") {
  N <- length(x) / 20
  freqs <- table(x) / 20
  pcts <- round(freqs / N * 100, 1)
  freqs_formatted <-
    paste0(round(freqs, 0), " (", as.numeric(pcts), "%)")
  names(freqs_formatted) <- paste0(value.prefix, names(freqs))
  
  out <- c("", freqs_formatted)
  return(out)
}

# function to apply Rubin's rule to imputed continuous variables

imp_cont_var_tbl1 <- function(df, cont_vars, grp) {
  # df <- tbl1_data_post_MI_5
  # cont_vars is a list that looks like:
  # cont_vars <- list("INCOME_PP" = 0, "EHR_HT_MEDIAN" = 1)
  # where the names of the items are variable names, and
  # the numbers are rounding digits
  # grp <- c("ETHNICITY_REV", "DIAB_DX5YR_FLAG")
  
  out <- df %>%
    group_by(across(all_of(grp))) %>%
    group_keys()
  
  for (var in names(cont_vars)) {
    temp <- df %>%
      group_by(across(all_of(c("imp", grp)))) %>%
      summarise(mean_imp = mean(get(var), na.rm = T),
                var_imp = var(get(var), na.rm = T)) %>%
      ungroup() %>%
      group_by(across(all_of(grp))) %>%
      summarise(mean = mean(mean_imp),
                sd = mean(var_imp) %>% sqrt()) %>%
      mutate(out = paste0(
        round_pad(mean, cont_vars[[var]]),
        " (",
        round_pad(sd, cont_vars[[var]]),
        ")"
      ))
    out[var] <- temp$out
  }
  
  return(t(out))
}


# ----creating actual table1 object----
tbl1_render <- function (data,
                         varlist,
                         exposure,
                         rowlabelhead,
                         cat_fx,
                         caption) {
  tbl1_df <-
    data %>%
    table1(
      as.formula(paste(varlist,
                       exposure)),
      data = .,
      overall = NULL,
      rowlabelhead = rowlabelhead,
      render.continuous = my.render.cont,
      render.categorical = cat_fx,
      caption = caption
    )
  
  return(tbl1_df)
  
}
