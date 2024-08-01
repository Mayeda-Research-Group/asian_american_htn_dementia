# FUNCTION TO BE PULLED INTO SCRIPT 3e
# Created by: Natalie Gradwohl
# Creation date: 2/13/2024


# ----function to run Cox models with multiple covariates----
multiple_models <- function(data, predictors, ethnicity, response) {
  # data <- imputed_tte_data_long
  # predictors <- "htn_dx5yr_flag + current_age_c75:htn_dx5yr_flag"
  # ethnicity <- 9
  # response <- eval(age_surv_long)
  
  m <- data %>% length()
  sep_models <- vector(mode = "list", length = m)
  
  for (i in 1:m) {
    data_i <- data[[i]]
    
    sep_models[[i]] <-
      coxph(
        as.formula(paste(response, "~", predictors)),
        subset = (ethnicity_rev == ethnicity),
        data = data_i
      )
  }
  
  coefs <- lapply(sep_models, function(x) {
    x$coefficient
  })
  vcov <- lapply(sep_models, function(x) {
    vcov(x)
  })
  
  coef_summary <- MIcombine(coefs, vcov)
  
  #age calculation
  
  # pool the models to get a combined summary
  linear_summary <-
    pool(sep_models) %>% broom::tidy(conf.int = TRUE, exp = TRUE)
  
  hr_summary <- linear_summary %>%
    mutate_if(is.numeric, round, 2) %>%
    mutate(
      estimate_table = paste0(estimate, " (", conf.low, ", " , conf.high, ")"),
      ethnicity = ethnicity
    ) %>%
    select(term, estimate_table, ethnicity)
  
  
  output_list <-
    list(coef_summary = coef_summary, hr_summary = hr_summary)
  return(output_list)
  
}

# ----age-specific HR calculation----
# using coefficients and variance-covariance matrix from fully-adjusted model
age_hr <- function(age, coef, vcov) {
  # age <- 80
  c <- c(1, 0, 0, 0, 0, 0, age-75) # subtracting 75 from age because htn*age coefficient is centered at age 75
  pe <- coef %*% c # log HR
  
  var <- t(c) %*% vcov %*% c
  se <- sqrt(var) # log HR
  
  dat <- data.frame(hr = exp(pe), 
                    lower_ci = exp(pe - 1.96 * se),  
                    upper_ci = exp(pe + 1.96 * se)) 
  
  return(dat)
  
}

                      