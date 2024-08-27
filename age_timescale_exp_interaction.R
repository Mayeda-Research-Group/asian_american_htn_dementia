


# path_to_box <- "~/Library/CloudStorage/Box-Box/"
# path_to_imputed <-
#   "Asian_Americans_dementia_data/htnetes_adrd/Imputed_data/"
# path_to_proj <-
#   "Asian_Americans_dementia/Manuscripts/AA_ADRD_htnetes/Code/Cleaned_Scripts/"
#
# imputed_tte_data <-
#   readRDS(file = paste0(
#     path_to_box,
#     path_to_imputed,
#     "aa_adrd_htnetes_imputed_tte_data.rds"
#   ))
#
# # prepare for analysis ----
# # run the script that preps the data according to R1 revision
# source(paste0(path_to_box, path_to_proj, "2.0.R1_data_prep.R"))
#
# # use one imputed dataset as an example ----


# set up paths and datasets ----------------------------------------------------
data <- readRDS(here("datasets", "wide_data_pre_mi_5.rds"))

m0 <- coxph(Surv(survey_age, end_age, event_dem) ~
              htn_dx5yr_flag * asian + female + edu_3,
            data = data)

# make data long
cuts <- c(min(data$survey_age), # the minimum start age
          seq(60, 102, 1), # all the intervals in the middle
          max(data$end_age)) # the maximum end of followup age )

longdata <- survSplit(
  data    = data,
  formula = Surv(survey_age, end_age, event_dem) ~ .,
  cut     = cuts,
  # rename the variables that indicate the start and end of each interval
  start   = "age_start",
  end     = "age_end"
)


# # take a look at the data to understand the structure
# longdata %>%
#   select(subjid,
#          htn_dx5yr_flag, asian, female,
#          age_start, age_end, main_dem_v1_end_dem_flag) %>%
#   View()
#
# # the original survey age and event age are preserved with decimal places
# longdata %>%
#   group_by(subjid) %>%
#   slice_head() %>%
#   # slice_tail() %>%
#   View()

# fit the same model on the long dataset
m0_long <- coxph(Surv(age_start, age_end, event_dem) ~
                   htn_dx5yr_flag * asian + female + edu_3,
                 data = longdata)

# compare the two models: one using the wide dataset, one using the long dataset
# they are equivalent
m0
m0_long


# models with interaction between current age and exposure ----

longdata <- longdata %>%
  mutate(current_age = age_start,
         current_agec75 = age_start - 75)

m1 <- coxph(
  Surv(age_start, age_end, event_dem) ~
    htn_dx5yr_flag * asian +
    female + edu_3 + htn_dx5yr_flag:current_agec75,
  data = longdata
)

m1

# do we also want the asian:htn association to also vary by age?
m2 <- coxph(
  Surv(age_start, age_end, main_dem_v1_end_dem_flag) ~
    htn_dx5yr_flag * asian +
    female + edu_3 +
    htn_dx5yr_flag:current_agec75 + htn_dx5yr_flag:asian:current_agec75,
  data = longdata
)

m2

# extract coefficients relevant to htn
m1coef <-
  m1$coefficients[str_detect(names(m1$coefficients), "htn")]
m2coef <-
  m2$coefficients[str_detect(names(m2$coefficients), "htn")]

# calculate and plot HR over current age
hr_tib <- expand_grid(htn = 1,
                      asian = c(0, 1),
                      current_age = 60:90) %>%
  mutate(
    current_agec75 = current_age - 75,
    htn_x_asian = htn * asian,
    htn_x_current_agec75 = htn * current_agec75,
    htn_x_current_agec75_x_asian = htn * asian * current_agec75
  ) %>%
  relocate(asian, current_age, current_agec75, .after = everything())

hr_tib$hr_m1 <- exp(as.matrix(hr_tib[, 1:3]) %*% m1coef)
hr_tib$hr_m2 <- exp(as.matrix(hr_tib[, 1:4]) %*% m2coef)

# I don't think this is the best way to do this -- very likely to multiply
# the wrong terms!
# I think the later example where I use model.matrix to achieve this is better

hr_tib %>%
  pivot_longer(
    cols = starts_with("hr_"),
    names_to = "model",
    names_prefix = "hr_",
    values_to = "hr"
  ) %>%
  mutate(asian = factor(asian)) %>%
  ggplot(aes(x = current_age, y = hr, group = asian)) +
  geom_line(aes(color = asian)) +
  # geom_hline(aes(yintercept = 1.53, linetype = "paper result: white"),
  #            color = "grey50") +
  # geom_hline(aes(yintercept = 1.43, linetype = "paper result: asian"),
  #            color = "grey50") +
  # geom_vline(aes(xintercept = 83), color = "grey50") +
  facet_wrap(~ model)

# compare with plr ----
library(splines)
library(gtools)

# produce table with HR/OR estimates side by side
compare <- function(coxm, plrm) {
  coxm %>%
    tidy(exponentiate = TRUE) %>%
    select(term, estimate, p.value) %>%
    rename(cox_estimate = estimate, cox_p.value = p.value) %>%
    left_join(
      plrm %>%
        tidy(exponentiate = TRUE) %>%
        select(term, estimate, p.value) %>%
        filter(str_detect(term, "Intercept|ns", negate = TRUE)) %>%
        rename(plr_estimate = estimate, plr_p.value = p.value),
      by = c("term")
    ) %>%
    mutate(
      across(ends_with("estimate"), function(x)
        round(x, 2)),
      across(ends_with("p.value"), stars.pval),
      cox_estimate = paste0(cox_estimate, cox_p.value),
      plr_estimate = paste0(plr_estimate, plr_p.value)
    ) %>%
    select(term, cox_estimate, plr_estimate)
}


plrm0 <- longdata %>%
  glm(
    main_dem_v1_end_dem_flag ~ ns(age_start, 5) +
      htn_dx5yr_flag * asian + female + edu_3,
    data = .,
    family = binomial()
  )

compare(m0, plrm0)

plrm1 <- longdata %>%
  glm(
    main_dem_v1_end_dem_flag ~ ns(age_start, 5) +
      htn_dx5yr_flag * asian + female +
      htn_dx5yr_flag:current_agec75,
    data = .,
    family = binomial()
  )

compare(m1, plrm1)

plrm2 <- longdata %>%
  glm(
    main_dem_v1_end_dem_flag ~ ns(age_start, 5) +
      htn_dx5yr_flag * asian + female +
      htn_dx5yr_flag:current_agec75 +
      htn_dx5yr_flag:asian:current_agec75,
    data = .,
    family = binomial()
  )

compare(m2, plrm2)

longdata <- longdata %>% mutate(ethnicity_rev = factor(
  ethnicity_rev,
  levels = c(9, 1, 2, 3, 5),
  labels = c("White", "South Asian", "Chinese", "Japanese", "Filipino")
))

# include specific asian ethnicities ----
m1_ethns <- coxph(
  Surv(age_start, age_end, main_dem_v1_end_dem_flag) ~
    htn_dx5yr_flag * ethnicity_rev +
    female + htn_dx5yr_flag:current_agec75,
  data = longdata
)

m2_ethns <- coxph(
  Surv(age_start, age_end, main_dem_v1_end_dem_flag) ~
    htn_dx5yr_flag * ethnicity_rev +
    female +
    htn_dx5yr_flag:current_agec75 + htn_dx5yr_flag:ethnicity_rev:current_agec75,
  data = longdata
)

m1_ethns
m2_ethns


plrm1_ethns <- longdata %>%
  glm(
    main_dem_v1_end_dem_flag ~ ns(age_start, 5) +
      htn_dx5yr_flag * ethnicity_rev + female + edu_3 +
      htn_dx5yr_flag:current_agec75,
    data = .,
    family = binomial()
  )

compare(m1_ethns, plrm1_ethns)

# extract coefficients relevant to htn
m1ethncoef <-
  m1_ethns$coefficients[str_detect(names(m1_ethns$coefficients), "htn")]
m2ethncoef <-
  m2_ethns$coefficients[str_detect(names(m2_ethns$coefficients), "htn")]

# calculate and plot HR over current age
hr_tib <- expand_grid(
  htn_dx5yr_flag = 1,
  ethnicity_rev = c(9, 1, 2, 3, 5) ,
  current_age = 60:90
) %>%
  mutate(
    current_agec75 = current_age - 75,
    ethnicity_rev = factor(
      ethnicity_rev,
      levels = c(9, 1, 2, 3, 5),
      labels = c("White", "South Asian", "Chinese", "Japanese", "Filipino")
    )
  )

hr_tib$hr_m1 <- exp(
  model.matrix(
    ~ htn_dx5yr_flag + htn_dx5yr_flag:ethnicity_rev +
      htn_dx5yr_flag:current_agec75,
    data = hr_tib
  )[, names(m1ethncoef)] %*% m1ethncoef
)

hr_tib$hr_m2 <- exp(
  model.matrix(
    ~ htn_dx5yr_flag + htn_dx5yr_flag:ethnicity_rev +
      htn_dx5yr_flag:current_agec75 +
      htn_dx5yr_flag:ethnicity_rev:current_agec75,
    data = hr_tib
  )[, names(m2ethncoef)] %*% m2ethncoef
)

hr_tib %>%
  pivot_longer(
    cols = starts_with("hr_"),
    names_to = "model",
    names_prefix = "hr_",
    values_to = "hr"
  ) %>%
  ggplot(aes(x = current_age, y = hr, group = ethnicity_rev)) +
  geom_line(aes(color = ethnicity_rev)) +
  facet_wrap(~ model) +
  facet_grid(rows = vars(model), cols = vars(ethnicity_rev))


# # another way is to use tt() to incorporate current age
# # but can't get it to run, with the following error message
# # Error: vector memory exhausted (limit reached?)
# # I suspect that the time intervals are being broken into pieces that are
# # too small and computation is not possible.
# m3 <- coxph(
#   Surv(SURVEY_AGE, MAIN_DEM_V1_END_AGE, MAIN_DEM_V1_END_DEM_FLAG) ~
#     htn_DX5YR_FLAG + tt(htn_DX5YR_FLAG) +
#     FEMALE,
#   data = data,
#   tt = function(x, t, ...) x * t
# )
