# Multiple imputation of wide dataset
# Created by: Natalie Gradwohl
# Adpated from Juliet Zhou's code for Asian Americans ADRD Diabetes project
# Creation date: 2/12/2024

if (!require("pacman"))
  install.packages("pacman", repos = 'http://cran.us.r-project.org')

p_load("here", "tidyverse", "mice", "mitools")

options(scipen = 9999)

# loading data
source(here("scripts", "0_paths.R"))

pre_mi <- readRDS(paste0(path_to_analytic_data, "wide_data_pre_mi.rds"))

# checking missingness----------------------------------------------------------

vars_to_impute  <- c(
  "survey_age",
  "female",
  "ethnicity_rev",
  "usaborn_rev",
  "usabornfather_rev",
  "usabornmother_rev",
  "education_rev",
  "sizeofhh",
  "income",
  "income_pp",
  "maritalstatus",
  "generalhealth", 
  "smoking_status",
  "sr_total_height_in",
  "ehr_ht_median",
  "sbp_closest_to_baseline"
)


missing_summary <- data.frame(varname = vars_to_impute, pctmiss = NA)
row.names(missing_summary) <- vars_to_impute

for (i in vars_to_impute) {
  missing_summary[i, "pctmiss"] <-
    100 * sum(is.na(pre_mi[, i])) / nrow(pre_mi)
  # print(i)
  # print(table(wide_data_pre_mi[, i], exclude = NULL))
}

missing_ordered <- missing_summary[order(missing_summary$pctmiss),]
missing_ordered

ordered_var_list <- missing_ordered$varname

# ordering variables based on least to most missingness
aa_impute <- pre_mi[, ordered_var_list]

# making sure variables are correct type before imputation
str(aa_impute)

factor_cols <- c(
  "female",
  "maritalstatus",
  "usabornmother_rev",
  "usabornfather_rev",
  "ethnicity_rev"
)

ordinal_cols <- c(
  "sizeofhh",
  "education_rev",
  "income",
  "smoking_status",
  "generalhealth"
)

aa_impute <- aa_impute %>%
  mutate_at(factor_cols, as.factor) %>%
  mutate_at(ordinal_cols, as.ordered)

# double checking variable type after converting to factors/ordered factors
str(aa_impute)

# initiating imputation
init = mice(
  aa_impute,
  maxit = 0,
  defaultMethod = c("pmm", "logreg", "polyreg", "polr"),
  seed = 12345
)
meth = init$method
predictM = init$predictorMatrix
init$loggedEvents #NULL
table(init$nmis)

# making sure it doesn't impute these vars
predictM[c("income", "sizeofhh"), "income_pp"] <- 0
predictM

# generating 20 imputed datasets with a maximum of 10 iterations
imputed = mice(
  aa_impute,
  method = meth,
  predictorMatrix = predictM,
  m = 20,
  # RUN 20 ON DESKTOP
  maxit = 10,
  defaultMethod = c("pmm", "logreg", "polyreg", "polr"),
  seed = 12345,
  print = T
)



# saving imputed datasets object
save(imputed, file = paste0(path_to_analytic_data, "imputed_raw.RData"))
