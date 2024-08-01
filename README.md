# asian_american_htn_dementia

These scripts are part of an R Project that use the "here" package to input data and export results. The code works in 4 main phases: data setup, exploratory analysis, main analysis, and secondary analysis. Each phase is summarized below. 


Data setup

1. 0_paths.R: Sets path to box and datasets to be used in dataset construction script
2. 1a_dataset_construction.R: Cleans cardiometabolic dataset (wide), derives median blood pressure in year closest to baseline in BP dataset (long), and merges the BP measure with the cardiometabolic data to produce a wide analytic dataset.
3. 1b_multiple_imputation.R: Looks at variable missingness and generates 20 imputed datasets with a maximum of 10 iterations
4. 1c_cleaning_imputed_data.R: Sets up data for table 1, mild object, and separate imputed datasets for 5/7/9 criteria
5. 1d_merging_in_rx_data.R: Uses Rx dataset, derives hypertension medication at 1 year prior to baseline variable to be merged into complete case data for table 1


Exploratory analysis

6. tbl1_functions.R: Functions to label table 1, render continuous variables in table1 object, and create table1 object (can input data/exposure so that it works on 5/7/9-year data)
7. 2a_table1_complete_case.R: Creates 3 different datasets based on membership criteria for each exposure definition (hypertension diagnosis flag 5/7/9 years before baseline), applies variable/value labels, generates table 1s for each of the datasets, calculates ethnicity x hypertension n’s and percentages for header, and exports table 1s in 1 file
8. 2b_table1_imputed.R Using imputed data, creates 3 different datasets based on membership criteria for each exposure definition (hypertension diagnosis flag 5/7/9 years before baseline), applies variable/value labels, generates table 1s for the categorical variables by averaging across imputations, pools continuous variables using Rubin’s rules


Main analysis

9. ir_irr_functions.R: Functions to calculate cases/person-years, calculate age-specific IRs, collect the age-specific IRs into a neater data frame, collect cases/PYs into dataframe (potentially for supplement), calculate age-adjusted IRs
10. 3a_incidence_rates.R: Uses functions from ir_irr_functions script to calculate age-specific IRs for each ethnicity, age-adjusted IRs for each ethnicity; creates a figure for the age-specific IRs
11. 3b_survival_curves.R Creates Kaplan Meier curves using dementia as the outcome and cumulative incidence curves looking at dementia and death using the survminer package and plotting using ggsurvplot.
12. cox_multiple_models_function.R: Function to run multiple cox models (inputting covariates, outcome, and ethnicity) and output tidy results.
13. 3c_cox_models_current_age.R: Makes imputed datasets into long datasets. Uses a custom function to run multiple cox proportional hazards models for each ethnicity/set of covariates, extracts hazard ratios and confidence intervals, calculates age-specific HRs, outputs table 2s and figures.

Secondary analysis

14. 4_cox_cumulative_incidence.R: Written by Juliet Zhou. Creates a function that fits a Cox model, and then generate cumulative incidence curve for one individual (with a set of covariate values) and calculate RR and RD comparing hypertension vs no hypertension.
