# mall_agwt_ipm
Code for formatting the data and running the integrated population model for mallard and green-winged teal

# Metadata

## Folders
The "raw_dat" folder contains all the data needed for analysis in their original form.

The subfolder named "environmental_covariates" under the "raw_dat" folder contains raw climatic data downloaded from the National Centers for Environmental Database for the analysis.

The "data" folder contains data that have been formatted and ready for analysis.

## Scripts

"01_mall_agwt_data_formatting.R" is for formatting the individual/population data used in the model (BPop, wing data, and dead-recovery data).

"02_envi_cov_data_formatting.R" is for formatting environmental covariate data.

"03_mall_agwt_run_script.R" contains the code for running the integrated population model and produce figures as results.

"04_mall_agwt_tLTRE.R"is for running the retrospective analysis and producing figures related to the analysis.

"mall_agwt_ipm_clean.jags" contains the model code for the integrated population model.





