# cell_quota_model_of_prostate_cancer
This is the core code for running the prostate cancer models in my PhD dissertation.
FIRST: the folder "Core fitting and forecasting" contains the codes that can be used for chapters 2 and 3.
"trial10" gives the fitting sample for patient 15, while "forecasting_62_trial10" gives the forecasting comparisons.
The "Core fitting and forecasting" folder only contains one fitting sample (patient 15); however, the best fit parameters for all patients are also included for both models in folders "parameters_newTwo_pop" and "parameters_two_pop". Thus, the "forecasting_62_trial10" can be used to produce all forecasting comparisons.
SECOND: the folder "Core identifiability" contains the code that is used in chapter 4.
Running "identifiability_fixed_mu" produces one profile likelihood plot with respect to parameter mu. 
The code can be extended similarly for all relevant parameters.
