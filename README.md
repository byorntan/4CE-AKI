# 4CE COVID-19 Consortium - Acute Kidney Injury

This is a works-in-progress - the script that I have written to process the Phase 2 data from the 4CE consortium into a table containing every AKI episode is still largely incomplete.

## Milestones achieved
1) Detection of AKI using serum Cr trends
   - Using a rolling time window and KDIGO criteria with the minimum serum Cr in the time window as the baseline Cr
   - Using pracma::findpeaks() to detect peaks in Cr and then matching these peaks against the earlier rolling-window-based method to detect true serum Cr peaks
2) Generating graphs of serum Cr trends with ggplot2

## Future Directions
1) Generate a model (starting with GLM-based models) to predict the following:
   - AKI severity
   - Time to AKI onset
   - COVID-19 severity using AKI characteristics

## System Prerequisites
1) R 3.6 or greater (ideally)
2) Following R packages must be installed:
   - dplyr, tidyr, purrr, ggplot2 (all are part of the tidyverse package)
   - data.table
   - zoo
   - pracma
   - RcppRoll

## Additional Files Required
1) ICD9/10 codes of common comorbidities (comorbid_icd_code.csv)
   * This file must be a .csv containing the following headers:
```
icd_version	icd_code	comorbid_type	description
icd_version: ICD version (DIAG-ICD9, DIAG-ICD10)
icd_code: ICD code
comorbid_type: comorbid type (ihd,htn,dm,asthma,copd,bronchiectasis,ild,ckd,pe,dvt,cancer)
description: String of diagnosis (for easier readability)
```

2) ICD9/10 codes of prothrombotic events (thromb_icd_code.csv)
   * This file must be a csv containing following headers:
```
icd_version	icd_code	type	description
icd_version: ICD version (DIAG-ICD9, DIAG-ICD10)
icd_code: ICD code
type: prothrombotic event type (dvt,vt,pe,mi)
description: String of diagnosis (for easier readability)
```

## How to Use
1) Ensure pre-requisite custom files are present alongside the 4CE data tables.
2) Edit the file import script "file_prep.R" to ensure that the correct path to the required CSV files are there
3) Ensure the following scripts are in the same directory:
   - file_prep.R
   - demog_comorbid_thromb_severe.R
   - aki_detection.R
   - medications.R
   - figure1_cr_trends.R
3) Run the main script "run_aki_analysis.R"