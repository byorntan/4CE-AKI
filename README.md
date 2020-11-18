# 4CE COVID-19 Consortium - Acute Kidney Injury

This is a works-in-progress - the script that I have written to process the Phase 2 data from the 4CE consortium into a table containing every AKI episode is still largely incomplete.

## Milestones achieved
1) Detection of AKI using serum Cr trends
   - Using a rolling time window and KDIGO criteria with the minimum serum Cr in the time window as the baseline Cr
2) Generating graphs of serum Cr trends with ggplot2
   - Plotting serum Cr trends of patients with/without AKI and severity
   - Plotting serum Cr trends of patients with/without remdesivir/lopinavir+ritonavir use and severity
3) Basic time-to-event analyses for (a) time to recovery to 1.25x baseline serum Cr, (b) time to death

## Future Directions
1) Generate a model (starting with GLM-based models) to predict the following:
   - AKI severity
   - Time to AKI onset
   - COVID-19 severity using AKI characteristics

## System Prerequisites
1) R 3.6 or greater (ideally - most of the code was written on R 4.0 although it has been known to work on 3.5.3)
2) Following R packages must be installed:
   - dplyr, tidyr, purrr, ggplot2 (all are part of the tidyverse package)
   - data.table
   - zoo
   - RcppRoll
   - table1
   - survival
   - survminer
I have also included a modified version of the original Dockerfile script in https://github.com/covidclinical/Phase2.0_Docker_Analysis that will generate a Docker image containing these additional packages for the AKI analyses. These files are under a separate folder 'aki_docker'.

## How to Use
Most of the code is now being updated directly in the covidclinical repository. For the latest version, please visit https://github.com/covidclinical/Phase2.1AKIRPackage
Note that your data must conform to the 4CE Consortium format and must be able to pass QC standards.
This repository will be occasionally updated with more stable versions of the code.
