library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(RcppRoll)
library(pracma)
library(zoo)

aki_kdigo_grade <- function(x) {
  creat = as.numeric(x[4])
  baseline_7d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  grade = 0
  diff = creat - baseline_48h
  ratio = creat/baseline_7d
  if(diff > 26.5 || ratio >= 1.5) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 353.7 || ratio >= 3) {
    grade = 3
  }
  grade
}

aki_kdigo_grade_retro <- function(x) {
  # Instead of using the pure KDIGO definition, this looks at future Cr values, and 
  # determines whether the current value, in retrospect, represents an episode of AKI
  creat = as.numeric(x[4])
  baseline_7d = as.numeric(x[7])
  baseline_48h = as.numeric(x[8])
  grade = 0
  diff = creat - baseline_48h
  ratio = creat/baseline_7d
  if(diff > 26.5 || ratio >= 1.5) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 353.7 || ratio >= 3) {
    grade = 3
  }
  grade
}

akd_grade_7d <- function(x) {
  creat = as.numeric(x[4])
  baseline_7d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  baseline_7d_retro = as.numeric(x[7])
  baseline_48h_retro = as.numeric(x[8])
  baseline = min(baseline_7d,baseline_48h,baseline_7d_retro,baseline_48h_retro)
  cr_7d = as.numeric(x[9])
  grade = 0
  ratio = creat/baseline
  # We will code grade B/C as 0.5
  if(ratio > 1.25) {
    grade = 0.5
  }
  if(ratio >= 1.5 & ratio < 2) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 353.7 || ratio >= 3) {
    grade = 3
  }
  grade
}

akd_grade_90d <- function(x) {
  creat = as.numeric(x[4])
  baseline_7d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  baseline_7d_retro = as.numeric(x[7])
  baseline_48h_retro = as.numeric(x[8])
  baseline = min(baseline_7d,baseline_48h,baseline_7d_retro,baseline_48h_retro)
  cr_90d = as.numeric(x[10])
  grade = 0
  ratio = creat/baseline
  # We will code grade B/C as 0.5
  if(ratio > 1.25) {
    grade = 0.5
  }
  if(ratio >= 1.5 & ratio < 2) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 353.7 || ratio >= 3) {
    grade = 3
  }
  grade
}

# ====================
# AKI Detection Code
# ====================
# For the purposes of generating a table of AKI events, we will have to create a new data.table
# with the serum creatinine levels.
# We will then use this table, labs_cr_aki, to generate a summary table containing details of 
# each AKI event and the post-AKI recovery labs

# Extract serum Cr levels
labs_cr_aki <- observations[observations$concept_code == '2160-0',] #LOINC code for Cr 2160-0
# Remove unnecessary columns
labs_cr_aki <- labs_cr_aki[,-c(4,5)]

# There are two possible scenarios which we have to consider when detecting each AKI event:
# (1) AKI occurs after admission
#   - easy to detect with the formal KDIGO definition as we only need to use older data points
#   - we can use an approach similar to what is used in the MIMIC-III data-set: use a rolling
#     time frame and detect the lowest Cr value to use as our baseline
# (2) Patient presents with an AKI at the time of admission 
#   - this makes it harder for us to determine what is the true baseline especially with limited
#     longitudinal data
#   - Hence one way to get around this is to generate a "retrospective" baseline (i.e. look into
#     future Cr values and find the minimum in a rolling timeframe) and use this as a surrogate
#     baseline
#   - Such a workaround is the only feasible way of dealing with missing retrospective data though
#     it is likely that we may miss quite a number of true AKIs using this method

# Scenario (1): AKI occurs during admission
# Find minimum Cr level in a rolling 7 day timeframe
labs_cr_aki <- data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
labs_cr_aki <- labs_cr_aki[CJ(unique(patient_id),seq(min(days_since_admission)-7,max(days_since_admission)))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
labs_cr_aki <- labs_cr_aki %>% group_by(patient_id) %>% mutate(min_cr_7day = roll_min(value,8,fill=NA,na.rm=TRUE,partial=TRUE,align="right")) %>% filter(!is.na(value)) %>% ungroup()
# Find minimum Cr level in a rolling 2 day timeframe (48h)
labs_cr_aki <- data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
labs_cr_aki <- labs_cr_aki[CJ(unique(patient_id),seq(min(days_since_admission)-2,max(days_since_admission)))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
labs_cr_aki <- labs_cr_aki %>% group_by(patient_id) %>% mutate(min_cr_48h = roll_min(value,3,fill=NA,na.rm=TRUE,partial=TRUE,align="right")) %>% filter(!is.na(value)) %>% ungroup()

# Scenario (2): Patient presents with an AKI already on board
# Find minimum Cr level in a rolling 7 day timeframe
labs_cr_aki <- data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
labs_cr_aki <- labs_cr_aki[CJ(unique(patient_id),seq(min(days_since_admission),max(days_since_admission)+7))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
labs_cr_aki <- labs_cr_aki %>% group_by(patient_id) %>% mutate(min_cr_7day_retro = roll_min(value,8,fill=NA,na.rm=TRUE,partial=TRUE,align="left")) %>% filter(!is.na(value)) %>% ungroup()
# Find minimum Cr level in a rolling 2 day timeframe (48h)
labs_cr_aki <- data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
labs_cr_aki <- labs_cr_aki[CJ(unique(patient_id),seq(min(days_since_admission),max(days_since_admission)+2))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
labs_cr_aki <- labs_cr_aki %>% group_by(patient_id) %>% mutate(min_cr_48h_retro = roll_min(value,3,fill=NA,na.rm=TRUE,partial=TRUE,align="left")) %>% filter(!is.na(value)) %>% ungroup()

# Another outcome we are interested in is to look at acute kidney disease, AKD (in between AKI and CKD)
# We will use the definitions proposed for AKD as described by Chawla et. al. 2017 (ref (1))
# We are interested in renal recovery at the 7-day and 90-day timepoint
# We will use a cutoff of recovery to 1.25x baseline Cr as recovery, as used in ref (2)
# References:
# 1. Chawla, L., Bellomo, R., Bihorac, A. et al. Acute kidney disease and renal recovery: consensus report
#    of the Acute Disease Quality Initiative (ADQI) 16 Workgroup. Nat Rev Nephrol 13, 241-257 (2017). 
#    https://doi.org/10.1038/nrneph.2017.2
# 2. Pannu, N., James, M., Hemmelgarn, B. & Klarenbach, S. Association between AKI, Recovery of Renal 
#    Function, and Long-Term Outcomes after Hospital Discharge. Clinical Journal of the American Society of
#    Nephrology 8, 194-202 (2013). https://doi.org/10.2215/CJN.06480612

# Generate sCr levels at +7d (cr_7d) and +90d (cr_90d) timepoints (for determining post-AKI recovery, AKD)
labs_cr_aki <- setDT(labs_cr_aki)[,':='(cr_7d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+7,incbounds = TRUE)],1),cr_90d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][between(labs_cr_aki$days_since_admission[labs_cr_aki_rv$patient_id==patient_id],days_since_admission,days_since_admission+90,incbounds = TRUE)],1)),by=c('patient_id','days_since_admission')][]

# At this point, our table has these headers:
# patient_id  site_id  days_since_admission  value min_cr_7day min_cr_48h  min_cr_7day_retro min_cr_48h_retro  cr_7d cr_90d

# Now we have to start grading AKI severity at each time point
# This approach is similar to how the MIMIC-III dataset generates AKI severity
# Generate two columns using both the formal KDIGO AKI definition and the modified retrospective AKI definition
labs_cr_aki$aki_kdigo <- apply(labs_cr_aki,1,aki_kdigo_grade)
labs_cr_aki$aki_kdigo_retro <- apply(labs_cr_aki,1,aki_kdigo_grade_retro)
labs_cr_aki <- labs_cr_aki %>% group_by(patient_id,days_since_admission) %>% mutate(aki_kdigo_final = max(aki_kdigo,aki_kdigo_retro)) %>% ungroup()

# Generate two columns grading AKD severity at 7d and 90d (grade 0B/C is coded as 0.5)
labs_cr_aki$akd_7d <- apply(labs_cr_aki,1,akd_grade_7d)
labs_cr_aki$akd_90d <- apply(labs_cr_aki,1,akd_grade_90d)

# Generate delta Cr (current Cr - baseline Cr<7day>)
labs_cr_aki <- labs_cr_aki %>% group_by(patient_id) %>% mutate(delta_cr = max(value-min_cr_7day, value-min_cr_7day_retro))%>% ungroup()

# Generate a separate peak table using the pracma::findpeaks function
# We will use these as a filter for detecting the real peaks in Cr and match these against those detected using the rolling window method
# This will also give us the start and end times of the peaks
labs_cr_peak <- data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
labs_cr_peak <- labs_cr_peak %>% group_by(patient_id) %>% complete(days_since_admission = full_seq(days_since_admission,1)) %>% mutate(value = na.approx(value,na.rm=FALSE))
labs_cr_peak <- labs_cr_peak %>% split(.$patient_id) %>% map(~findpeaks(.$value)) %>% map_df(~data_frame(peak_cr = .x[,1],index_peak=.x[,2],index_start = .x[,3],index_end = .x[,4]),.id='patient_id') %>% arrange(patient_id,index_peak)
# The problem with the above method is that it does not interpolate leading NAs - we need to offset the indexes obtained using pracma::findpeaks()
labs_peak_offset <- labs_cr_aki %>% group_by(patient_id) %>% summarise(day_start = min(days_since_admission))
labs_peak_offset$day_offset <- labs_peak_offset$day_start - 1
labs_peak_offset <- labs_peak_offset[,c(1,3)]
labs_cr_peak <- merge(labs_cr_peak,labs_peak_offset,by="patient_id",all.x=TRUE)
labs_cr_peak[3:5] <- labs_cr_peak[3:5]+labs_cr_peak[,6]
colnames(labs_cr_peak)[3] <- "days_since_admission"
# Now we have the peak Cr, start and end times of the peaks, and the time of the maximum Cr all stored in labs_cr_peak

# At this point, our table has these headers:
# patient_id  site_id  days_since_admission  value min_cr_7day min_cr_48h  min_cr_7day_retro min_cr_48h_retro  cr_7d cr_90d aki_kdigo aki_kdigo_retro aki_kdigo_final akd_7d  akd_90d delta_cr

# Now to generate a table with the AKI peak data (absolute Cr values for peak, baseline, increment, Cr in 7 days and 
# 90 days post-peak, and start and end of each AKI episode)
labs_aki_summ <- labs_cr_aki
labs_aki_summ$baseline_cr <- labs_aki_summ$value - labs_aki_summ$delta_cr
# Clean up the data to have fewer columns
labs_aki_summ <- labs_aki_summ[,-c(5:8,11,12)]
labs_aki_summ <- labs_aki_summ[,c(1:3,11,4,10,5:9)]
colnames(labs_aki_summ)[5] <- "peak_cr"
# Merge in peak start and end data and filter those with valid peak data 
# - these are more likely to reflect true AKIs
labs_aki_summ <- merge(labs_aki_summ,labs_cr_peak,by=c("patient_id","peak_cr","days_since_admission"),all.x=TRUE)
# patient_id  peak_cr days_since_admission  site_id baseline_cr  delta_cr cr_7d cr_90d aki_kdigo_final akd_7d  akd_90d  index_start index_end
labs_aki_summ <- labs_aki_summ[,c(1,4,3,2,5:13)]
labs_aki_summ <- labs_aki_summ[!is.na(labs_aki_summ$index_start),]
labs_aki_summ <- labs_aki_summ[order(labs_aki_summ$patient_id,labs_aki_summ$days_since_admission),]
# Apply filter to remove events which do not fulfill KDIGO criteria
labs_aki_summ <- labs_aki_summ[labs_aki_summ$aki_kdigo_final >= 1,]
colnames(labs_aki_summ)[12] <- "aki_start"
colnames(labs_aki_summ)[13] <- "aki_end"
# Final headers for labs_aki_summ:
# patient_id site_id days_since_admission  peak_cr  baseline_cr delta_cr  cr_7d cr_90d  aki_kdigo_final akd_7d  akd_90d  aki_start aki_end
# days_since_admission - time at which peak Cr is achieved
# aki_start - time at which Cr begins to rise
# aki_end - time at which Cr goes back to baseline / before the start of next AKI

# Save the generated AKI table for future reference
write.csv(labs_aki_summ,"PatientAKIEvents.csv",row.names=FALSE)

# We also want to generate tables to determine (1) if AKIs occurred before/after severe disease onset 
# (2) how long before/after disease severity
# These tables will help in segregating the populations for analysis later
severe_time <- demographics_filt[,c(1,8)]
labs_aki_severe <- merge(labs_aki_summ,severe_time,by="patient_id",all.x=TRUE)
labs_aki_severe <- labs_aki_severe %>% group_by(patient_id) %>% mutate(severe_to_aki = ifelse(time_to_severe != -999, time_to_severe - aki_start,-999)) %>% ungroup()
labs_aki_severe <- labs_aki_severe %>% group_by(patient_id) %>% mutate(severe_before_aki = ifelse(severe_to_aki < 0,1,0)) %>% ungroup()
# Final headers for labs_aki_severe:
# patient_id site_id days_since_admission  peak_cr  baseline_cr delta_cr  cr_7d cr_90d  aki_kdigo_final akd_7d  akd_90d  aki_start aki_end  time_to_severe  severe_to_aki severe_before_aki

# Save the generated AKI tables for future reference
write.csv(labs_aki_summ,"PatientAKIEvents.csv",row.names=FALSE)
write.csv(labs_aki_severe,"PatientAKIEvents_Severe.csv",row.names=FALSE)