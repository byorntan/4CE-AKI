library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
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

# ==================
# File Preparation
# ==================

# Read in compiled CSV files (these must already contain ALL sites - if not yet done please use rbind to join all site data)
demographics <- read.csv("PatientSummary.csv")
observations <- read.csv("PatientObservations.csv")
course <- read.csv("PatientClinicalCourse.csv")

# Read in file containing ICD9/10 codes of prothrombotic events.
# This file must be a csv containing following headers:
# icd_version icd_code  type  description
#
# icd_version: ICD version (DIAG-ICD9, DIAG-ICD10)
# icd_code: ICD code
# type: prothrombotic event type (dvt,vt,pe,mi)
# description: String of diagnosis (for easier readability)
thromb_ref <- read.csv("thromb_icd_code.csv")
# Split table into ICD9 and ICD10 codes for easier coding later
thromb_icd9_ref <- thromb_ref[thromb_ref$icd_version == "DIAG-ICD9",-c(1,4)]
thromb_icd10_ref <- thromb_ref[thromb_ref$icd_version == "DIAG-ICD10",-c(1,4)]

# Read in file containing mapping of ICD9/10 codes to comorbidities. 
# This file must be a csv containing following headers:
# icd_version icd_code  comorbid_type description
#
# icd_version: ICD version (DIAG-ICD9, DIAG-ICD10)
# icd_code: ICD code
# comorbid_type: comorbid type (ihd,htn,dm,asthma,copd,bronchiectasis,ild,ckd,pe,dvt,cancer)
# description: String of diagnosis (for easier readability)
comorbid_ref <- read.csv("comorbid_icd_code.csv")
# Split the table into ICD-9 and ICD-10 codes for easier comparison and merging later
comorbid_icd9_ref <- comorbid_ref[comorbid_ref$icd_version == "DIAG-ICD9",-c(1,4)]
comorbid_icd10_ref <- comorbid_ref[comorbid_ref$icd_version == "DIAG-ICD10",-c(1,4)]

# first generate a unique ID for each patient
demographics <- demographics %>% mutate(patient_id=paste(site_id,patient_num,sep="_"))
course <- course %>% mutate(patient_id=paste(site_id,patient_num,sep="_"))
observations <- observations %>% mutate(patient_id=paste(site_id,patient_num,sep="_"))

# From this point on, we will be using our custom-generated patient_id as a unique patient identifier
# Reorder the columns in each table to bring patient_id to the first column and remove patient_num
demographics <- demographics[,c(15,1,3:14)]
course <- course[,c(8,1,3:7)]
observations <- observations[,c(7,1,3:6)]

# Generate a diagnosis table
diagnosis <- observations[observations$concept_type %in% c("DIAG-ICD9","DIAG-ICD10"),-6]
colnames(diagnosis) <- c("patient_id","siteid","days_since_admission","concept_type","icd_code")

# Generate a procedures table
procedures <- observations[observations$concept_type %in% c("PROC-ICD9", "PROC-ICD10"),-6]
colnames(procedures) <- c("patient_id","siteid","days_since_admission","icd_version","procedure_code")

# ==============
# Demographics
# ==============
# We want to calculate (1) time to severe event (2) time to death
# We can use Kaplan-Meier curves to determine associated with these severe events
demographics_filt <- demographics %>% mutate(time_to_severe = ifelse(severe == 1, severe_date - admission_date,-999))
demographics_filt <- demographics_filt %>% mutate(time_to_death = ifelse(deceased == 1, death_date - admission_date,-1))
demographics_filt <- demographics_filt %>% mutate(length_stay = ifelse(still_in_hospital==1,days_since_admission,last_discharge_date - admission_date))

# Reorder the columns to be more readable
demographics_filt <- demographics_filt[,c(1,2,11:13,17,8,15,10,16)]
# Final headers for demographics_filt
# patient_id  site_id sex age_group race  length_stay severe  time_to_severe  deceased  time_to_death

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
write.csv(labs_aki_severe,"PatientAKIEvents_Severe.csv",row.names=FALSE)

# ======================================
# Comorbidities & Prothrombotic Events
# ======================================
diag_icd9 <- diagnosis[diagnosis$concept_type == "DIAG-ICD9",]
diag_icd10 <- diagnosis[diagnosis$concept_type == "DIAG-ICD10",]

# Filter comorbids from all diagnoses
comorbid_icd9 <- diag_icd9[diag_icd9$concept_code %in% comorbid_icd9_ref$icd_code,]
comorbid_icd10 <- diag_icd10[diag_icd10$concept_code %in% comorbid_icd10_ref$icd_code,]
# Filter comorbids to be restricted to diagnosis codes made -15days
comorbid_icd9 <- comorbid_icd9[comorbid_icd9$days_since_admission < -15,-c(2:4)]
comorbid_icd10 <- comorbid_icd10[comorbid_icd10$days_since_admission < -15,-c(2:4)]
# Map the comorbid codes
comorbid_icd9 <- merge(comorbid_icd9,comorbid_icd9_ref,by="icd_code",all.x=TRUE)
comorbid_icd10 <- merge(comorbid_icd10,comorbid_icd10_ref,by="icd_code",all.x=TRUE)
comorbid <- rbind(comorbid_icd9,comorbid_icd10)
comorbid <- comorbid[,-2]
comorbid$present <- 1
comorbid <- spread(comorbid,comorbid_type,present)
comorbid[is.na(comorbid)] <- 0

# Filter prothrombotic events from all diagnoses
thromb_icd9 <- diag_icd9[diag_icd9$concept_code %in% thromb_icd9_ref$icd_code,]
thromb_icd10 <- diag_icd10[diag_icd10$concept_code %in% thromb_icd10_ref$icd_code,]
# Filter prothrombotic diagnoses to be restricted to diagnosis codes made after -15days
thromb_icd9 <- thromb_icd9[thromb_icd9$days_since_admission >= -15,-c(2,4)]
thromb_icd10 <- thromb_icd10[thromb_icd10$days_since_admission >= -15,-c(2,4)]
# Map the prothrombotic codes - store day diagnosed
thromb_icd9 <- merge(thromb_icd9,thromb_icd9_ref,by="icd_code",all.x=TRUE)
thromb_icd10 <- merge(thromb_icd10,thromb_icd10_ref,by="icd_code",all.x=TRUE)
thromb_diag <- rbind(thromb_icd9,thromb_icd10)
thromb_diag <- thromb_diag[,c(1,4,2)]
thromb_diag <- spread(thromb_diag,type,days_since_admission)
thromb_diag[is.na(thromb_diag)] <- -999

# =============
# Medications
# =============
medications <- observations[observations$concept_type == "MED-CLASS",]
medications <- medications[,-c(4,6)]
medications <- medications %>% arrange(patient_id,days_since_admission,concept_code)
# Use 15 days as the cutoff for chronic medications
med_chronic <- medications[medications$days_since_admission < -15,]
med_new <- medications[medications$days_since_admission >= -15,]

# Re-code chronic medications into wide format
med_chronic <- med_chronic[!duplicated(med_chronic[,c(1,2,4)]),]
med_chronic <- med_chronic[,-c(2,3)]
med_chronic$concept_code <- paste(med_chronic$concept_code,"old",sep="_")
med_chronic$present <- 1
med_chronic <- spread(med_chronic,concept_code,present)
med_chronic[is.na(med_chronic)] <- 0

# For simplicity of initial analysis, we will use the earliest date where each new medication class is
# used. However, if we are to incorporate a recurrent neural network model to account for temporal changes
# in medications, this approach cannot be used.
med_new <- med_new[!duplicated(med_new[,c(1,2,4)]),]
med_new <- med_new[,c(1,4,3)]
colnames(med_new)[3] <- "start_day"
# Compute the new medications as a time difference from the start of AKI
# If the value is < 0 (and presumably > -365) then the medication was initiated before the start of AKI
# Otherwise it means the medication had been started after the peak Cr had been achieved
# This will give insight into which medications may be potentially nephrotoxic
aki_start_time <- labs_aki_summ[,c(1,12)]
med_new_aki <- merge(med_new,aki_start_time,by="patient_id",all.x=TRUE)
med_new_aki$offset_aki <- med_new_aki$start_day - med_new_aki$aki_start
# Re-code whether medication was given before AKI - 1 = yes, 0 = no
med_new_aki <- mutate(med_before_aki = ifelse(offset_aki <=0,1,0))
med_new_aki <- med_new_aki[,-c(3,4,5)]
med_new_aki <- spread(med_new_aki,concept_code,med_before_aki)
med_new_aki[is.na(med_new_aki)] <- -999

# Generate another table with the start date of the new medications
med_new <- spread(med_new,concept_code,start_day)
med_new[is.na(med_new)] <- -999

# ====================
# Time to Intubation
# ====================
# (1) Determine intubation from procedure code
intubation_code <- c("0BH17EZ","5A093*", "5A094*", "5A095*","96.04","96.70","96.71","96.72")
intubation <- procedures[procedures$procedure_code %in% intubation_code,-c(2,4)]
intubation <- intubation[,c(1,2)]
intubation <- intubation[order(intubation$patient_id,intubation$days_since_admission)]
intubation <- intubation[!duplicated(intubation$patient_id),]
# (2) In some cases intubation may not be coded as a procedure. Hence surrogate way is to determine if 
#     patient had been diagnosed with ARDS and/or VAP
vap_ards_codes <- c("J80","J95.851","518.82","997.31")
vap_ards_diag <- diagnosis[diagnosis$icd_code %in% vap_ards_codes,]
intubation <- rbind(vap_ards_diag[,c(1,3)],intubation)
intubation <- intubation[order(intubation$patient_id,intubation$days_since_admission)]
intubation <- intubation[!duplicated(intubation$patient_id),]


# Compile all the tables together


# Statistical analysis

# =======================================================================================
# Figure 1: Cr trends from start of AKI / after peak Cr for severe vs non-severe groups
# =======================================================================================
# This analysis uses the index AKI episode
# Feel free to modify the code to look at other types of AKI episodes, e.g. you can filter for the most severe AKI episode
# First, filter the labs_aki_severe table to show only the index AKI episodes
aki_index <- labs_aki_severe %>% group_by(patient_id) %>% filter(days_since_admission == min(days_since_admission)) %>% ungroup()
# patient_id site_id days_since_admission  peak_cr  baseline_cr delta_cr  cr_7d cr_90d  aki_kdigo_final akd_7d  akd_90d  aki_start aki_end  time_to_severe  severe_to_aki severe_before_aki
# Generate the patient list including (1) severity indices from this filtered table (2) day of peak Cr
# severe - 0 = never severe, 1 = severe before AKI, 2 = severe after AKI
aki_index <- aki_index %>% group_by(patient_id) %>% mutate(severe = ifelse(time_to_severe != -999,ifelse(severe_before_aki == 1,1,2),0)) %>% ungroup()
aki_index <- aki_index[,c(1,3,17,12)]
colnames(aki_index)[2] <- "peak_cr_time"
# Headers of aki_index: patient_id  peak_cr_time  severe  aki_start

# Now bind this list of days of peak Cr to a new table containing the Cr trends
labs_aki_timeplot <- labs_cr_aki[,c(1,3,4,5,7)]
labs_aki_timeplot <- merge(labs_aki_timeplot,aki_index,by="patient_id",all.x=TRUE)
labs_aki_timeplot <- labs_aki_timeplot[,c(1,7,2,6,3,4,5)]
# patient_id  severe  days_since_admission  peak_cr_time  value min_cr_7day min_cr_7day_retro

# Calculate the day from peak Cr
labs_aki_timeplot <- labs_aki_timeplot %>% group_by(patient_id) %>% mutate(time_from_peak = days_since_admission - peak_cr_time) %>% ungroup()
# Filter this table for Cr values that fall within the 7 day window
labs_aki_timeplot <- labs_aki_timeplot %>% group_by(patient_id) %>% filter(between(time_from_peak,0,7)) %>% ungroup()
# Normalise to baseline values used for AKI calculation
labs_aki_timeplot <- labs_aki_timeplot %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_7day,min_cr_7day_retro)) %>% ungroup()
labs_aki_timeplot <- labs_aki_timeplot %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
labs_aki_timeplot <- labs_aki_timeplot[,c(1,2,8,10)]
# Headers: patient_id  severe  time_from_peak  ratio

# Calculate mean and SD each for severe and non-severe groups
aki_cr_plot_summ <- labs_aki_timeplot %>% group_by(severe,time_from_peak) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()

# Plot the graphs
cr_timeplot <- ggplot(aki_cr_plot_summ,aes(x=time_from_peak,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity")
print(cr_timeplot)

# Plot from start of admission to 10 days
# index_cr <- labs_cr_aki[,c(1,3,4,5,7)]
# index_cr <- merge(index_cr,aki_index,by="patient_id",all.x=TRUE)
# This line removes all patients who do not have AKI
# - if we want to keep them replace it with the commented line below which replaces the null aki_start times with 0
## index_cr[is.na(index_cr)] <- 0
# index_cr <- index_cr[!is.na(index_cr$peak_cr_time),]
# index_cr <- index_cr %>% group_by(patient_id) %>% filter(between(days_since_admission,0,peak_cr_time+7)) %>% ungroup()
# index_cr <- index_cr[order(index_cr$patient_id,index_cr$days_since_admission),]
# index_cr <- index_cr %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_7day,min_cr_7day_retro)) %>% ungroup()
# index_cr <- index_cr %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
# aki_cr_plot_all <- index_cr %>% group_by(severe,days_since_admission) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()
# cr_timeplot_all <- ggplot(aki_cr_plot_all[which(aki_cr_plot_all$days_since_admission < 10),],aes(x=days_since_admission,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity")
# print(cr_timeplot_all)

# Plot from start of AKI to 7 days later 
index_cr <- labs_cr_aki[,c(1,3,4,5,7)]
index_cr <- merge(index_cr,aki_index,by="patient_id",all.x=TRUE)
# This line removes all patients who do not have AKI
index_cr <- index_cr[!is.na(index_cr$aki_start),]
index_cr <- index_cr %>% group_by(patient_id) %>% mutate(time_from_start = days_since_admission - aki_start) %>% ungroup()
index_cr <- index_cr[order(index_cr$patient_id,index_cr$days_since_admission),]
# Change the desired duration as required
index_cr <- index_cr %>% group_by(patient_id) %>% filter(between(time_from_start,0,7)) %>% ungroup()
index_cr <- index_cr %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_7day,min_cr_7day_retro)) %>% ungroup()
index_cr <- index_cr %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
aki_cr_plot_all <- index_cr %>% group_by(severe,time_from_start) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()
cr_timeplot_all <- ggplot(aki_cr_plot_all,aes(x=time_from_start,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity")
print(cr_timeplot_all)
