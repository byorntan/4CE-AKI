library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(RcppRoll)
library(zoo)

aki_kdigo_grade <- function(x) {
  creat = as.numeric(x[4])
  baseline_7d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  grade = 0
  diff = creat - baseline_48h
  ratio = round(creat/baseline_7d,2)
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
  ratio = round(creat/baseline_7d,2)
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
  ratio = round(creat/baseline,2)
  diff = creat - baseline
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
  ratio = round(creat/baseline,2)
  diff = creat - baseline
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

pos_min <- function(cr,day,lag=TRUE,gap=7) {
  len = length(cr)
  day_pos = day
  for(i in 1:len) {
    if(lag) {
      j = max(1,i-gap)
      pos = j-1+which.min(cr[j:i])
      day_pos[i] = day[pos]
    } else {
      j = min(i+gap,len)
      pos = i-1+which.min(cr[i:j])
      day_pos[i] = day[pos]
    }
  }
  day_pos
}

which.peaks <- function(x,partial=TRUE,decreasing=FALSE) {
  if(decreasing) {
    if(partial) {
      which(diff(c(FALSE,diff(x) > 0,TRUE)) > 0)
    } else {
      which(diff(diff(x)>0)>0) + 1
    }
  } else {
    if(partial) {
      which(diff(c(TRUE,diff(x) >= 0,FALSE)) < 0)
    } else {
      which(diff(diff(x)>=0)<0) + 1
    }
  }
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


# Now we are going to generate the start days of each AKI
labs_cr_aki_tmp <- labs_cr_aki
labs_cr_aki_tmp$valid = 1

# Find the day of the minimum Cr used for grading AKIs (taken as baseline)
labs_cr_aki_tmp <- labs_cr_aki_tmp %>% group_by(patient_id) %>% complete(days_since_admission = full_seq(days_since_admission,1)) %>% mutate(value = na.fill(value,Inf))
labs_cr_aki_tmp2 <- labs_cr_aki_tmp
labs_cr_aki_tmp3 <- labs_cr_aki_tmp
labs_cr_aki_tmp2 <- labs_cr_aki_tmp2 %>% split(.$patient_id) %>% map(~pos_min(.$value,.$days_since_admission)) %>% map_df(~data_frame(.x),.id='patient_id')
colnames(labs_cr_aki_tmp2)[2] <- "day_min"
labs_cr_aki_tmp3 <- labs_cr_aki_tmp3 %>% split(.$patient_id) %>% map(~pos_min(.$value,.$days_since_admission,lag=FALSE)) %>% map_df(~data_frame(.x),.id='patient_id')
colnames(labs_cr_aki_tmp3)[2] <- "day_min_retro"
labs_cr_aki_tmp4 <- cbind(labs_cr_aki_tmp,"day_min" = labs_cr_aki_tmp2$day_min,"day_min_retro" = labs_cr_aki_tmp3$day_min_retro)
labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[!is.na(labs_cr_aki_tmp4$valid),]

# Generate delta_cr
labs_cr_aki_tmp4 <- labs_cr_aki_tmp4 %>% group_by(patient_id,days_since_admission) %>% mutate(min_cr_7d_final = min(min_cr_7day,min_cr_7day_retro)) %>% mutate(delta_cr = value - min_cr_7d_final)

# Use the largest delta_cr to find the peak of each AKI
labs_cr_aki_delta_maxima <- labs_cr_aki_tmp4 %>% group_by(patient_id) %>% summarise(days_since_admission=days_since_admission[which.peaks(delta_cr,decreasing=FALSE)],delta_maxima = delta_cr[which.peaks(delta_cr,decreasing=FALSE)])
labs_cr_aki_delta_maxima$delta_is_max = 1
labs_cr_aki_tmp4 <- merge(labs_cr_aki_tmp4,labs_cr_aki_delta_maxima,by=c("patient_id","days_since_admission"),all.x=TRUE)

# Filter for KDIGO grades > 0
labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[labs_cr_aki_tmp5$aki_kdigo_final > 0,]
labs_cr_aki_tmp4[is.na(labs_cr_aki_tmp4)] <- 0
# Filter for maxima of delta_cr (which should give us the peaks)
labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$delta_is_max > 0,]

# Filter and reorder columns to generate our final table of all AKI events
labs_aki_summ <- labs_cr_aki_tmp4 %>% select(patient_id,site_id,days_since_admission,value,day_min,day_min_retro,min_cr_7day,min_cr_48h,min_cr_7day_retro,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)

# Remove our temporary tables (comment these out to check output)
rm(labs_cr_aki_tmp)
rm(labs_cr_aki_tmp2)
rm(labs_cr_aki_tmp3)
rm(labs_cr_aki_tmp4)

# Final headers for labs_aki_summ:
# patient_id,site_id,days_since_admission,value,day_min,day_min_retro,min_cr_7day,min_cr_48h,min_cr_7day_retro,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d
# days_since_admission - time at which peak Cr is achieved
# day_min - time at which Cr begins to rise

# We also want to generate tables to determine (1) if AKIs occurred before/after severe disease onset 
# (2) how long before/after disease severity
# These tables will help in segregating the populations for analysis later
severe_time <- demographics_filt[,c(1,8)]
labs_aki_severe <- merge(labs_aki_summ,severe_time,by="patient_id",all.x=TRUE)
labs_aki_severe <- labs_aki_severe %>% group_by(patient_id) %>% mutate(severe_to_aki = ifelse(time_to_severe != -999, time_to_severe - aki_start,-999)) %>% ungroup()
labs_aki_severe <- labs_aki_severe %>% group_by(patient_id) %>% mutate(severe_before_aki = ifelse(severe_to_aki < 0,1,0)) %>% ungroup()
# Final headers for labs_aki_severe:
# patient_id,site_id,days_since_admission,value,day_min,day_min_retro,min_cr_7day,min_cr_48h,min_cr_7day_retro,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d  time_to_severe  severe_to_aki severe_before_aki

# Save the generated AKI tables for future reference
write.csv(labs_aki_summ,"PatientAKIEvents.csv",row.names=FALSE)
write.csv(labs_aki_severe,"PatientAKIEvents_Severe.csv",row.names=FALSE)