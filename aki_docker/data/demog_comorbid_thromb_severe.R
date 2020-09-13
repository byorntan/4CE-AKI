library(dplyr)
library(tidyr)
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
comorbid <- comorbid %>% distinct()
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
thromb_diag <- thromb_diag %>% distinct()
thromb_diag <- spread(thromb_diag,type,days_since_admission)
thromb_diag[is.na(thromb_diag)] <- -999

# Final headers for comorbid table
# Note: order of columns may depend on the overall characteristics of your patient population
# patient_id	ihd,htn,dm,asthma,copd,bronchiectasis,ild,ckd,pe,dvt,cancer
# Values stored are in binary, 1 = present, 0 = absent

# Final headers for thromb_diag table
# Note: order of columns may depend on the overall characteristics of your patient population
# patient_id	dvt,vt,pe,mi
# Values stored are the number of days since admission when the event was first recorded
# (stored as -999 if no such event occurred)

# ====================
# Time to Intubation
# ====================
# (1) Determine intubation from procedure code
intubation_code <- c("0BH13EZ","0BH17EZ","0BH18EZ","0B21XEZ","5A09357","5A09358","5A09359","5A0935B","5A0935Z","5A09457","5A09458","5A09459","5A0945B","5A0945Z","5A09557","5A09558","5A09559","5A0955B","5A0955Z","96.7","96.04","96.70","96.71","96.72")
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

# Final table headers:
# patient_id	days_since_admission
# days_since_admission = time to first intubation event

# ==================
# Time to RRT
# ==================
rrt_code <- c("5A1D70Z","5A1D80Z","5A1D90Z","3E1.M39Z","549.8","399.5")
#hd_code <- c("5A1D70Z","5A1D80Z","5A1D90Z","399.5")
#pd_code <- c("3E1.M39Z","549.8")
rrt <- procedures[procedures$procedure_code %in% rrt_code,-c(2,4)]
rrt <- rrt[,c(1,2)]
rrt <- rrt[order(rrt$patient_id,rrt$days_since_admission)]
rrt <- rrt[!duplicated(rrt$patient_id),]

# Generate list of patients already on RRT prior to admission
# This list can be used to exclude ESRF patients in subsequent analyses
patients_already_rrt <- rrt$patient_id[rrt$days_since_admission < 0]

# Generate list of patients who were only initiated on RRT for the first time ever
# during admission for COVID-19
rrt_new <- rrt[!(rrt$patient_id %in% patients_already_rrt),]

# Final table headers for rrt and rrt_new:
# patient_id	days_since_admission
# days_since_admission = time to first RRT event