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