library(dplyr)
library(tidyr)
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
med_chronic$concept_code <- paste("old",med_chronic$concept_code,sep="_")
med_chronic$present <- 1
med_chronic <- spread(med_chronic,concept_code,present)
med_chronic[is.na(med_chronic)] <- 0

# Create subtable for ACE-i/ARB pre-exposure
med_acearb_chronic <- med_chronic %>% select(patient_id,old_ACEI,old_ARB)
med_acearb_chronic <- med_acearb_chronic %>% group_by(patient_id) %>% mutate(acei_arb_preexposure = ifelse(old_ACEI + old_ARB > 0,1,0))
med_acearb_chronic <- med_acearb_chronic %>% select(patient_id,acei_arb_preexposure)

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

# Generate simplified table for determining who were started on COAGB near admission
med_coagb_new <- med_new %>% select(patient_id,COAGB)
med_coagb_new$COAGB[med_coagb_new$COAGB < -15] <- 0
med_coagb_new$COAGB[med_coagb_new$COAGB >= -15] <- 1

# Generate simplified table for determining who were started on novel antivirals
med_covid19_new <- med_new %>% select(patient_id,COVIDVIRAL)
# Uncomment following line to select for patients who were started on novel antivirals less than 72h from admission
#med_covid19_new <- med_covid19_new[med_covid19_new$COVIDVIRAL <= 3,]
med_covid19_new <- med_covid19_new %>% group_by(patient_id) %>% mutate(covid_rx = ifelse(COVIDVIRAL >= 0,1,0))
med_covid19_new_date <- med_covid19_new
colnames(med_covid19_new_date)[2] <- "covid_rx_start"
med_covid19_new <- med_covid19_new %>% select(patient_id,covid_rx)