library(dplyr)
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