
#' Runs the analytic workflow for the AKI project
#'
#' @keywords 4CE
#' @export

runAnalysis <- function(is_obfuscated=TRUE,obfuscation_value=3) {

    ## make sure this instance has the latest version of the quality control and data wrangling code available
    devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)

    ## get the site identifier associated with the files stored in the /4ceData/Input directory that 
    ## is mounted to the container
    currSiteId = FourCePhase2.1Data::getSiteId()

    ## run the quality control
    FourCePhase2.1Data::runQC(currSiteId)

    ## DO NOT CHANGE ANYTHING ABOVE THIS LINE

    ## To Do: implement analytic workflow, saving results to a site-specific 
    ## file to be sent to the coordinating site later via submitAnalysis()

    ## ========================================
    ## PART 1: Read in Data Tables
    ## ========================================
    
    demographics <- read.csv("Input/LocalPatientSummary.csv")
    observations <- read.csv("Input/LocalPatientObservations.csv")
    data(thromb_ref)
    data(comorbid_ref)
    data(thromb_icd9_ref)
    data(comorbid_icd9_ref)
    data(thromb_icd10_ref)
    data(comorbid_icd10_ref)
    
    # first generate a unique ID for each patient
    demographics <- demographics %>% dplyr::mutate(patient_id=paste(currSiteId,patient_num,sep="_"))
    # course <- course %>% mutate(patient_id=paste(currSiteId,patient_num,sep="_"))
    observations <- observations %>% dplyr::mutate(patient_id=paste(currSiteId,patient_num,sep="_"))
    
    # From this point on, we will be using our custom-generated patient_id as a unique patient identifier
    # Reorder the columns in each table to bring patient_id to the first column and remove patient_num
    demographics <- demographics %>% dplyr::select(patient_id,siteid,admission_date,days_since_admission,last_discharge_date,still_in_hospital,severe_date,severe,death_date,deceased,sex,age_group,race,race_collected)
    # course <- course[,c(8,1,3:7)]
    observations <- observations %>% dplyr::select(patient_id,siteid,days_since_admission,concept_type,concept_code,value)
    
    # Generate a diagnosis table
    diagnosis <- observations[observations$concept_type %in% c("DIAG-ICD9","DIAG-ICD10"),-6]
    colnames(diagnosis) <- c("patient_id","siteid","days_since_admission","concept_type","icd_code")
    
    # Generate a procedures table
    procedures <- observations[observations$concept_type %in% c("PROC-ICD9", "PROC-ICD10"),-6]
    colnames(procedures) <- c("patient_id","siteid","days_since_admission","icd_version","procedure_code")
    
    demographics_filt <- demographics %>% dplyr::mutate(time_to_severe = ifelse(severe == 1, as.numeric(as.Date(severe_date) - as.Date(admission_date)),NA))
    demographics_filt <- demographics_filt %>% dplyr::mutate(time_to_death = ifelse(deceased == 1, as.numeric(as.Date(death_date) - as.Date(admission_date)),NA))
    demographics_filt <- demographics_filt %>% dplyr::mutate(length_stay = ifelse(still_in_hospital==1,days_since_admission,as.numeric(as.Date(last_discharge_date) - as.Date(admission_date))))
    
    # Reorder the columns to be more readable
    demographics_filt <- demographics_filt %>% dplyr::select(patient_id,siteid,sex,age_group,race,length_stay,severe,time_to_severe,deceased,time_to_death)
    
    # Final headers for demographics_filt
    # patient_id  siteid sex age_group race  length_stay severe  time_to_severe  deceased  time_to_death
    
    # Comorbidities & Prothrombotic Events
    # ======================================
    diag_icd9 <- diagnosis[diagnosis$concept_type == "DIAG-ICD9",]
    diag_icd10 <- diagnosis[diagnosis$concept_type == "DIAG-ICD10",]
    # 
    # Filter comorbids from all diagnoses
    comorbid_icd9 <- diag_icd9[diag_icd9$icd_code %in% comorbid_icd9_ref$icd_code,]
    comorbid_icd10 <- diag_icd10[diag_icd10$icd_code %in% comorbid_icd10_ref$icd_code,]
    # Filter comorbids to be restricted to diagnosis codes made -15days
    comorbid_icd9 <- comorbid_icd9[comorbid_icd9$days_since_admission < -15,-c(2:4)]
    comorbid_icd10 <- comorbid_icd10[comorbid_icd10$days_since_admission < -15,-c(2:4)]
    # Map the comorbid codes
    comorbid_icd9 <- merge(comorbid_icd9,comorbid_icd9_ref,by="icd_code",all.x=TRUE)
    comorbid_icd10 <- merge(comorbid_icd10,comorbid_icd10_ref,by="icd_code",all.x=TRUE)
    comorbid <- rbind(comorbid_icd9,comorbid_icd10)
    comorbid <- comorbid %>% dplyr::select(patient_id,comorbid_type)
    comorbid$present <- 1
    comorbid <- comorbid %>% dplyr::distinct()
    comorbid <- comorbid %>% tidyr::spread(comorbid_type,present)
    comorbid[is.na(comorbid)] <- 0
    
    # 
    # # Filter prothrombotic events from all diagnoses
    # thromb_icd9 <- diag_icd9[diag_icd9$icd_code %in% thromb_icd9_ref$icd_code,]
    # thromb_icd10 <- diag_icd10[diag_icd10$icd_code %in% thromb_icd10_ref$icd_code,]
    # # Filter prothrombotic diagnoses to be restricted to diagnosis codes made after -15days
    # thromb_icd9 <- thromb_icd9[thromb_icd9$days_since_admission >= -15,-c(2,4)]
    # thromb_icd10 <- thromb_icd10[thromb_icd10$days_since_admission >= -15,-c(2,4)]
    # # Map the prothrombotic codes - store day diagnosed
    # thromb_icd9 <- merge(thromb_icd9,thromb_icd9_ref,by="icd_code",all.x=TRUE)
    # thromb_icd10 <- merge(thromb_icd10,thromb_icd10_ref,by="icd_code",all.x=TRUE)
    # thromb_diag <- rbind(thromb_icd9,thromb_icd10)
    # thromb_diag <- thromb_diag %>% dplyr::select(patient_id,type,days_since_admission) %>% dplyr::distinct()
    # thromb_diag <- thromb_diag %>% tidyr::spread(type,days_since_admission)
    # thromb_diag[is.na(thromb_diag)] <- NA
    
    # Final headers for comorbid table
    # Note: order of columns may depend on the overall characteristics of your patient population
    # patient_id	ihd,htn,dm,asthma,copd,bronchiectasis,ild,ckd,pe,dvt,cancer
    # Values stored are in binary, 1 = present, 0 = absent
    
    # Final headers for thromb_diag table
    # Note: order of columns may depend on the overall characteristics of your patient population
    # patient_id	dvt,vt,pe,mi
    # Values stored are the number of days since admission when the event was first recorded
    # (stored as NA if no such event occurred)
    
    # Time to Intubation
    # ====================
    # (1) Determine intubation from procedure code
    intubation_code <- c("0BH13EZ","0BH17EZ","0BH18EZ","0B21XEZ","5A09357","5A09358","5A09359","5A0935B","5A0935Z","5A09457","5A09458","5A09459","5A0945B","5A0945Z","5A09557","5A09558","5A09559","5A0955B","5A0955Z","96.7","96.04","96.70","96.71","96.72")
    intubation <- procedures[procedures$procedure_code %in% intubation_code,-c(2,4)]
    intubation <- intubation[,c(1,2)]
    intubation <- intubation[order(intubation$patient_id,intubation$days_since_admission),]
    intubation <- intubation[!duplicated(intubation$patient_id),]
    # (2) In some cases intubation may not be coded as a procedure. Hence surrogate way is to determine if 
    #     patient had been diagnosed with ARDS and/or VAP
    vap_ards_codes <- c("J80","J95.851","518.82","997.31")
    vap_ards_diag <- diagnosis[diagnosis$icd_code %in% vap_ards_codes,]
    intubation <- rbind(vap_ards_diag[,c(1,3)],intubation)
    intubation <- intubation[order(intubation$patient_id,intubation$days_since_admission),]
    intubation <- intubation[!duplicated(intubation$patient_id),]
    
    # Final table headers:
    # patient_id	days_since_admission
    # days_since_admission = time to first intubation event
    
    # Time to RRT
    # ==================
    rrt_code <- c("5A1D70Z","5A1D80Z","5A1D90Z","3E1.M39Z","549.8","399.5")
    #hd_code <- c("5A1D70Z","5A1D80Z","5A1D90Z","399.5")
    #pd_code <- c("3E1.M39Z","549.8")
    rrt <- procedures[procedures$procedure_code %in% rrt_code,-c(2,4)]
    rrt <- rrt[,c(1,2)]
    rrt <- rrt[order(rrt$patient_id,rrt$days_since_admission),]
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
    
    # ====================
    # PART 2: AKI Detection Code
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
    # Find minimum Cr level in a rolling 90 day timeframe
    
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    #labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission)-90,max(days_since_admission)))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(labs_cr_aki$patient_id),seq(min(labs_cr_aki$days_since_admission)-90,max(labs_cr_aki$days_since_admission)))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_90d = RcppRoll::roll_min(value,91,fill=NA,na.rm=TRUE,partial=TRUE,align="right")) %>% dplyr::filter(!is.na(value))
    # Find minimum Cr level in a rolling 2 day timeframe (48h)
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    #labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission)-2,max(days_since_admission)))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(labs_cr_aki$patient_id),seq(min(labs_cr_aki$days_since_admission)-2,max(labs_cr_aki$days_since_admission)))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_48h = RcppRoll::roll_min(value,3,fill=NA,na.rm=TRUE,partial=TRUE,align="right")) %>% dplyr::filter(!is.na(value))
    
    # Scenario (2): Patient presents with an AKI already on board
    # Find minimum Cr level in a rolling 7 day timeframe
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    #labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission),max(days_since_admission)+7))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(labs_cr_aki$patient_id),seq(min(labs_cr_aki$days_since_admission),max(labs_cr_aki$days_since_admission)+7))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_retro_7day = RcppRoll::roll_min(value,8,fill=NA,na.rm=TRUE,partial=TRUE,align="left")) %>% dplyr::filter(!is.na(value))
    # Find minimum Cr level in a rolling 2 day timeframe (48h)
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    #labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission),max(days_since_admission)+2))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(labs_cr_aki$patient_id),seq(min(labs_cr_aki$days_since_admission),max(labs_cr_aki$days_since_admission)+2))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_48h_retro = RcppRoll::roll_min(value,3,fill=NA,na.rm=TRUE,partial=TRUE,align="left")) %>% dplyr::filter(!is.na(value))
    
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
    labs_cr_aki <- data.table::setDT(labs_cr_aki)[,':='(cr_7d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+7,incbounds = TRUE)],1),cr_90d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+90,incbounds = TRUE)],1)),by=c('patient_id','days_since_admission')][]
    
    # At this point, our table has these headers:
    # patient_id  siteid  days_since_admission  value min_cr_90d min_cr_48h  min_cr_retro_7day min_cr_48h_retro  cr_7d cr_90d
    
    # Now we have to start grading AKI severity at each time point
    # This approach is similar to how the MIMIC-III dataset generates AKI severity
    # Generate two columns using both the formal KDIGO AKI definition and the modified retrospective AKI definition
    labs_cr_aki$aki_kdigo <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::aki_kdigo_grade)
    labs_cr_aki$aki_kdigo_retro <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::aki_kdigo_grade_retro)
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(aki_kdigo_final = max(aki_kdigo,aki_kdigo_retro))
    
    # Generate two columns grading AKD severity at 7d and 90d (grade 0B/C is coded as 0.5)
    labs_cr_aki$akd_7d <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::akd_grade_7d)
    labs_cr_aki$akd_90d <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::akd_grade_90d)
    
    # Now we are going to generate the start days of each AKI
    labs_cr_aki_tmp <- labs_cr_aki
    labs_cr_aki_tmp$valid = 1
    
    # Find the day of the minimum Cr used for grading AKIs (taken as baseline)
    labs_cr_aki_tmp <- labs_cr_aki_tmp %>% dplyr::group_by(patient_id) %>% tidyr::complete(days_since_admission = tidyr::full_seq(days_since_admission,1)) %>% dplyr::mutate(value = zoo::na.fill(value,Inf))
    labs_cr_aki_tmp2 <- labs_cr_aki_tmp
    labs_cr_aki_tmp3 <- labs_cr_aki_tmp
    labs_cr_aki_tmp2 <- labs_cr_aki_tmp2 %>% split(.$patient_id) %>% purrr::map(~pos_min(.$value,.$days_since_admission)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_id')
    colnames(labs_cr_aki_tmp2)[2] <- "day_min"
    labs_cr_aki_tmp3 <- labs_cr_aki_tmp3 %>% split(.$patient_id) %>% purrr::map(~pos_min(.$value,.$days_since_admission,lag=FALSE)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_id')
    colnames(labs_cr_aki_tmp3)[2] <- "day_min_retro"
    labs_cr_aki_tmp4 <- cbind(labs_cr_aki_tmp,"day_min" = labs_cr_aki_tmp2$day_min,"day_min_retro" = labs_cr_aki_tmp3$day_min_retro)
    labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[!is.na(labs_cr_aki_tmp4$valid),]
    
    # Generate delta_cr
    labs_cr_aki_tmp4 <- labs_cr_aki_tmp4 %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(min_cr_7d_final = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::mutate(delta_cr = value - min_cr_7d_final)
    
    # Use the largest delta_cr to find the peak of each AKI
    labs_cr_aki_delta_maxima <- labs_cr_aki_tmp4 %>% dplyr::group_by(patient_id) %>% dplyr::filter(delta_cr %in% delta_cr[which.peaks(delta_cr,decreasing=FALSE)])
    labs_cr_aki_delta_maxima$delta_is_max = 1
    labs_cr_aki_delta_maxima <- labs_cr_aki_delta_maxima %>% dplyr::rename(delta_maxima = delta_cr) %>% dplyr::select(patient_id,days_since_admission,delta_maxima,delta_is_max)
    labs_cr_aki_tmp4 <- merge(labs_cr_aki_tmp4,labs_cr_aki_delta_maxima,by=c("patient_id","days_since_admission"),all.x=TRUE)
    
    # Generate a separate table (for reference) of all creatinine peaks not fulfilling KDIGO AKI criteria
    labs_cr_nonaki <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$aki_kdigo_final == 0,]
    labs_cr_nonaki[is.na(labs_cr_nonaki)] <- 0
    labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$delta_is_max > 0,]
    labs_cr_nonaki <- labs_cr_nonaki %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)
    
    # Filter for KDIGO grades > 0
    labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$aki_kdigo_final > 0,]
    labs_cr_aki_tmp4[is.na(labs_cr_aki_tmp4)] <- 0
    # Filter for maxima of delta_cr (which should give us the peaks)
    labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$delta_is_max > 0,]
    
    # Filter and reorder columns to generate our final table of all AKI events
    labs_aki_summ <- labs_cr_aki_tmp4 %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)
    
    labs_aki_summ <- labs_aki_summ %>% dplyr::distinct(patient_id,days_since_admission,.keep_all=TRUE)
    
    # Final headers for labs_aki_summ:
    # patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d
    # days_since_admission - time at which peak Cr is achieved
    # day_min - time at which Cr begins to rise
    
    # Generate the highest Cr peak for non-AKI peaks detected
    labs_nonaki_summ <- labs_cr_nonaki %>% dplyr::group_by(patient_id) %>% dplyr::slice(which.max(delta_cr))
    
    # We also want to generate tables to determine (1) if AKIs occurred before/after severe disease onset 
    # (2) how long before/after disease severity
    # These tables will help in segregating the populations for analysis later
    severe_time <- demographics_filt %>% dplyr::select(patient_id,severe,time_to_severe)
    labs_aki_severe <- merge(labs_aki_summ,severe_time,by="patient_id",all.x=TRUE)
    labs_aki_severe <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_to_aki = ifelse(!is.na(time_to_severe), time_to_severe - day_min,NA))
    labs_aki_severe <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_before_aki = ifelse(severe_to_aki < 0,1,0))
    # Final headers for labs_aki_severe:
    # patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d  severe time_to_severe  severe_to_aki severe_before_aki
    
    labs_nonaki_severe <- merge(labs_nonaki_summ,severe_time,by="patient_id",all.x=TRUE)
    labs_nonaki_severe$severe_to_aki <- NA 
    labs_nonaki_severe$severe_before_aki <- 1 
    
    ## Save the generated AKI tables for future reference / debugging (note: these will NOT be uploaded!!)
    #write.csv(labs_aki_summ,"PatientAKIEvents.csv",row.names=FALSE)
    #write.csv(labs_aki_severe,"PatientAKIEvents_Severe.csv",row.names=FALSE)
    
    # =============
    # Medications
    # =============
    medications <- observations[observations$concept_type == "MED-CLASS",]
    medications <- medications[,-c(4,6)]
    medications <- medications %>% dplyr::arrange(patient_id,days_since_admission,concept_code)
    # Use 15 days as the cutoff for chronic medications
    med_chronic <- medications[medications$days_since_admission < -15,]
    med_new <- medications[medications$days_since_admission >= -15,]
    
    # Re-code chronic medications into wide format
    med_chronic <- med_chronic[!duplicated(med_chronic[,c(1,2,4)]),]
    med_chronic <- med_chronic[,-c(2,3)]
    med_chronic$concept_code <- paste("old",med_chronic$concept_code,sep="_")
    med_chronic$present <- 1
    med_chronic <- med_chronic %>% tidyr::spread(concept_code,present)
    med_chronic[is.na(med_chronic)] <- 0
    
    # Create subtable for ACE-i/ARB pre-exposure
    acei_present = ("old_ACEI" %in% colnames(med_chronic))
    arb_present = ("old_ARB" %in% colnames(med_chronic))
    if(acei_present == TRUE && arb_present == TRUE) {
        med_acearb_chronic <- med_chronic %>% dplyr::select(patient_id,old_ACEI,old_ARB)
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = ifelse(old_ACEI + old_ARB > 0,1,0))
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::select(patient_id,acei_arb_preexposure)
    } else if (acei_present == TRUE) {
        med_acearb_chronic <- med_chronic %>% dplyr::select(patient_id,old_ACEI)
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = ifelse(old_ACEI > 0,1,0))
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::select(patient_id,acei_arb_preexposure)
    } else if (arb_present == TRUE) {
        med_acearb_chronic <- med_chronic %>% dplyr::select(patient_id,old_ARB)
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = ifelse(old_ARB > 0,1,0))
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::select(patient_id,acei_arb_preexposure)
    }
    
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
    aki_start_time <- labs_aki_summ[,c(1,5)]
    med_new_aki <- merge(med_new,aki_start_time,by="patient_id",all.x=TRUE)
    med_new_aki <- med_new_aki[!is.na(med_new_aki$day_min),]
    med_new_aki <- med_new_aki %>% dplyr::distinct()
    med_new_aki <- med_new_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(offset_aki = start_day - day_min) %>% dplyr::filter(offset_aki == min(offset_aki))
    # Re-code whether medication was given before AKI - 1 = yes, 0 = no
    med_new_aki <- med_new_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(med_before_aki = ifelse(offset_aki <=0,1,0))
    med_new_aki <- med_new_aki[,c(1,2,6)]
    med_new_aki <- med_new_aki %>% tidyr::spread(concept_code,med_before_aki)
    #med_new_aki[is.na(med_new_aki)] <- -999
    
    # Generate another table with the start date of the new medications
    med_new <- med_new %>% tidyr::spread(concept_code,start_day)
    #med_new[is.na(med_new)] <- -999
    
    # Generate simplified table for determining who were started on COAGB near admission
    med_coagb_new <- med_new %>% dplyr::select(patient_id,COAGB)
    med_coagb_new$COAGB[med_coagb_new$COAGB < -15] <- 0
    med_coagb_new$COAGB[med_coagb_new$COAGB >= -15] <- 1
    
    # Generate simplified table for determining who were started on novel antivirals
    covid19antiviral_present = ("COVIDVIRAL" %in% colnames(med_new))
    if(covid19antiviral_present == TRUE){
        med_covid19_new <- med_new %>% dplyr::select(patient_id,COVIDVIRAL)
        # Uncomment following line to select for patients who were started on novel antivirals less than 72h from admission
        #med_covid19_new <- med_covid19_new[med_covid19_new$COVIDVIRAL <= 3,]
        med_covid19_new <- med_covid19_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covid_rx = ifelse(COVIDVIRAL >= 0,1,0))
        med_covid19_new_date <- med_covid19_new
        colnames(med_covid19_new_date)[2] <- "covid_rx_start"
        med_covid19_new <- med_covid19_new %>% dplyr::select(patient_id,covid_rx)
    }
    
    # =====================
    # Demographics Table
    # =====================
    
    demog_summ <- demographics_filt %>% dplyr::select(patient_id,sex,age_group,race,severe,deceased,time_to_severe,time_to_death)
    demog_summ <- merge(demog_summ,comorbid,by="patient_id",all.x=TRUE)
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$aki <- 0
    demog_summ$aki[demog_summ$patient_id %in% labs_aki_summ$patient_id] <- 1
    demog_summ$severe <- factor(demog_summ$severe,levels=c(0,1),labels=c("Non-severe","Severe"))
    demog_summ$deceased <- factor(demog_summ$deceased,levels=c(0,1),labels=c("Alive","Deceased"))
    demog_summ$aki <- factor(demog_summ$aki,levels=c(0,1),labels=c("No AKI","AKI"))
    label(demog_summ$sex) <- "Sex"
    label(demog_summ$age_group) <- "Age Group"
    label(demog_summ$race) <- "Race"
    label(demog_summ$severe) <- "Severity"
    label(demog_summ$deceased) <- "Survival"
    label(demog_summ$time_to_severe) <- "Time to Severity Onset"
    label(demog_summ$time_to_death) <- "Time to Death"
    units(demog_summ$time_to_severe) <- "days"
    units(demog_summ$time_to_death) <- "days"
    
    demog_table <- table1::table1(~ sex + age_group + race + severe + deceased + time_to_severe + time_to_death | aki,data=demog_summ,overall="Total",render.continuous=FourCePhase2.1AKI:::my.render.cont,render.categorical=FourCePhase2.1AKI:::my.render.cat,export=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_Demographics_Summary.csv")))
    
    ## ==================================================================================
    ## PART 3: Serum Creatinine Trends - Plots against Time from Peak Serum Creatinine
    ## ==================================================================================

    # Severity indices
    # (1) non-severe, no AKI
    # (2) non-severe, AKI
    # (3) severe, no AKI
    # (4) severe, AKI before severity onset
    # (5) severe, AKI after severity onset
    
    # Novel antivirals indices
    # (1) non-severe, no novel anti-COVID-19 agents (i.e. severe == 1 or 2 && covidrx == 0)
    # (2) non-severe, with novel anti-COVID-19 agents (i.e. severe == 1 or 2 && covidrx == 1)
    # (3) severe, no novel anti-COVID-19 agents (i.e. severe == 3, 4 or 5 && covidrx == 0)
    # (4) severe, with novel anti-COVID-19 agents (i.e. severe == 3, 4 or 5 && covidrx == 1)
    
    # This analysis uses the index AKI episode
    # Feel free to modify the code to look at other types of AKI episodes, e.g. you can dplyr::filter for the most severe AKI episode
    # First, dplyr::filter the labs_aki_severe table to show only the index AKI episodes
    aki_only_index <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(days_since_admission == min(days_since_admission)) %>% dplyr::ungroup()
    # patient_id	site_id	days_since_admission	value	day_min	day_min_retro	min_cr_90d	min_cr_48h	min_cr_retro_7day	min_cr_48h_retro	min_cr_7d_final	cr_7d	cr_90d	delta_cr	aki_kdigo	aki_kdigo_retro	aki_kdigo_final	akd_7d	akd_90d	severe  time_to_severe	severe_to_aki	severe_before_aki
    
    # Generate the patient list including (1) severity indices from this dplyr::filtered table (2) day of peak Cr
    # severe - 2 = never severe, 4 = severe, AKI before severity onset, 5 = severe, AKI after severity onset
    aki_only_index <- aki_only_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(!is.na(time_to_severe),ifelse(severe_before_aki == 1,5,4),2)) %>% dplyr::ungroup()
    aki_only_index <- aki_only_index %>% dplyr::select(patient_id,days_since_admission,severe,day_min,severe_to_aki)
    colnames(aki_only_index)[2] <- "peak_cr_time"
    colnames(aki_only_index)[4] <- "aki_start"
    # Headers of aki_only_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki
    
    no_aki_list <- demographics_filt %>% dplyr::select(patient_id,severe)
    no_aki_list <- no_aki_list[!(no_aki_list$patient_id %in% aki_only_index$patient_id),]
    no_aki_list <- no_aki_list %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(severe == 1,3,1))
    
    labs_nonaki_summ <- labs_nonaki_summ[labs_nonaki_summ$patient_id %in% no_aki_list$patient_id,]
    labs_nonaki_severe <- labs_nonaki_severe[labs_nonaki_severe$patient_id %in% no_aki_list$patient_id,]
    labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$patient_id %in% no_aki_list$patient_id,]
    
    # Create a non-AKI equivalent for aki_only_index - except that this takes the largest delta_cr (and the earliest occurence of such a delta_cr)
    no_aki_index <- labs_nonaki_severe %>% dplyr::group_by(patient_id) %>% dplyr::filter(delta_cr == max(delta_cr)) %>% dplyr::filter(days_since_admission == min(days_since_admission)) %>% dplyr::ungroup()
    no_aki_index <- no_aki_index %>% dplyr::select(patient_id,days_since_admission,severe,day_min,severe_to_aki)
    no_aki_index <- no_aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(severe == 1,3,1))
    colnames(no_aki_index)[2] <- "peak_cr_time"
    colnames(no_aki_index)[4] <- "aki_start"
    #no_aki_index$severe_to_aki <- -999
    
    aki_index <- dplyr::bind_rows(aki_only_index,no_aki_index)
    if(covid19antiviral_present == TRUE) {
        aki_index <- merge(aki_index,med_covid19_new,by="patient_id",all.x=TRUE)
        aki_index$covid_rx[is.na(aki_index$covid_rx)] <- 0
        aki_index <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covidrx_grp = ifelse(severe <= 2, ifelse(covid_rx == 0,1,2),ifelse(covid_rx == 0,3,4))) %>% dplyr::ungroup()
        aki_index <- aki_index %>% dplyr::select(patient_id,peak_cr_time,severe,aki_start,severe_to_aki,covidrx_grp)
    } else {
        aki_index <- aki_index %>% dplyr::select(patient_id,peak_cr_time,severe,aki_start,severe_to_aki)
    }
    # Headers of aki_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki  covidrx_grp
    
    # Uncomment the following line to remove patients who were previously on RRT prior to admission
    # aki_index <- aki_index[!(aki_index$patient_id %in% patients_already_rrt),]
    
    # Create a common labs_cr_all table containing the serum Cr values, the severity groupings and anti-viral groupings
    labs_cr_aki_tmp <- labs_cr_aki %>% dplyr::select(patient_id,days_since_admission,value,min_cr_90d,min_cr_retro_7day)
    labs_cr_nonaki_tmp <- labs_cr_nonaki %>% dplyr::select(patient_id,days_since_admission,value,min_cr_90d,min_cr_retro_7day)
    labs_cr_all <- dplyr::bind_rows(labs_cr_aki_tmp,labs_cr_nonaki_tmp)
    labs_cr_all <- merge(labs_cr_all,aki_index,by="patient_id",all.x=TRUE)
    
    # Now, generate a table containing lab values with timepoints calculated from time of peak cr
    peak_trend <- labs_cr_all
    if(covid19antiviral_present == TRUE) {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,covidrx_grp,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_retro_7day)
    } else {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_retro_7day)
    }
    # patient_id  severe  covidrx_grp  days_since_admission  peak_cr_time  value min_cr_90d min_cr_retro_7day
    
    # Calculate the day from peak Cr
    peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_peak = days_since_admission - peak_cr_time) %>% dplyr::ungroup()
    ## dplyr::filter this table for Cr values that fall within the 7 day window
    # peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(between(time_from_peak,0,7)) %>% dplyr::ungroup()
    # Normalise to baseline values used for AKI calculation
    peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
    
    # In the event longitudinal data becomes very long, we will create a column where the very first baseline Cr for the index AKI is generated for each patient
    first_baseline <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(time_from_peak == 0) %>% dplyr::select(patient_id,baseline_cr) %>% dplyr::distinct()
    colnames(first_baseline)[2] <- "first_baseline_cr"
    peak_trend <- merge(peak_trend,first_baseline,by="patient_id",all.x=TRUE)
    if(covid19antiviral_present == TRUE) {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,covidrx_grp,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_retro_7day,time_from_peak,baseline_cr,first_baseline_cr)
    } else {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_retro_7day,time_from_peak,baseline_cr,first_baseline_cr)
    }
    peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/first_baseline_cr) %>% dplyr::ungroup()
    #peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
    
    # peak_trend will now be a common table to plot from the dplyr::selected AKI peak
    
    # =======================================================================================
    # Figure 1(a): Cr trends from start of AKI / after peak Cr for severe vs non-severe groups
    # =======================================================================================
    
    # First create a plot of the creatinine trends of AKI vs non-AKI patients from the day of the first AKI peak (or the highest Cr peak for non-AKI patients)
    peak_aki_vs_non_aki <- peak_trend %>% dplyr::select(patient_id,severe,time_from_peak,ratio)
    colnames(peak_aki_vs_non_aki) <- c("patient_id","aki","time_from_peak","ratio")
    peak_aki_vs_non_aki <- peak_aki_vs_non_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki = ifelse((aki == 2 | aki == 4 | aki == 5),1,0))
    peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki %>% dplyr::group_by(aki,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    aki_label <- data.table::data.table(c(0,1),c("Non-AKI","AKI"))
    colnames(aki_label) <- c("aki","aki_label")
    peak_aki_vs_non_aki_summ <- merge(peak_aki_vs_non_aki_summ,aki_label,by="aki",all.x=TRUE)
    if(is_obfuscated==TRUE) {
        peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki_summ %>% dplyr::group_by(aki,time_from_peak) %>% dplyr::filter(n >= obfuscation_value)
    }
    write.csv(peak_aki_vs_non_aki_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_CrFromPeak_AKI_vs_NonAKI.csv")),row.names=FALSE)
    peak_aki_vs_non_aki_timeplot <- ggplot2::ggplot(peak_aki_vs_non_aki_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=aki))+ggplot2::geom_line(ggplot2::aes(color = factor(aki))) + ggplot2::geom_point(ggplot2::aes(color = factor(aki))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "AKI Group") + ggplot2::xlim(-30,30)
    print(peak_aki_vs_non_aki_timeplot)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_AKI_vs_NonAKI.png")),plot=peak_aki_vs_non_aki_timeplot,width=12,height=9,units="cm")
    

    # Now, derive our first table peak_trend_severe to compare across the different severity groups
    peak_trend_severe <- peak_trend %>% dplyr::select(patient_id,severe,time_from_peak,ratio)
    # Headers: patient_id  severe  time_from_peak  ratio
    peak_trend_severe <- peak_trend_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
    
    # Calculate mean and SD each for severe and non-severe groups
    peak_cr_summ <- peak_trend_severe %>% dplyr::group_by(severe,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    severe_label <- data.table::data.table(c(1,2,3,4),c("Non-severe, no AKI","Non-severe, AKI","Severe, no AKI","Severe, AKI"))
    colnames(severe_label) <- c("severe","severe_label")
    peak_cr_summ <- merge(peak_cr_summ,severe_label,by="severe",all.x=TRUE)
    if(is_obfuscated==TRUE) {
        peak_cr_summ <- peak_cr_summ %>% dplyr::group_by(severe) %>% dplyr::filter(n >= obfuscation_value)
    }
    write.csv(peak_cr_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_CrfromPeak_Severe_AKI.csv")),row.names=FALSE)
    # Plot the graphs
    peak_cr_timeplot <- ggplot2::ggplot(peak_cr_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30)
    print(peak_cr_timeplot)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromPeak_Severe_AKI.png")),plot=peak_cr_timeplot,width=12,height=9,units="cm")
    
    # Plot from start of admission to 30 days post-peak AKI (if no AKI, then from peak Cr)
    adm_to_aki_cr <- labs_cr_all
    #adm_to_aki_cr$peak_cr_time[is.na(adm_to_aki_cr$peak_cr_time)] <- 0
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::filter(between(days_since_admission,0,peak_cr_time+30)) %>% dplyr::ungroup()
    adm_to_aki_cr <- adm_to_aki_cr[order(adm_to_aki_cr$patient_id,adm_to_aki_cr$days_since_admission),]
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
    adm_to_aki_summ <- adm_to_aki_cr %>% dplyr::group_by(severe,days_since_admission) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    adm_to_aki_summ <- merge(adm_to_aki_summ,severe_label,by="severe",all.x=TRUE)
    if(is_obfuscated==TRUE) {
        adm_to_aki_summ <- adm_to_aki_summ %>% dplyr::group_by(severe) %>% dplyr::filter(n >= obfuscation_value)
    }
    write.csv(peak_cr_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromAdmToPeak+10D_Severe_AKI.csv")),row.names=FALSE)
    adm_to_aki_timeplot <- ggplot2::ggplot(adm_to_aki_summ[which(adm_to_aki_summ$days_since_admission <= 30),],ggplot2::aes(x=days_since_admission,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity")
    print(adm_to_aki_timeplot)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromAdmToPeak+10D_Severe_AKI.png")),plot=adm_to_aki_timeplot,width=12,height=9,units="cm")
    
    # Plot from start of AKI to 30 days later 
    
    aki_30d_cr <- labs_cr_all
    # Uncomment the following line to restrict analysis to AKI patients only
    #aki_30d_cr <- aki_30d_cr[aki_30d_cr$severe %in% c(2,4,5),]
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_start = days_since_admission - aki_start) %>% dplyr::ungroup()
    aki_30d_cr <- aki_30d_cr[order(aki_30d_cr$patient_id,aki_30d_cr$days_since_admission),]
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
    aki_30d_cr_summ <- aki_30d_cr %>% dplyr::group_by(severe,time_from_start) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    aki_30d_cr_summ <- merge(aki_30d_cr_summ,severe_label,by="severe",all.x=TRUE)
    if(is_obfuscated==TRUE) {
        aki_30d_cr_summ <- aki_30d_cr_summ %>% dplyr::group_by(severe) %>% dplyr::filter(n >= obfuscation_value)
    }
    write.csv(aki_30d_cr_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_Severe_AKI.csv")),row.names=FALSE)
    aki_30d_cr_timeplot <- ggplot2::ggplot(aki_30d_cr_summ,ggplot2::aes(x=time_from_start,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30)
    print(aki_30d_cr_timeplot)
    ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_Severe_AKI.png")),plot=aki_30d_cr_timeplot,width=12,height=9,units="cm")
    
    # ================================================================================================================================
    # Figure 1(b) Comparing serum creatinine trends of severe and non-severe patients, with or without remdesivir/lopinavir+ritonavir
    # ================================================================================================================================
    if(covid19antiviral_present == TRUE) {
        # Plotting from peak
        peak_trend_covidviral <- peak_trend %>% dplyr::select(patient_id,covidrx_grp,time_from_peak,ratio)
        peak_cr_covidviral_summ <- peak_trend_covidviral %>% dplyr::group_by(covidrx_grp,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n())) %>% dplyr::ungroup()
        if(is_obfuscated==TRUE) {
            peak_cr_covidviral_summ  <- peak_cr_covidviral_summ  %>% dplyr::group_by(covidrx_grp) %>% dplyr::filter(n >= obfuscation_value)
        }
        write.csv(peak_cr_covidviral_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromPeak_CovidViral.csv")),row.names=FALSE)
        peak_cr_covidviral_timeplot <- ggplot2::ggplot(peak_cr_covidviral_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=covidrx_grp))+ggplot2::geom_line(ggplot2::aes(color = factor(covidrx_grp))) + ggplot2::geom_point(ggplot2::aes(color = factor(covidrx_grp))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity + COVID-19 Treatment") + ggplot2::xlim(-30,30)
        print(peak_cr_covidviral_timeplot)
        ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromPeak_CovidViral.png")),plot=peak_cr_covidviral_timeplot,width=12,height=9,units="cm")
        
        # Plotting from initiation of novel anti-virals
        cr_from_covidrx_trend <- merge(peak_trend,med_covid19_new_date,by="patient_id",all.x=TRUE)
        # dplyr::filter out patients who have never received any of the novel antivirals
        #cr_from_covidrx_trend$covid_rx_start[is.na(cr_from_covidrx_trend$covid_rx_start)] <- -999
        cr_from_covidrx_trend$covid_rx[is.na(cr_from_covidrx_trend$covid_rx)] <- 0
        cr_from_covidrx_trend <- cr_from_covidrx_trend[cr_from_covidrx_trend$covid_rx == 1,]
        
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::select(patient_id,severe,covidrx_grp,days_since_admission,covid_rx_start,peak_cr_time,value,min_cr_90d,min_cr_retro_7day)
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(severe <= 2,0,1)) %>% dplyr::mutate(covidrx_grp = ifelse((covidrx_grp == 2 || covidrx_grp == 4),1,0)) %>% dplyr::ungroup()
        # Calculate the day from initiation of novel anti-virals
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_covidrx = days_since_admission - ifelse(covidrx_grp == 1,covid_rx_start,0)) %>% dplyr::ungroup()
        # Filter this table for Cr values that fall within the 7 day window
        # cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(between(time_from_covidrx,0,30)) %>% dplyr::ungroup()
        # Normalise to baseline values used for AKI calculation
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
        cr_from_covidrx_trend_severe <- cr_from_covidrx_trend %>% dplyr::select(patient_id,severe,time_from_covidrx,ratio)
        # Headers: patient_id  severe (coded as 0/1)  time_from_covidrx  ratio
        cr_from_covidrx_summ <- cr_from_covidrx_trend_severe %>% dplyr::group_by(severe,time_from_covidrx) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n())) %>% dplyr::ungroup()
        if(is_obfuscated==TRUE) {
            cr_from_covidrx_summ  <- cr_from_covidrx_summ %>% dplyr::group_by(severe) %>% dplyr::filter(n >= obfuscation_value)
        }
        
        write.csv(peak_cr_covidviral_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromCovidRx_Severe.csv")),row.names=FALSE)
        cr_from_covidrx_timeplot <- ggplot2::ggplot(cr_from_covidrx_summ,ggplot2::aes(x=time_from_covidrx,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from COVIDVIRAL Start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30)
        print(cr_from_covidrx_timeplot)
        ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromCovidRx_Severe.png")),plot=cr_from_covidrx_timeplot,width=12,height=9,units="cm")
    }
    
    ## ====================================
    ## PART 4: Time To Event Analysis
    ## ====================================
    
    labs_cr_recovery <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(time_from_peak >= 0)
    labs_cr_recovery_tmp <- labs_cr_recovery %>% dplyr::group_by(patient_id) %>% tidyr::complete(time_from_peak = tidyr::full_seq(time_from_peak,1)) %>% dplyr::mutate(ratio = zoo::na.fill(ratio,Inf))
    time_to_ratio1.25 <- labs_cr_recovery_tmp %>% split(.$patient_id) %>% purrr::map(~get_day(.$ratio,.$time_from_peak,target=1.25)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_id')
    colnames(time_to_ratio1.25)[2] <- "time_to_ratio1.25"
    
    labs_aki_summ_index <- labs_aki_summ %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(days_since_admission == min(days_since_admission))
    index_aki_grade <- labs_aki_summ_index %>% dplyr::select(patient_id,aki_kdigo_final)
    aki_index_recovery <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::filter(severe %in% c(2,4,5)) %>% dplyr::mutate(severe=ifelse(severe==2,0,1))
    aki_index_recovery <- merge(aki_index_recovery,time_to_ratio1.25,by="patient_id",all.x=TRUE)
    aki_index_recovery <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(recover_1.25x = ifelse(is.na(time_to_ratio1.25),0,1))
    
    discharge_day <- demographics %>% dplyr::select(patient_id,admission_date,last_discharge_date,deceased) %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_death_km = as.numeric(as.Date(last_discharge_date)-as.Date(admission_date))) %>% dplyr::select(patient_id,deceased,time_to_death_km)
    aki_index_recovery <- merge(aki_index_recovery,discharge_day,by="patient_id",all.x=TRUE)
    aki_index_recovery <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_ratio1.25 = dplyr::if_else(recover_1.25x == 0,time_to_death_km,time_to_ratio1.25))
    aki_index_recovery <- merge(aki_index_recovery,labs_aki_summ_index[,c(1,17)],by="patient_id",all.x=TRUE)
    
    comorbid_list <- colnames(comorbid)[-1]
    aki_index_recovery <- merge(aki_index_recovery,comorbid,by="patient_id",all.x=TRUE)
    
    surv_recover <- survival::Surv(time=aki_index_recovery$time_to_ratio1.25,event=aki_index_recovery$recover_1.25x)
    fit_km_recover <- survival::survfit(surv_recover ~ severe, data=aki_index_recovery)
    plot_recover <- survminer::ggsurvplot(fit_km_recover,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event")
    plot_recover_summ <- survminer::surv_summary(fit_km_recover,data=aki_index_recovery)
    plot_recover_summ_table <- plot_recover$data.survtable
    write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_Plot.csv")),row.names=FALSE)
    write.csv(plot_recover_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_Table.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe.png")),plot=print(plot_recover),width=9,height=12,units="cm")
    coxph_recover <- survival::coxph(as.formula(paste("surv_recover ~ ",paste(c("severe","aki_kdigo_final",comorbid_list),collapse="+"))), data=aki_index_recovery)
    coxph_recover_plot <- survminer::ggforest(coxph_recover,data=aki_index_recovery)
    coxph_recover_summ <- summary(coxph_recover) 
    write.csv(coxph_recover_summ$coefficients,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH.png")),plot=print(coxph_recover_plot),width=10,height=5,units="cm")
    
    surv_death_aki_only<- survival::Surv(time=aki_index_recovery$time_to_death_km,event=aki_index_recovery$deceased)
    fit_death_aki_only <- survival::survfit(surv_death_aki_only ~ severe, data=aki_index_recovery)
    plot_death_aki_only <- survminer::ggsurvplot(fit_death_aki_only,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_death_aki_only_summ <- survminer::surv_summary(fit_death_aki_only,data=aki_index_recovery)
    plot_death_aki_only_summ_table <- plot_death_aki_only$data.survtable
    write.csv(plot_death_aki_only_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe_Plot.csv")),row.names=FALSE)
    write.csv(plot_death_aki_only_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe_Table.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe.png")),plot=print(plot_death_aki_only),width=9,height=12,units="cm")
    coxph_death_aki_only <- survival::coxph(as.formula(paste("surv_death_aki_only ~ ",paste(c("severe","aki_kdigo_final",comorbid_list),collapse="+"))), data=aki_index_recovery)
    coxph_death_aki_only_plot <- survminer::ggforest(coxph_death_aki_only,data=aki_index_recovery)
    coxph_death_aki_only_summ <- summary(coxph_death_aki_only) 
    write.csv(coxph_death_aki_only_summ$coefficients,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH.png")),plot=print(coxph_death_aki_only_plot),width=10,height=5,units="cm")
    
    aki_index_death <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(is_aki=ifelse(severe %in% c(3,4,5),1,0)) %>% dplyr::mutate(severe=ifelse(severe %in% c(2,4,5),1,0))
    aki_index_death <- merge(aki_index_death,discharge_day,by="patient_id",all.x=TRUE)
    aki_index_death <- merge(aki_index_death,labs_aki_summ_index[,c(1,17)],by="patient_id",all.x=TRUE)
    aki_index_death <- merge(aki_index_death,comorbid,by="patient_id",all.x=TRUE)
    
    surv_death <- survival::Surv(time=aki_index_death$time_to_death_km,event=aki_index_death$deceased)
    fit_death <- survival::survfit(surv_death ~ severe, data=aki_index_death)
    plot_death <- survminer::ggsurvplot(fit_death,data=aki_index_death,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_death_summ <- survminer::surv_summary(fit_death,data=aki_index_death)
    plot_death_summ_table <- plot_death$data.survtable
    write.csv(plot_death_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_Severe_Plot.csv")),row.names=FALSE)
    write.csv(plot_death_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_Severe_Table.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_Severe.png")),plot=print(plot_death),width=9,height=12,units="cm")
    coxph_death <- survival::coxph(as.formula(paste("surv_death ~ ",paste(c("severe","is_aki",comorbid_list),collapse="+"))), data=aki_index_death)
    coxph_death_plot <- survminer::ggforest(coxph_death,data=aki_index_death)
    coxph_death_summ <- summary(coxph_death) 
    write.csv(coxph_death_summ$coefficients,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_CoxPH.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_CoxPH.png")),plot=print(coxph_death_plot),width=10,height=5,units="cm")
    
}

