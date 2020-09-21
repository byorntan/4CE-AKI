library(dplyr)

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
# Feel free to modify the code to look at other types of AKI episodes, e.g. you can filter for the most severe AKI episode
# First, filter the labs_aki_severe table to show only the index AKI episodes
aki_only_index <- labs_aki_severe %>% group_by(patient_id) %>% filter(days_since_admission == min(days_since_admission)) %>% ungroup()
# patient_id	site_id	days_since_admission	value	day_min	day_min_retro	min_cr_90d	min_cr_48h	min_cr_7day_retro	min_cr_48h_retro	min_cr_7d_final	cr_7d	cr_90d	delta_cr	aki_kdigo	aki_kdigo_retro	aki_kdigo_final	akd_7d	akd_90d	severe  time_to_severe	severe_to_aki	severe_before_aki

# Generate the patient list including (1) severity indices from this filtered table (2) day of peak Cr
# severe - 2 = never severe, 4 = severe, AKI before severity onset, 5 = severe, AKI after severity onset
aki_only_index <- aki_only_index %>% group_by(patient_id) %>% mutate(severe = ifelse(time_to_severe != -999,ifelse(severe_before_aki == 1,5,4),2)) %>% ungroup()
aki_only_index <- aki_only_index %>% select(patient_id,days_since_admission,severe,day_min,severe_to_aki)
colnames(aki_only_index)[2] <- "peak_cr_time"
colnames(aki_only_index)[4] <- "aki_start"
# Headers of aki_only_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki

no_aki_list <- demographics_filt %>% select(patient_id,severe)
no_aki_list <- no_aki_list[!(no_aki_list$patient_id %in% aki_only_index$patient_id),]
no_aki_list <- no_aki_list %>% group_by(patient_id) %>% mutate(severe = ifelse(severe == 1,3,1))

labs_nonaki_summ <- labs_nonaki_summ[labs_nonaki_summ$patient_id %in% no_aki_list$patient_id,]
labs_nonaki_severe <- labs_nonaki_severe[labs_nonaki_severe$patient_id %in% no_aki_list$patient_id,]
labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$patient_id %in% no_aki_list$patient_id,]

# Create a non-AKI equivalent for aki_only_index - except that this takes the largest delta_cr (and the earliest occurence of such a delta_cr)
no_aki_index <- labs_nonaki_severe %>% group_by(patient_id) %>% filter(delta_cr == max(delta_cr)) %>% filter(days_since_admission == min(days_since_admission)) %>% ungroup()
no_aki_index <- no_aki_index %>% select(patient_id,days_since_admission,severe,day_min,severe_to_aki)
no_aki_index <- no_aki_index %>% group_by(patient_id) %>% mutate(severe = ifelse(severe == 1,3,1))
colnames(no_aki_index)[2] <- "peak_cr_time"
colnames(no_aki_index)[4] <- "aki_start"
no_aki_index$severe_to_aki <- -999

aki_index <- bind_rows(aki_only_index,no_aki_index)

aki_index <- merge(aki_index,med_covid19_new,by="patient_id",all.x=TRUE)
aki_index$covid_rx[is.na(aki_index$covid_rx)] <- 0
aki_index <- aki_index %>% group_by(patient_id) %>% mutate(covidrx_grp = ifelse(severe <= 2, ifelse(covid_rx == 0,1,2),ifelse(covid_rx == 0,3,4))) %>% ungroup()
aki_index <- aki_index %>% select(patient_id,peak_cr_time,severe,aki_start,severe_to_aki,covidrx_grp)
# Headers of aki_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki  covidrx_grp

# Uncomment the following line to remove patients who were previously on RRT prior to admission
# aki_index <- aki_index[!(aki_index$patient_id %in% patients_already_rrt),]

# Create a common labs_cr_all table containing the serum Cr values, the severity groupings and anti-viral groupings
labs_cr_aki_tmp <- labs_cr_aki %>% select(patient_id,days_since_admission,value,min_cr_90d,min_cr_7day_retro)
labs_cr_nonaki_tmp <- labs_cr_nonaki %>% select(patient_id,days_since_admission,value,min_cr_90d,min_cr_7day_retro)
labs_cr_all <- bind_rows(labs_cr_aki_tmp,labs_cr_nonaki_tmp)
labs_cr_all <- merge(labs_cr_all,aki_index,by="patient_id",all.x=TRUE)

# Now, generate a table containing lab values with timepoints calculated from time of peak cr
peak_trend <- labs_cr_all
peak_trend <- peak_trend %>% select(patient_id,severe,covidrx_grp,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_7day_retro)
# patient_id  severe  covidrx_grp  days_since_admission  peak_cr_time  value min_cr_90d min_cr_7day_retro

# Calculate the day from peak Cr
peak_trend <- peak_trend %>% group_by(patient_id) %>% mutate(time_from_peak = days_since_admission - peak_cr_time) %>% ungroup()
## Filter this table for Cr values that fall within the 7 day window
# peak_trend <- peak_trend %>% group_by(patient_id) %>% filter(between(time_from_peak,0,7)) %>% ungroup()
# Normalise to baseline values used for AKI calculation
peak_trend <- peak_trend %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_90d,min_cr_7day_retro)) %>% ungroup()

# In the event longitudinal data becomes very long, we will create a column where the very first baseline Cr for the index AKI is generated for each patient
first_baseline <- peak_trend %>% group_by(patient_id) %>% filter(time_from_peak == 0) %>% select(patient_id,baseline_cr) %>% distinct()
colnames(first_baseline)[2] <- "first_baseline_cr"
peak_trend <- merge(peak_trend,first_baseline,by="patient_id",all.x=TRUE)
peak_trend <- peak_trend %>% select(patient_id,severe,covidrx_grp,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_7day_retro,time_from_peak,baseline_cr,first_baseline_cr)

peak_trend <- peak_trend %>% group_by(patient_id) %>% mutate(ratio = value/first_baseline_cr) %>% ungroup()
#peak_trend <- peak_trend %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()

# peak_trend will now be a common table to plot from the selected AKI peak

# =======================================================================================
# Figure 1(a): Cr trends from start of AKI / after peak Cr for severe vs non-severe groups
# =======================================================================================

# First create a plot of the creatinine trends of AKI vs non-AKI patients from the day of the first AKI peak (or the highest Cr peak for non-AKI patients)
peak_aki_vs_non_aki <- peak_trend %>% select(patient_id,severe,time_from_peak,ratio)
colnames(peak_aki_vs_non_aki) <- c("patient_id","aki","time_from_peak","ratio")
peak_aki_vs_non_aki <- peak_aki_vs_non_aki %>% group_by(patient_id) %>% mutate(aki = ifelse((aki == 2 | aki == 4 | aki == 5),1,0))
peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki %>% group_by(aki,days_since_admission) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n()),n=n()) %>% ungroup()
aki_label <- data.table(c(0,1),c("Non-AKI","AKI"))
colnames(aki_label) <- c("aki","aki_label")
peak_aki_vs_non_aki_summ <- merge(peak_aki_vs_non_aki_summ,aki_label,by="aki",all.x=TRUE)
write.csv(peak_aki_vs_non_aki_summ,"labs_cr_fromPeak_aki_vs_non_aki.csv",row.names=FALSE)
peak_aki_vs_non_aki_timeplot <- ggplot(peak_aki_vs_non_aki_summ,aes(x=time_from_peak,y=mean_ratio,group=aki))+geom_line(aes(color = factor(aki))) + geom_point(aes(color = factor(aki))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "AKI Group")
print(peak_aki_vs_non_aki_timeplot)
ggsave("peak_Cr_AKI_vs_NonAKI.png",plot=peak_aki_vs_non_aki_timeplot)

# Now, derive our first table peak_trend_severe to compare across the different severity groups
peak_trend_severe <- peak_trend %>% select(patient_id,severe,time_from_peak,ratio)
# Headers: patient_id  severe  time_from_peak  ratio
peak_trend_severe <- peak_trend_severe %>% group_by(patient_id) %>% mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))

# Calculate mean and SD each for severe and non-severe groups
peak_cr_summ <- peak_trend_severe %>% group_by(severe,time_from_peak) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n()),n=n()) %>% ungroup()
severe_label <- data.table(c(1,2,3,4),c("Non-severe, no AKI","Non-severe, AKI","Severe, no AKI","Severe, AKI"))
colnames(severe_label) <- c("severe","severe_label")
peak_cr_summ <- merge(peak_cr_summ,severe_label,by="severe",all.x=TRUE)
write.csv(peak_cr_summ,"labs_cr_fromPeak_Severe+AKI.csv",row.names=FALSE)
# Plot the graphs
peak_cr_timeplot <- ggplot(peak_cr_summ,aes(x=time_from_peak,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity")
print(peak_cr_timeplot)

# Plot from start of admission to 10 days post-peak AKI (if no AKI, then from peak Cr)
adm_to_aki_cr <- labs_cr_all
#adm_to_aki_cr$peak_cr_time[is.na(adm_to_aki_cr$peak_cr_time)] <- 0
adm_to_aki_cr <- adm_to_aki_cr %>% group_by(patient_id) %>% filter(between(days_since_admission,0,peak_cr_time+10)) %>% ungroup()
adm_to_aki_cr <- adm_to_aki_cr[order(adm_to_aki_cr$patient_id,adm_to_aki_cr$days_since_admission),]
adm_to_aki_cr <- adm_to_aki_cr %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_90d,min_cr_7day_retro)) %>% ungroup()
adm_to_aki_cr <- adm_to_aki_cr %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
adm_to_aki_summ <- adm_to_aki_cr %>% group_by(severe,days_since_admission) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n()),n=n()) %>% ungroup()
adm_to_aki_summ <- adm_to_aki_summ %>% group_by(patient_id) %>% mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
adm_to_aki_summ <- merge(adm_to_aki_summ,severe_label,by="severe",all.x=TRUE)
write.csv(peak_cr_summ,"labs_cr_fromAdmToPeak+10D_Severe+AKI.csv",row.names=FALSE)
adm_to_aki_timeplot <- ggplot(adm_to_aki_summ[which(adm_to_aki_summ$days_since_admission <= 10),],aes(x=days_since_admission,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity")
print(adm_to_aki_timeplot)

# Plot from start of AKI to 10 days later 

aki_10d_cr <- labs_cr_all
# Uncomment the following line to restrict analysis to AKI patients only
#aki_10d_cr <- aki_10d_cr[aki_10d_cr$severe %in% c(2,4,5),]
aki_10d_cr <- aki_10d_cr %>% group_by(patient_id) %>% mutate(time_from_start = days_since_admission - aki_start) %>% ungroup()
aki_10d_cr <- aki_10d_cr[order(aki_10d_cr$patient_id,aki_10d_cr$days_since_admission),]
aki_10d_cr <- aki_10d_cr %>% group_by(patient_id) %>% filter(between(time_from_start,-10,10)) %>% ungroup()
aki_10d_cr <- aki_10d_cr %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_90d,min_cr_7day_retro)) %>% ungroup()
aki_10d_cr <- aki_10d_cr %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
aki_10d_cr_summ <- aki_10d_cr %>% group_by(severe,time_from_start) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n()),n=n()) %>% ungroup()
aki_10d_cr_summ <- aki_10d_cr_summ %>% group_by(patient_id) %>% mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
aki_10d_cr_summ <- merge(aki_10d_cr_summ,severe_label,by="severe",all.x=TRUE)
write.csv(aki_10d_cr_summ,"labs_cr_fromStart_Severe+AKI.csv",row.names=FALSE)
aki_10d_cr_timeplot <- ggplot(aki_10d_cr_summ,aes(x=time_from_start,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity")
print(aki_10d_cr_timeplot)

# ================================================================================================================================
# Figure 1(b) Comparing serum creatinine trends of severe and non-severe patients, with or without remdesivir/lopinavir+ritonavir
# ================================================================================================================================
# Plotting from peak
peak_trend_covidviral <- peak_trend %>% select(patient_id,covidrx_grp,time_from_peak,ratio)
peak_cr_covidviral_summ <- peak_trend_covidviral %>% group_by(covidrx_grp,time_from_peak) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()
peak_cr_covidviral_timeplot <- ggplot(peak_cr_covidviral_summ,aes(x=time_from_peak,y=mean_ratio,group=covidrx_grp))+geom_line(aes(color = factor(covidrx_grp))) + geom_point(aes(color = factor(covidrx_grp))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity + COVID-19 Treatment")
print(peak_cr_covidviral_timeplot)                                             

# Plotting from initiation of novel anti-virals
cr_from_covidrx_trend <- merge(peak_trend,med_covid19_new_date,by="patient_id",all.x=TRUE)
# Filter out patients who have never received any of the novel antivirals
cr_from_covidrx_trend$covid_rx_start[is.na(cr_from_covidrx_trend$covid_rx_start)] <- -999
cr_from_covidrx_trend$covid_rx[is.na(cr_from_covidrx_trend$covid_rx)] <- 0
cr_from_covidrx_trend <- cr_from_covidrx_trend[cr_from_covidrx_trend$covid_rx == 1,]

cr_from_covidrx_trend <- cr_from_covidrx_trend %>% select(patient_id,severe,covidrx_grp,days_since_admission,covid_rx_start,peak_cr_time,value,min_cr_90d,min_cr_7day_retro)
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% group_by(patient_id) %>% mutate(severe = ifelse(severe <= 2,0,1)) %>% mutate(covidrx_grp = ifelse((covidrx_grp == 2 || covidrx_grp == 4),1,0)) %>% ungroup()
# Calculate the day from initiation of novel anti-virals
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% group_by(patient_id) %>% mutate(time_from_covidrx = days_since_admission - ifelse(covidrx_grp == 1,covid_rx_start,0)) %>% ungroup()
# Filter this table for Cr values that fall within the 7 day window
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% group_by(patient_id) %>% filter(between(time_from_covidrx,0,7)) %>% ungroup()
# Normalise to baseline values used for AKI calculation
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_90d,min_cr_7day_retro)) %>% ungroup()
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
cr_from_covidrx_trend_severe <- cr_from_covidrx_trend %>% select(patient_id,severe,time_from_covidrx,ratio)
# Headers: patient_id  severe (coded as 0/1)  time_from_covidrx  ratio
cr_from_covidrx_summ <- cr_from_covidrx_trend_severe %>% group_by(severe,time_from_covidrx) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()
cr_from_covidrx_timeplot <- ggplot(cr_from_covidrx_summ,aes(x=time_from_covidrx,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from COVIDVIRAL Start",y = "Serum Cr/Baseline Cr", color = "Severity")
print(cr_from_covidrx_timeplot)

# =================================================================================================
# Figure 1(c) Kaplan-Meier plot for time to renal recovery (defined as time to reach ratio < 1.25)
# =================================================================================================

# ============================================================================================
# Figure 1(d) Comparing AKI plots for patients with CKD vs non-CKD (eGFR >= 60 vs eGFR < 60)
# ============================================================================================
# Works in progress as exact age is not available in data extracted for Phase 2.0 (only broad age groups)
egdr_ckdepi <- function(sex,age,race=0,creatinine) {
  creatinine <- creatinine * 0.0113 #convert to mg/dL for formula
  race_factor = 1
  sex_factor = 1
  k_factor = 0.9
  a_factor = -0.411
  if(race == 1) {
    race_factor <- 1.159
  }
  if(sex == 1) {
    sex_factor <- 1.018
    k_factor <- 0.7
    a_factor <- -0.329
  }
  egfr <- 141 * (min(creatinine/k_factor,1))^a_factor * (max(creatinine/k_factor,1))^(-1.209) * 0.993^age * sex_factor * race_factor
  return(egfr)
}

