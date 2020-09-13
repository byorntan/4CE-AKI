library(dplyr)
# =======================================================================================
# Figure 1: Cr trends from start of AKI / after peak Cr for severe vs non-severe groups
# =======================================================================================
# Severity indices to use here - five patient groups
# (1) non-severe, no AKI
# (2) non-severe, AKI
# (3) severe, no AKI
# (4) severe, AKI before severity onset
# (5) severe, AKI after severity onset

# This analysis uses the index AKI episode
# Feel free to modify the code to look at other types of AKI episodes, e.g. you can filter for the most severe AKI episode
# First, filter the labs_aki_severe table to show only the index AKI episodes
aki_only_index <- labs_aki_severe %>% group_by(patient_id) %>% filter(days_since_admission == min(days_since_admission)) %>% ungroup()
# patient_id	site_id	days_since_admission	value	day_min	day_min_retro	min_cr_7day	min_cr_48h	min_cr_7day_retro	min_cr_48h_retro	min_cr_7d_final	cr_7d	cr_90d	delta_cr	aki_kdigo	aki_kdigo_retro	aki_kdigo_final	akd_7d	akd_90d	time_to_severe	severe_to_aki	severe_before_aki

# Generate the patient list including (1) severity indices from this filtered table (2) day of peak Cr
# severe - 2 = never severe, 4 = severe, AKI before severity onset, 5 = severe, AKI after severity onset
aki_only_index <- aki_only_index %>% group_by(patient_id) %>% mutate(severe = ifelse(time_to_severe != -999,ifelse(severe_before_aki == 1,5,4),2)) %>% ungroup()
aki_only_index <- aki_only_index[,c(1,3,23,5,21)]
colnames(aki_only_index)[2] <- "peak_cr_time"
colnames(aki_only_index)[4] <- "aki_start"
# Headers of aki_only_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki

no_aki_list <- demographics_filt[,c(1,7)]
no_aki_list <- no_aki_list[!(no_aki_list$patient_id %in% aki_only_index$patient_id),]
no_aki_list <- no_aki_list %>% group_by(patient_id) %>% mutate(severe = ifelse(severe == 1,3,1))
no_aki_list$peak_cr_time <- 0
no_aki_list$aki_start <- 0
no_aki_list$severe_to_aki <- -999
no_aki_list <- no_aki_list[,c(1,3,2,4,5)]

aki_index <- bind_rows(aki_only_index,no_aki_list)

# Uncomment the following line to remove patients who were previously on RRT prior to admission
# aki_index <- aki_index[!(aki_index$patient_id %in% patients_already_rrt),]

# Now bind this list of days of peak Cr to a new table containing the Cr trends
peak_aki_trend <- labs_cr_aki[,c(1,3,4,5,7)]
peak_aki_trend <- merge(peak_aki_trend,aki_index,by="patient_id",all.x=TRUE)
peak_aki_trend <- peak_aki_trend %>% select(patient_id,severe,days_since_admission,peak_cr_time,value,min_cr_7day,min_cr_7day_retro)
# patient_id  severe  days_since_admission  peak_cr_time  value min_cr_7day min_cr_7day_retro severe_to_aki

# Calculate the day from peak Cr
peak_aki_trend <- peak_aki_trend %>% group_by(patient_id) %>% mutate(time_from_peak = days_since_admission - peak_cr_time) %>% ungroup()
# Filter this table for Cr values that fall within the 7 day window
peak_aki_trend <- peak_aki_trend %>% group_by(patient_id) %>% filter(between(time_from_peak,0,7)) %>% ungroup()
# Normalise to baseline values used for AKI calculation
peak_aki_trend <- peak_aki_trend %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_7day,min_cr_7day_retro)) %>% ungroup()
peak_aki_trend <- peak_aki_trend %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
peak_aki_trend <- peak_aki_trend %>% select(patient_id,severe,time_from_peak,ratio)
# Headers: patient_id  severe  time_from_peak  ratio

# Calculate mean and SD each for severe and non-severe groups
peak_cr_summ <- peak_aki_trend %>% group_by(severe,time_from_peak) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()

# Plot the graphs
peak_cr_timeplot <- ggplot(peak_cr_summ,aes(x=time_from_peak,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity")
print(peak_cr_timeplot)

# Plot from start of admission to 10 days post-peak AKI (if no peak, then up to 10 days)
adm_to_aki_cr <- labs_cr_aki[,c(1,3,4,5,7)]
adm_to_aki_cr <- merge(adm_to_aki_cr,aki_index,by="patient_id",all.x=TRUE)
#adm_to_aki_cr$peak_cr_time[is.na(adm_to_aki_cr$peak_cr_time)] <- 0
adm_to_aki_cr <- adm_to_aki_cr %>% group_by(patient_id) %>% filter(between(days_since_admission,0,peak_cr_time+10)) %>% ungroup()
adm_to_aki_cr <- adm_to_aki_cr[order(adm_to_aki_cr$patient_id,adm_to_aki_cr$days_since_admission),]
adm_to_aki_cr <- adm_to_aki_cr %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_7day,min_cr_7day_retro)) %>% ungroup()
adm_to_aki_cr <- adm_to_aki_cr %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
adm_to_aki_summ <- adm_to_aki_cr %>% group_by(severe,days_since_admission) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()
adm_to_aki_timeplot <- ggplot(adm_to_aki_summ[which(adm_to_aki_summ$days_since_admission <= 10),],aes(x=days_since_admission,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity")
print(adm_to_aki_timeplot)

# Plot from start of AKI to 10 days later 
aki_10d_cr <- labs_cr_aki[,c(1,3,4,5,7)]
aki_10d_cr <- merge(aki_10d_cr,aki_index,by="patient_id",all.x=TRUE)
# Uncomment the following line to restrict analysis to AKI patients only
#aki_10d_cr <- aki_10d_cr[aki_10d_cr$severe %in% c(2,4,5),]
aki_10d_cr <- aki_10d_cr %>% group_by(patient_id) %>% mutate(time_from_start = days_since_admission - aki_start) %>% ungroup()
aki_10d_cr <- aki_10d_cr[order(aki_10d_cr$patient_id,aki_10d_cr$days_since_admission),]
aki_10d_cr <- aki_10d_cr %>% group_by(patient_id) %>% filter(between(time_from_start,0,10)) %>% ungroup()
aki_10d_cr <- aki_10d_cr %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_7day,min_cr_7day_retro)) %>% ungroup()
aki_10d_cr <- aki_10d_cr %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
aki_10d_cr_summ <- aki_10d_cr %>% group_by(severe,time_from_start) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()
aki_10d_cr_timeplot <- ggplot(aki_10d_cr_summ,aes(x=time_from_start,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity")
print(aki_10d_cr_timeplot)