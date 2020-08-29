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
# patient_id	site_id	days_since_admission	value	day_min	day_min_retro	min_cr_7day	min_cr_48h	min_cr_7day_retro	min_cr_48h_retro	min_cr_7d_final	cr_7d	cr_90d	delta_cr	aki_kdigo	aki_kdigo_retro	aki_kdigo_final	akd_7d	akd_90d	severe  time_to_severe	severe_to_aki	severe_before_aki

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
aki_index$covidrx[is.na(aki_index$covid_rx)] <- 0
aki_index <- aki_index %>% group_by(patient_id) %>% mutate(covidrx_grp = ifelse(severe <= 2, ifelse(covid_rx == 0,1,2),ifelse(covid_rx == 0,3,4))) %>% ungroup()
aki_index <- aki_index %>% select(patient_id,peak_cr_time,severe,aki_start,severe_to_aki,covidrx_grp)
# Headers of aki_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki  covidrx_grp

# Uncomment the following line to remove patients who were previously on RRT prior to admission
# aki_index <- aki_index[!(aki_index$patient_id %in% patients_already_rrt),]

# Create a common labs_cr_all table containing the serum Cr values, the severity groupings and anti-viral groupings
labs_cr_aki_tmp <- labs_cr_aki %>% select(patient_id,days_since_admission,value,day_min,min_cr_7day)
labs_cr_nonaki_tmp <- labs_cr_nonaki %>% select(patient_id,days_since_admission,value,day_min,min_cr_7day)
labs_cr_all <- bind_rows(labs_cr_aki_tmp,labs_cr_nonaki_tmp)
labs_cr_all <- merge(labs_cr_all,aki_index,by="patient_id",all.x=TRUE)

# Now, generate a table containing lab values with timepoints calculated from time of peak cr
peak_trend <- labs_cr_all
peak_trend <- peak_trend %>% select(patient_id,severe,covidrx_grp,days_since_admission,peak_cr_time,value,min_cr_7day,min_cr_7day_retro)
# patient_id  severe  covidrx_grp  days_since_admission  peak_cr_time  value min_cr_7day min_cr_7day_retro severe_to_aki

# Calculate the day from peak Cr
peak_trend <- peak_trend %>% group_by(patient_id) %>% mutate(time_from_peak = days_since_admission - peak_cr_time) %>% ungroup()
# Filter this table for Cr values that fall within the 7 day window
peak_trend <- peak_trend %>% group_by(patient_id) %>% filter(between(time_from_peak,0,7)) %>% ungroup()
# Normalise to baseline values used for AKI calculation
peak_trend <- peak_trend %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_7day,min_cr_7day_retro)) %>% ungroup()
peak_trend <- peak_trend %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()

# peak_trend will now be a common table to plot from the selected AKI peak

# =======================================================================================
# Figure 1(a): Cr trends from start of AKI / after peak Cr for severe vs non-severe groups
# =======================================================================================

# Now, derive our first table peak_trend_severe to compare across the different severity groups
peak_trend_severe <- peak_trend %>% select(patient_id,severe,time_from_peak,ratio)
# Headers: patient_id  severe  time_from_peak  ratio

# Calculate mean and SD each for severe and non-severe groups
peak_cr_summ <- peak_trend_severe %>% group_by(severe,time_from_peak) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()

# Plot the graphs
peak_cr_timeplot <- ggplot(peak_cr_summ,aes(x=time_from_peak,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity")
print(peak_cr_timeplot)

# Plot from start of admission to 10 days post-peak AKI (if no peak, then up to 10 days)
adm_to_aki_cr <- labs_cr_all
#adm_to_aki_cr$peak_cr_time[is.na(adm_to_aki_cr$peak_cr_time)] <- 0
adm_to_aki_cr <- adm_to_aki_cr %>% group_by(patient_id) %>% filter(between(days_since_admission,0,peak_cr_time+10)) %>% ungroup()
adm_to_aki_cr <- adm_to_aki_cr[order(adm_to_aki_cr$patient_id,adm_to_aki_cr$days_since_admission),]
adm_to_aki_cr <- adm_to_aki_cr %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_7day,min_cr_7day_retro)) %>% ungroup()
adm_to_aki_cr <- adm_to_aki_cr %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
adm_to_aki_summ <- adm_to_aki_cr %>% group_by(severe,days_since_admission) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()
adm_to_aki_timeplot <- ggplot(adm_to_aki_summ[which(adm_to_aki_summ$days_since_admission <= 10),],aes(x=days_since_admission,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity")
print(adm_to_aki_timeplot)

# Plot from start of AKI to 10 days later 

aki_10d_cr <- labs_cr_all
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
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% select(patient_id,severe,covidrx_grp,days_since_admission,covid_rx_start,peak_cr_time,value,min_cr_7day,min_cr_7day_retro)
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% group_by(patient_id) %>% mutate(severe = ifelse(severe <= 2,0,1)) %>% mutate(covidrx_grp = ifelse((covidrx_grp == 2 || covidrx_grp == 4),1,0)) %>% ungroup()
# Calculate the day from initiation of novel anti-virals
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% group_by(patient_id) %>% mutate(time_from_covidrx = days_since_admission - ifelse(covidrx_grp == 1,covid_rx_start,0)) %>% ungroup()
# Filter this table for Cr values that fall within the 7 day window
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% group_by(patient_id) %>% filter(between(time_from_covidrx,0,7)) %>% ungroup()
# Normalise to baseline values used for AKI calculation
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% group_by(patient_id) %>% mutate(baseline_cr = min(min_cr_7day,min_cr_7day_retro)) %>% ungroup()
cr_from_covidrx_trend <- cr_from_covidrx_trend %>% group_by(patient_id) %>% mutate(ratio = value/baseline_cr) %>% ungroup()
cr_from_covidrx_trend_severe <- cr_from_covidrx_trend %>% select(patient_id,severe,time_from_covidrx,ratio)
# Headers: patient_id  severe (coded as 0/1)  time_from_covidrx  ratio
cr_from_covidrx_summ <- cr_from_covidrx_trend_severe %>% group_by(severe,time_from_covidrx) %>% summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(n())) %>% ungroup()
cr_from_covidrx_timeplot <- ggplot(cr_from_covidrx_summ,aes(x=time_from_covidrx,y=mean_ratio,group=severe))+geom_line(aes(color = factor(severe))) + geom_point(aes(color = factor(severe))) + geom_errorbar(aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=position_dodge(0.05))+ theme(legend.position="right") + labs(x = "Days from COVIDVIRAL Start",y = "Serum Cr/Baseline Cr", color = "Severity")
print(cr_from_covidrx_timeplot)

