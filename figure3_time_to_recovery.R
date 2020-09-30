get_day <- function(ratio,time_from_peak,target=1.25) {
  index = detect_index(ratio,function(x) x <= 1.25)
  day = time_from_peak[index]
  day
}

labs_cr_recovery <- peak_trend %>% group_by(patient_id) %>% filter(time_from_peak >= 0)
labs_cr_recovery_tmp <- labs_cr_recovery %>% group_by(patient_id) %>% complete(time_from_peak = full_seq(time_from_peak,1)) %>% mutate(ratio = na.fill(ratio,Inf))

time_to_ratio1.25 <- labs_cr_recovery_tmp %>% split(.$patient_id) %>% map(~get_day(.$ratio,.$time_from_peak,target=1.25)) %>% map_df(~data_frame(.x),.id='patient_id')
colnames(time_to_ratio1.25)[2] <- "time_to_ratio1.25"

labs_aki_summ_index <- labs_aki_summ %>% group_by(patient_id) %>% filter(days_since_admission >= 0) %>% filter(days_since_admission == min(days_since_admission))
index_aki_grade <- labs_aki_summ_index %>% select(patient_id,aki_kdigo_final)
aki_index_recovery <- aki_index %>% group_by(patient_id) %>% filter(severe %in% c(2,4,5)) %>% mutate(severe=ifelse(severe==2,0,1))
aki_index_recovery <- merge(aki_index_recovery,time_to_ratio1.25,by="patient_id",all.x=TRUE)
aki_index_recovery <- aki_index_recovery %>% group_by(patient_id) %>% mutate(recover_1.25x = ifelse(is.na(time_to_ratio1.25),0,1))

discharge_day <- demographics %>% select(patient_id,admission_date,last_discharge_date,deceased) %>% group_by(patient_id) %>% mutate(time_to_death_km = as.numeric(as.Date(last_discharge_date)-as.Date(admission_date))) %>% select(patient_id,deceased,time_to_death_km)
aki_index_recovery <- merge(aki_index_recovery,discharge_day,by="patient_id",all.x=TRUE)
aki_index_recovery <- aki_index_recovery %>% group_by(patient_id) %>% mutate(time_to_ratio1.25 = if_else(recover_1.25x == 0,time_to_death_km,time_to_ratio1.25))
aki_index_recovery <- merge(aki_index_recovery,labs_aki_summ_index[,c(1,17)],by="patient_id",all.x=TRUE)

library(survival)
fit_km_recover <- survfit(Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ severe, data=aki_index_recovery)
plot(fit_km_recover,conf.int=TRUE,mark.time=TRUE,col=c("red","blue"),fun="event",xlab="Time to Recovery (25% of Baseline sCr) (days)",ylab="Proportion Recovered")
legend(x=80,y=0.2,legend=c("Non-severe","Severe"),col=c("red","blue"),fill=c("red","blue"))
coxph_recover <- coxph(Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ severe + aki_kdigo_final, data=aki_index_recovery)
summary(coxph_recover) 

# surv_rrt <- Surv(time=km_no_prior_rrt$aki_admit_to_rrt,event=km_no_prior_rrt$rrt)
# fit_rrt <- survfit(surv_rrt ~ km_no_prior_rrt$aki_kdigo_final,data=km_no_prior_rrt)
# plot_rrt <- ggsurvplot(fit_rrt,data=km_no_prior_rrt,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = theme_bw(),fun="event")
# print(plot_rrt)
# fit.coxph <- coxph(surv_rrt ~ agegroup + sex + race + ihd + htn + dm + aki_kdigo_final,data=km_no_prior_rrt)
# forest_coxph <- ggforest(fit.coxph,data=km_no_prior_rrt)
# print(forest_coxph)
# 
# surv_recovery <- Surv(time=km_no_prior_rrt$time_to_ratio1.25,event=km_no_prior_rrt$recover_1.25x)
# fit_recovery <- survfit(surv_recovery ~ km_no_prior_rrt$agegroup,data=km_no_prior_rrt)
# plot_recovery <- ggsurvplot(fit_recovery,data=km_no_prior_rrt,pval=TRUE)
# print(plot_recovery)