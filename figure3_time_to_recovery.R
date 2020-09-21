labs_cr_recovery <- peak_trend %>% group_by(patient_id) %>% filter(time_from_peak >= 0)
labs_cr_recovery_tmp <- labs_cr_recovery %>% group_by(patient_id) %>% complete(time_from_peak = full_seq(time_from_peak,1)) %>% mutate(ratio = na.fill(ratio,Inf))

get_day <- function(ratio,time_from_peak,target=1.25) {
  index = detect_index(ratio,function(x) x <= 1.25)
  day = time_from_peak[index]
  day
}

time_to_ratio1.25 <- labs_cr_recovery_tmp %>% split(.$patient_id) %>% map(~get_day(.$ratio,.$time_from_peak,target=1.25)) %>% map_df(~data_frame(.x),.id='patient_id')
colnames(time_to_ratio1.25)[2] <- "time_to_ratio1.25"

# get_agegroup <- function(age) {
#   agegroup = "0 to 25"
#   if(age >= 26 && age <= 40) {
#     agegroup = "26 to 40"
#   } else if(age > 40 && age <= 60) {
#     agegroup = "41 to 60"
#   } else if (age > 60 && age <= 80) {
#     agegroup = "61 to 80"
#   } else if (age > 80) {
#     agegroup = "80 and above"
#   }
#   agegroup
# }
# 
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