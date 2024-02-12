req_lib <- c('epiDisplay','lubridate','moonBook','ztable','survminer','knitr','kableExtra','tidyr','survivalAnalysis','rms','survival','dplyr','ggplot2','VIM')
for (pkg in req_lib) {
  if (!(pkg %in% rownames(installed.packages()))) {install.packages(pkg)}
  else {update.packages(pkg)
    library(pkg, character.only = TRUE)}}

getwd()
setwd("/Users/seochangok/Documents/NSTEMI_CKD")

rm(list=ls())
df <- read.csv("./NSTEMI_CKD_0128.csv", fileEncoding = 'euc-kr') # 최종 정제된 데이터 (n = 11910)
df2 <- data.frame(df)

df$acd12 <- ifelse(df$acd12 == 1 & df$acdfu < 366, 1, 0)
df$cd12 <- ifelse(df$cd12 == 1 & df$cdfu < 366, 1, 0)
df$mi12 <- ifelse(df$mi12 == 1 & df$mifu < 366, 1, 0)
df$cva12<- ifelse(df$cva12 == 1 & df$cvafu < 366, 1, 0)
df$PCI_group <- ifelse(df$PCI_time_2group == 0, 'EIS',
                       ifelse(df$PCI_time_2group == 1, 'DIS', NA))
df$GFR_group <- ifelse(df$GFR_2group == 0, 'PRF',
                       ifelse(df$GFR_2group == 1, 'DRF', NA))
df$GRS_group <- ifelse(!is.na(df$grace), ifelse(df$grace <= 140, 'lowGRS', 'highGRS'), NA)

# PCI_time_2group_early 0: delayed, 1: early
df$PCI_time_2group_early <- 1 - df$PCI_time_2group # HR 구할 때 필요

df_PRF <- subset(df, GFR_2group == 0)
df_DRF <- subset(df, GFR_2group == 1)
df_lowGRS <- subset(df, GRS_group == 'lowGRS')
df_highGRS <- subset(df, GRS_group == 'highGRS')
df_PRF_lowGRS <- subset(df_PRF, GRS_group == 'lowGRS')
df_PRF_highGRS <- subset(df_PRF, GRS_group == 'highGRS')
df_DRF_lowGRS <- subset(df_DRF, GRS_group == 'lowGRS')
df_DRF_highGRS <- subset(df_DRF, GRS_group == 'highGRS')

# 기술통계량
print(ztable(mytable( ~ .-X, data = df)), type="viewer")

print(ztable(mytable(GFR_2group~ .-X-PCI_time_2group, data = df)), type="viewer")
print(ztable(mytable(GFR_2group + PCI_time_2group ~ .-X-diseased_vessels, data = df)), type="viewer")
print(ztable(mytable(GFR_2group + PCI_time_2group ~ acd12 + cd12 + mi12 + cva12, data = df)), type="viewer")
print(ztable(mytable(GFR_2group + GRS_group ~ acd12 + cd12 + mi12 + cva12, data = df)), type="viewer")
print(ztable(mytable(GFR_2group + PCI_time_2group ~ acd12 + cd12 + mi12 + cva12, data = df_lowGRS)), type="viewer")
print(ztable(mytable(GFR_2group + PCI_time_2group ~ acd12 + cd12 + mi12 + cva12, data = df_highGRS)), type="viewer")

# Cox proportional hazard ratio
print(ztable(mytable(GFR_2group + PCI_time_2group ~ acd12 + cd12 + mi12 + cva12 + In_death + bleeding, data = df_complete)), type="viewer")

# ==============================================================================
# Crude & Adjusted Hazard ratio
# ==============================================================================
## All cause death
acd12_PRF <- coxph(Surv(acdfu, acd12) ~ PCI_time_2group_early+grace+Vessel_multi+Revasc+atri_fibril+Hb+male+age+BMI+Current_smoker+HTN+DM+Dyslipidemia+Prev_MI+Prev_PCI+Prev_CABG+Prev_CVA+SBP+Cardiogenic_shock+LVEF, data = df_PRF)
cox.display(acd12_PRF, crude.p.value = T, decimal=2)
ggforest(acd12_PRF, data=df_PRF, main="Cox proportional hazard ratio of all cause death in GFR >= 60")

acd12_DRF <- coxph(Surv(acdfu, acd12) ~ PCI_time_2group_early+grace+Vessel_multi+Revasc+atri_fibril+Hb+male+age+BMI+Current_smoker+HTN+DM+Dyslipidemia+Prev_MI+Prev_PCI+Prev_CABG+Prev_CVA+SBP+Cardiogenic_shock+LVEF, data = df_DRF)
cox.display(acd12_DRF, crude.p.value = T, decimal=2)
ggforest(acd12_DRF, data=df_DRF, main="Cox proportional hazard ratio of all cause death in GFR < 60")

## Cardiac death
cd12_PRF <- coxph(Surv(cdfu, cd12) ~ PCI_time_2group_early+grace+Vessel_multi+Revasc+atri_fibril+Hb+male+age+BMI+Current_smoker+HTN+DM+Dyslipidemia+Prev_MI+Prev_PCI+Prev_CABG+Prev_CVA+SBP+Cardiogenic_shock+LVEF, data = df_PRF)
cox.display(cd12_PRF, crude.p.value = T, decimal=2)
ggforest(cd12_PRF, data=df_PRF, main="Cox proportional hazard ratio of cardiac death in GFR >= 60")

cd12_DRF <- coxph(Surv(cdfu, cd12) ~ PCI_time_2group_early+grace+Vessel_multi+Revasc+atri_fibril+Hb+male+age+BMI+Current_smoker+HTN+DM+Dyslipidemia+Prev_MI+Prev_PCI+Prev_CABG+Prev_CVA+SBP+Cardiogenic_shock+LVEF, data = df_DRF)
cox.display(cd12_DRF, crude.p.value = T, decimal=2)
ggforest(cd12_DRF, data=df_DRF, main="Cox proportional hazard ratio of cardiac death in GFR < 60")

# Non fatal MI
mi12_PRF <- coxph(Surv(mifu, mi12) ~ PCI_time_2group_early+grace+Vessel_multi+Revasc+atri_fibril+Hb+male+age+BMI+Current_smoker+HTN+DM+Dyslipidemia+Prev_MI+Prev_PCI+Prev_CABG+Prev_CVA+SBP+Cardiogenic_shock+LVEF, data = df_PRF)
cox.display(mi12_PRF, crude.p.value = T, decimal=2)
ggforest(mi12_PRF, data=df_PRF, main="Cox proportional hazard ratio of non fatal MI in GFR >= 60")

mi12_DRF <- coxph(Surv(mifu, mi12) ~ PCI_time_2group_early+grace+Vessel_multi+Revasc+atri_fibril+Hb+male+age+BMI+Current_smoker+HTN+DM+Dyslipidemia+Prev_MI+Prev_PCI+Prev_CABG+Prev_CVA+SBP+Cardiogenic_shock+LVEF, data = df_DRF)
cox.display(mi12_DRF, crude.p.value = T, decimal=2)
ggforest(mi12_DRF, data=df_DRF, main="Cox proportional hazard ratio of non fatal MI in GFR < 60")

# Stroke event
cva12_PRF <- coxph(Surv(cvafu, cva12) ~ PCI_time_2group_early+grace+Vessel_multi+Revasc+atri_fibril+Hb+male+age+BMI+Current_smoker+HTN+DM+Dyslipidemia+Prev_MI+Prev_PCI+Prev_CABG+Prev_CVA+SBP+Cardiogenic_shock+LVEF, data = df_PRF)
cox.display(cva12_PRF, crude.p.value = T, decimal=2)
ggforest(cva12_PRF, data=df_PRF, main="Cox proportional hazard ratio of stroke event in GFR >= 60")

cva12_DRF <- coxph(Surv(cvafu, cva12) ~ PCI_time_2group_early+grace+Vessel_multi+Revasc+atri_fibril+Hb+male+age+BMI+Current_smoker+HTN+DM+Dyslipidemia+Prev_MI+Prev_PCI+Prev_CABG+Prev_CVA+SBP+Cardiogenic_shock+LVEF, data = df_DRF)
cox.display(cva12_DRF, crude.p.value = T, decimal=2)
ggforest(cva12_DRF, data=df_DRF, main="Cox proportional hazard ratio of stroke event in GFR < 60")

# ==============================================================================
# KM plot
# ============================================================================== 
## Log-rank : survdiff(Surv(time,status)~sex,data=lung)
acd12 <- survfit(Surv(acdfu, acd12) ~ GFR_group + PCI_group, data = df)
acd12_km <- ggsurvplot(acd12, title = "All cause Death", legend = c(0.15, 0.2), legend.title = "", 
                       legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                       censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                       conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                       break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                       xlab = "Time(days)", ylab = "Probability",
                       font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cd12 <- survfit(Surv(cdfu, cd12) ~ GFR_group + PCI_group, data = df)
cd12_km <- ggsurvplot(cd12, title = "Cardiac Death", legend = c(0.15, 0.2), legend.title = "", 
                      legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                      censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                      conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                      break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                      xlab = "Time(days)", ylab = "Probability",
                      font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
mi12 <- survfit(Surv(mifu, mi12) ~ GFR_group + PCI_group, data = df)
mi12_km <- ggsurvplot(mi12, title = "Non-fatal MI", legend = c(0.15, 0.2), legend.title = "", 
                      legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                      censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                      conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                      break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                      xlab = "Time(days)", ylab = "Probability",
                      font.x = c(15, "bold"),font.y = c(15, "bold"), font.tickslab = c(15))
cva12 <- survfit(Surv(cvafu, cva12) ~ GFR_group + PCI_group, data = df)
cva12_km <- ggsurvplot(cva12, title = "Stroke", legend = c(0.15, 0.2), legend.title = "", 
                       legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                       censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                       conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                       break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                       xlab = "Time(days)", ylab = "Probability",
                       font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))

acd12_PRF <- survfit(Surv(acdfu, acd12) ~ PCI_group, data = df_PRF)
acd12_PRF_km <- ggsurvplot(acd12_PRF, title = "All cause Death (eGFR >= 60)", legend = c(0.15, 0.4), legend.title = "",
                           legend.labs = c('DIS', 'EIS'),
                           censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.2, fontsize = 5,
                           conf.int = F, pval = T, pval.size = 7, pval.coord = c(1, 0.83),
                           break.time.by = 90, xlim = c(0, 365), ylim = c(0.8, 1.0),
                           xlab = "Time(days)", ylab = "Probability",
                           font.main = c(20, "bold"), font.x = c(18, "bold"), font.xlab = c(18), font.y = c(18, "bold"), font.ylab = c(18),
                           font.legend = c(16, "bold"), font.tickslab = c(18))
acd12_PRF_km$table <- acd12_PRF_km$table + theme_cleantable()
acd12_PRF_km$plot <- acd12_PRF_km$plot + annotate("text", x = 1, y = 0.81, size = 6, hjust = 0,
                                                  label = "Adjusted hazard ratio 0.95 (95% CI 0.70 - 1.29)")
acd12_PRF_km
cd12_PRF <- survfit(Surv(cdfu, cd12) ~ PCI_group, data = df_PRF)
cd12_PRF_km <- ggsurvplot(cd12_PRF, title = "Cardiac Death (eGFR >= 60)", legend = c(0.15, 0.2), legend.title = "", 
                          legend.labs = c('DIS', 'EIS'),
                          censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                          conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                          break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                          xlab = "Time(days)", ylab = "Probability",
                          font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
mi12_PRF <- survfit(Surv(mifu, mi12) ~ PCI_group, data = df_PRF)
mi12_PRF_km <- ggsurvplot(mi12_PRF, title = "Non-fatal MI (eGFR >= 60)", legend = c(0.15, 0.2), legend.title = "", 
                          legend.labs = c('DIS', 'EIS'),
                          censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                          conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                          break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                          xlab = "Time(days)", ylab = "Probability",
                          font.x = c(15, "bold"),font.y = c(15, "bold"), font.tickslab = c(15))
cva12_PRF <- survfit(Surv(cvafu, cva12) ~ PCI_group, data = df_PRF)
cva12_PRF_km <- ggsurvplot(cva12_PRF, title = "Stroke (eGFR >= 60)", legend = c(0.15, 0.2), legend.title = "", 
                           legend.labs = c('DIS', 'EIS'),
                           censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                           conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                           break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                           xlab = "Time(days)", ylab = "Probability",
                           font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
acd12_DRF <- survfit(Surv(acdfu, acd12) ~ PCI_group, data = df_DRF)
acd12_DRF_km <- ggsurvplot(acd12_DRF, title = "All cause Death (eGFR < 60)", legend = c(0.15, 0.4), legend.title = "",
                           legend.labs = c('DIS', 'EIS'),
                           censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.2, fontsize = 5,
                           conf.int = F, pval = T, pval.size = 7, pval.coord = c(1, 0.83),
                           break.time.by = 90, xlim = c(0, 365), ylim = c(0.8, 1.0),
                           xlab = "Time(days)", ylab = "Probability",
                           font.main = c(20, "bold"), font.x = c(18, "bold"), font.xlab = c(18), font.y = c(18, "bold"), font.ylab = c(18),
                           font.legend = c(16, "bold"), font.tickslab = c(18))
acd12_DRF_km$table <- acd12_DRF_km$table + theme_cleantable()
acd12_DRF_km$plot <- acd12_DRF_km$plot + annotate("text", x = 1, y = 0.81, size = 6, hjust = 0,
                                                  label = "Adjusted hazard ratio 1.43 (95% CI 1.12 - 1.83)")
acd12_DRF_km
cd12_DRF <- survfit(Surv(cdfu, cd12) ~ PCI_group, data = df_DRF)
cd12_DRF_km <- ggsurvplot(cd12_DRF, title = "Cardiac Death (eGFR < 60)", legend = c(0.15, 0.2), legend.title = "", 
                          legend.labs = c('DIS', 'EIS'),
                          censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                          conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                          break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                          xlab = "Time(days)", ylab = "Probability",
                          font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
mi12_DRF <- survfit(Surv(mifu, mi12) ~ PCI_group, data = df_DRF)
mi12_DRF_km <- ggsurvplot(mi12_DRF, title = "Non-fatal MI (eGFR < 60)", legend = c(0.15, 0.2), legend.title = "", 
                          legend.labs = c('DIS', 'EIS'),
                          censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                          conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                          break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                          xlab = "Time(days)", ylab = "Probability",
                          font.x = c(15, "bold"),font.y = c(15, "bold"), font.tickslab = c(15))
cva12_DRF <- survfit(Surv(cvafu, cva12) ~ PCI_group, data = df_DRF)
cva12_DRF_km <- ggsurvplot(cva12_DRF, title = "Stroke (eGFR < 60)", legend = c(0.15, 0.2), legend.title = "", 
                           legend.labs = c('DIS', 'EIS'),
                           censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                           conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                           break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                           xlab = "Time(days)", ylab = "Probability",
                           font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))

acd12_lowGRS <- survfit(Surv(acdfu, acd12) ~ GFR_group + PCI_group, data = df_lowGRS)
acd12_lowGRS_km <- ggsurvplot(acd12_lowGRS, title = "All cause Death (lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                              legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                              censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                              conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                              break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                              xlab = "Time(days)", ylab = "Probability",
                              font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
acd12_highGRS <- survfit(Surv(acdfu, acd12) ~ GFR_group + PCI_group, data = df_highGRS)
acd12_highGRS_km <- ggsurvplot(acd12_highGRS, title = "All cause Death (highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                               legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                               censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                               conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                               break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                               xlab = "Time(days)", ylab = "Probability",
                               font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cd12_lowGRS <- survfit(Surv(acdfu, cd12) ~ GFR_group + PCI_group, data = df_lowGRS)
cd12_lowGRS_km <- ggsurvplot(cd12_lowGRS, title = "Cardiac Death (lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                             legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                             censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                             conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                             break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                             xlab = "Time(days)", ylab = "Probability",
                             font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cd12_highGRS <- survfit(Surv(acdfu, cd12) ~ GFR_group + PCI_group, data = df_highGRS)
cd12_highGRS_km <- ggsurvplot(cd12_highGRS, title = "Cardiac Death (highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                              legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                              censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                              conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                              break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                              xlab = "Time(days)", ylab = "Probability",
                              font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
mi12_lowGRS <- survfit(Surv(acdfu, mi12) ~ GFR_group + PCI_group, data = df_lowGRS)
mi12_lowGRS_km <- ggsurvplot(mi12_lowGRS, title = "Non-fatal MI (lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                             legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                             censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                             conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                             break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                             xlab = "Time(days)", ylab = "Probability",
                             font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
mi12_highGRS <- survfit(Surv(acdfu, mi12) ~ GFR_group + PCI_group, data = df_highGRS)
mi12_highGRS_km <- ggsurvplot(mi12_highGRS, title = "Non-fatal MI (highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                              legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                              censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                              conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                              break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                              xlab = "Time(days)", ylab = "Probability",
                              font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cva12_lowGRS <- survfit(Surv(acdfu, cva12) ~ GFR_group + PCI_group, data = df_lowGRS)
cva12_lowGRS_km <- ggsurvplot(cva12_lowGRS, title = "Stroke (lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                              legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                              censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                              conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                              break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                              xlab = "Time(days)", ylab = "Probability",
                              font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cva12_highGRS <- survfit(Surv(acdfu, cva12) ~ GFR_group + PCI_group, data = df_highGRS)
cva12_highGRS_km <- ggsurvplot(cva12_highGRS, title = "Stroke (highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                               legend.labs = c('D + D', 'D + E', 'P + D', 'P + E'),
                               censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                               conf.int = F, pval = F, pval.size = 4, pval.coord = c(1, 0.8),
                               break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                               xlab = "Time(days)", ylab = "Probability",
                               font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))

acd12_PRF_lowGRS <- survfit(Surv(acdfu, acd12) ~ PCI_group, data = df_PRF_lowGRS)
acd12_PRF_lowGRS_km <- ggsurvplot(acd12_PRF_lowGRS, title = "All cause Death (eGFR >= 60 & lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                  legend.labs = c('DIS', 'EIS'),
                                  censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                  conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                  break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                  xlab = "Time(days)", ylab = "Probability",
                                  font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
acd12_PRF_highGRS <- survfit(Surv(acdfu, acd12) ~ PCI_group, data = df_PRF_highGRS)
acd12_PRF_highGRS_km <- ggsurvplot(acd12_PRF_highGRS, title = "All cause Death (eGFR >= 60 & highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                   legend.labs = c('DIS', 'EIS'),
                                   censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                   conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                   break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                   xlab = "Time(days)", ylab = "Probability",
                                   font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cd12_PRF_lowGRS <- survfit(Surv(acdfu, cd12) ~ PCI_group, data = df_PRF_lowGRS)
cd12_PRF_lowGRS_km <- ggsurvplot(cd12_PRF_lowGRS, title = "Cardiac Death (eGFR >= 60 & lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                 legend.labs = c('DIS', 'EIS'),
                                 censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                 conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                 break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                 xlab = "Time(days)", ylab = "Probability",
                                 font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cd12_PRF_highGRS <- survfit(Surv(acdfu, cd12) ~ PCI_group, data = df_PRF_highGRS)
cd12_PRF_highGRS_km <- ggsurvplot(cd12_PRF_highGRS, title = "Cardiac Death (eGFR >= 60 & highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                  legend.labs = c('DIS', 'EIS'),
                                  censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                  conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                  break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                  xlab = "Time(days)", ylab = "Probability",
                                  font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
mi12_PRF_lowGRS <- survfit(Surv(acdfu, mi12) ~ PCI_group, data = df_PRF_lowGRS)
mi12_PRF_lowGRS_km <- ggsurvplot(mi12_PRF_lowGRS, title = "Non-fatal MI (eGFR >= 60 & lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                 legend.labs = c('DIS', 'EIS'),
                                 censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                 conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                 break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                 xlab = "Time(days)", ylab = "Probability",
                                 font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
mi12_PRF_highGRS <- survfit(Surv(acdfu, mi12) ~ PCI_group, data = df_PRF_highGRS)
mi12_PRF_highGRS_km <- ggsurvplot(mi12_PRF_highGRS, title = "Non-fatal MI (eGFR >= 60 & highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                  legend.labs = c('DIS', 'EIS'),
                                  censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                  conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                  break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                  xlab = "Time(days)", ylab = "Probability",
                                  font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cva12_PRF_lowGRS <- survfit(Surv(acdfu, cva12) ~ PCI_group, data = df_PRF_lowGRS)
cva12_PRF_lowGRS_km <- ggsurvplot(cva12_PRF_lowGRS, title = "Stroke (eGFR >= 60 & lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                  legend.labs = c('DIS', 'EIS'),
                                  censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                  conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                  break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                  xlab = "Time(days)", ylab = "Probability",
                                  font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cva12_PRF_highGRS <- survfit(Surv(acdfu, cva12) ~ PCI_group, data = df_PRF_highGRS)
cva12_PRF_highGRS_km <- ggsurvplot(cva12_PRF_highGRS, title = "Stroke (eGFR >= 60 & highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                   legend.labs = c('DIS', 'EIS'),
                                   censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                   conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                   break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                   xlab = "Time(days)", ylab = "Probability",
                                   font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))

acd12_DRF_lowGRS <- survfit(Surv(acdfu, acd12) ~ PCI_group, data = df_DRF_lowGRS)
acd12_DRF_lowGRS_km <- ggsurvplot(acd12_DRF_lowGRS, title = "All cause Death (eGFR < 60 & lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                  legend.labs = c('DIS', 'EIS'),
                                  censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                  conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                  break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                  xlab = "Time(days)", ylab = "Probability",
                                  font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
acd12_DRF_highGRS <- survfit(Surv(acdfu, acd12) ~ PCI_group, data = df_DRF_highGRS)
acd12_DRF_highGRS_km <- ggsurvplot(acd12_DRF_highGRS, title = "All cause Death (eGFR < 60 & highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                   legend.labs = c('DIS', 'EIS'),
                                   censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                   conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                   break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                   xlab = "Time(days)", ylab = "Probability",
                                   font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cd12_DRF_lowGRS <- survfit(Surv(acdfu, cd12) ~ PCI_group, data = df_DRF_lowGRS)
cd12_DRF_lowGRS_km <- ggsurvplot(cd12_DRF_lowGRS, title = "Cardiac Death (eGFR < 60 & lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                 legend.labs = c('DIS', 'EIS'),
                                 censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                 conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                 break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                 xlab = "Time(days)", ylab = "Probability",
                                 font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cd12_DRF_highGRS <- survfit(Surv(acdfu, cd12) ~ PCI_group, data = df_DRF_highGRS)
cd12_DRF_highGRS_km <- ggsurvplot(cd12_DRF_highGRS, title = "Cardiac Death (eGFR < 60 & highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                  legend.labs = c('DIS', 'EIS'),
                                  censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                  conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                  break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                  xlab = "Time(days)", ylab = "Probability",
                                  font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
mi12_DRF_lowGRS <- survfit(Surv(acdfu, mi12) ~ PCI_group, data = df_DRF_lowGRS)
mi12_DRF_lowGRS_km <- ggsurvplot(mi12_DRF_lowGRS, title = "Non-fatal MI (eGFR < 60 & lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                 legend.labs = c('DIS', 'EIS'),
                                 censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                 conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                 break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                 xlab = "Time(days)", ylab = "Probability",
                                 font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
mi12_DRF_highGRS <- survfit(Surv(acdfu, mi12) ~ PCI_group, data = df_DRF_highGRS)
mi12_DRF_highGRS_km <- ggsurvplot(mi12_DRF_highGRS, title = "Non-fatal MI (eGFR < 60 & highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                  legend.labs = c('DIS', 'EIS'),
                                  censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                  conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                  break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                  xlab = "Time(days)", ylab = "Probability",
                                  font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cva12_DRF_lowGRS <- survfit(Surv(acdfu, cva12) ~ PCI_group, data = df_DRF_lowGRS)
cva12_DRF_lowGRS_km <- ggsurvplot(cva12_DRF_lowGRS, title = "Stroke (eGFR < 60 & lowGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                  legend.labs = c('DIS', 'EIS'),
                                  censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                  conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                  break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                  xlab = "Time(days)", ylab = "Probability",
                                  font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))
cva12_DRF_highGRS <- survfit(Surv(acdfu, cva12) ~ PCI_group, data = df_DRF_highGRS)
cva12_DRF_highGRS_km <- ggsurvplot(cva12_DRF_highGRS, title = "Stroke (eGFR < 60 & highGRS)", legend = c(0.15, 0.2), legend.title = "", 
                                   legend.labs = c('DIS', 'EIS'),
                                   censor = F, risk.table.title = "No. at Risk", risk.table = T, risk.table.height = 0.3,
                                   conf.int = F, pval = T, pval.size = 6, pval.coord = c(1, 0.9),
                                   break.time.by = 120, xlim = c(0, 365), ylim = c(0.8, 1.0),
                                   xlab = "Time(days)", ylab = "Probability",
                                   font.x = c(15, "bold"), font.y = c(15, "bold"), font.tickslab = c(15))

# ============================================================================== 
# KM plot print
# ============================================================================== 
acd12_km
summary(acd12, times = 365)$surv
cd12_km
summary(cd12, times = 365)$surv
mi12_km
summary(mi12, times = 365)$surv
cva12_km
summary(cva12, times = 365)$surv
acd12_PRF_km
summary(acd12_PRF, times = 365)$surv
cd12_PRF_km
summary(cd12_PRF, times = 365)$surv
mi12_PRF_km
summary(mi12_PRF, times = 365)$surv
cva12_PRF_km
summary(cva12_PRF, times = 365)$surv
acd12_DRF_km
summary(acd12_DRF, times = 365)$surv
cd12_DRF_km
summary(cd12_DRF, times = 365)$surv
mi12_DRF_km
summary(mi12_DRF, times = 365)$surv
cva12_DRF_km
summary(cva12_DRF, times = 365)$surv
acd12_lowGRS_km
summary(acd12_lowGRS, times = 365)$surv
acd12_highGRS_km
summary(acd12_highGRS, times = 365)$surv
cd12_lowGRS_km
summary(cd12_lowGRS, times = 365)$surv
cd12_highGRS_km
summary(cd12_highGRS, times = 365)$surv
mi12_lowGRS_km
summary(mi12_lowGRS, times = 365)$surv
mi12_highGRS_km
summary(mi12_highGRS, times = 365)$surv
cva12_lowGRS_km
summary(cva12_lowGRS, times = 365)$surv
cva12_highGRS_km
summary(cva12_highGRS, times = 365)$surv
acd12_PRF_lowGRS_km
summary(acd12_PRF_lowGRS, times = 365)$surv
acd12_PRF_highGRS_km
summary(acd12_PRF_highGRS, times = 365)$surv
cd12_PRF_lowGRS_km
summary(cd12_PRF_lowGRS, times = 365)$surv
cd12_PRF_highGRS_km
summary(cd12_PRF_highGRS, times = 365)$surv
mi12_PRF_lowGRS_km
summary(mi12_PRF_lowGRS, times = 365)$surv
mi12_PRF_highGRS_km
summary(mi12_PRF_highGRS, times = 365)$surv
cva12_PRF_lowGRS_km
summary(cva12_PRF_lowGRS, times = 365)$surv
cva12_PRF_highGRS_km
summary(cva12_PRF_highGRS, times = 365)$surv
acd12_DRF_lowGRS_km
summary(acd12_DRF_lowGRS, times = 365)$surv
acd12_DRF_highGRS_km
summary(acd12_DRF_highGRS, times = 365)$surv
cd12_DRF_lowGRS_km
summary(cd12_DRF_lowGRS, times = 365)$surv
cd12_DRF_highGRS_km
summary(cd12_DRF_highGRS, times = 365)$surv
mi12_DRF_lowGRS_km
summary(mi12_DRF_lowGRS, times = 365)$surv
mi12_DRF_highGRS_km
summary(mi12_DRF_highGRS, times = 365)$surv
cva12_DRF_lowGRS_km
summary(cva12_DRF_lowGRS, times = 365)$surv
cva12_DRF_highGRS_km
summary(cva12_DRF_highGRS, times = 365)$surv

# ==============================================================================
# Range GFR model
# ==============================================================================
tmp <- data.frame(df)
lst_acd <- list()
lst_cd <- list()
lst_mi <- list()
lst_cva <- list()

seq_min <- c('<30', '30<= <45', '45<= <60', '60<=')
seq_max <- c('30<= <45', '45<= <60', '60<= <200', '200<=')
seq_min_int <- c(0,30,45,60)
seq_max_int <- c(30,45,60,200)

for (i in 1:length(seq_min)){
  if (i < 60){
    #print(paste('if',i,seq_min_int[i],seq_max_int[i]))
    tmp_sub = subset(tmp, GFR >= seq_min_int[i] & GFR < seq_max_int[i])
  } else {
    print(paste('else',i))
    tmp_sub = subset(tmp, GFR >= seq_min_int[i])}
  
  fit_acd <- coxph(Surv(acdfu, acd12) ~ PCI_time_2group_early, data = tmp_sub)
  fit_cd <- coxph(Surv(cdfu, cd12) ~ PCI_time_2group_early, data = tmp_sub)
  fit_mi <- coxph(Surv(mifu, mi12) ~ PCI_time_2group_early, data = tmp_sub)
  fit_cva <- coxph(Surv(cvafu, cva12) ~ PCI_time_2group_early, data = tmp_sub)
  fit_in_death <- coxph(Surv(cvafu, cva12) ~ PCI_time_2group_early, data = tmp_sub)
  summ_acd <- summary(fit_acd)
  summ_cd <- summary(fit_cd)
  summ_mi <- summary(fit_mi)
  summ_cva <- summary(fit_cva)
  res_acd <- data.frame(GFR = seq_min[i],
                        Variable = names(fit_acd$coefficients),
                        HR = round(summ_acd$conf.int[, 1], 3),
                        lwr.95 = round(summ_acd$conf.int[, 3], 3),
                        upr.95 = round(summ_acd$conf.int[, 4], 3),
                        p.val = round(summ_acd$coefficients[, 5], 3))
  res_cd <- data.frame(GFR=seq_min[i],Variable=names(fit_cd$coefficients),HR=round(summ_cd$conf.int[, 1], 3),lwr.95=round(summ_cd$conf.int[, 3], 3),upr.95=round(summ_cd$conf.int[, 4], 3),p.val=round(summ_cd$coefficients[, 5], 3))
  res_mi <- data.frame(GFR=seq_min[i],Variable=names(fit_mi$coefficients),HR=round(summ_mi$conf.int[, 1], 3),lwr.95=round(summ_mi$conf.int[, 3], 3),upr.95=round(summ_mi$conf.int[, 4], 3),p.val=round(summ_mi$coefficients[, 5], 3))
  res_cva <- data.frame(GFR=seq_min[i],Variable=names(fit_cva$coefficients),HR=round(summ_cva$conf.int[, 1], 3),lwr.95=round(summ_cva$conf.int[, 3], 3),upr.95=round(summ_cva$conf.int[, 4], 3),p.val=round(summ_cva$coefficients[, 5], 3))
  lst_acd[[i]] <- res_acd
  lst_cd[[i]] <- res_cd
  lst_mi[[i]] <- res_mi
  lst_cva[[i]] <- res_cva}
tmp_sub90 = subset(tmp, GFR >= 60)
tmp_sub60 = subset(tmp, GFR >= 45 & GFR < 60)
tmp_sub45 = subset(tmp, GFR >= 30 & GFR < 45)
tmp_sub30 = subset(tmp, GFR < 30)
print(ztable(mytable(PCI_time_2group ~ acd12 + cd12 + mi12 + cva12, data = tmp_sub90)), type = 'viewer')
print(ztable(mytable(PCI_time_2group ~ acd12 + cd12 + mi12 + cva12, data = tmp_sub60)), type = 'viewer')
print(ztable(mytable(PCI_time_2group ~ acd12 + cd12 + mi12 + cva12, data = tmp_sub45)), type = 'viewer')
print(ztable(mytable(PCI_time_2group ~ acd12 + cd12 + mi12 + cva12, data = tmp_sub30)), type = 'viewer')
res_acd <- do.call('rbind', lst_acd)
res_cd <- do.call('rbind', lst_cd)
res_mi <- do.call('rbind', lst_mi)
res_cva <- do.call('rbind', lst_cva)

# print(ztable(res_acd, digits = 3), type="viewer")
# print(ztable(res_cd, digits = 3), type="viewer")
# print(ztable(res_mi, digits = 3), type="viewer")
# print(ztable(res_cva, digits = 3), type="viewer")

## graph
# size는 1200*900
ggplot(res_acd, aes(x = HR, y = GFR))+
  xlim(0.1, 3.3)+
  theme_bw()+
  ggtitle("All cause death")+
  xlab("Hazard Ratio")+ylab("eGFR (mL/min/1.73m^2)")+
  theme(title=element_text(size=20, face="bold"), axis.title=element_text(size=18, face="bold"), axis.text.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=14, face="bold"))+
  geom_errorbar(aes(xmin = lwr.95, xmax = upr.95), width=0.2, size=1)+
  geom_point(size = 3)+
  geom_vline(xintercept = 1, col='red', linewidth = 1.5, linetype = 1)+
  scale_y_discrete(limits = seq_min)
res_acd
ggplot(res_cd, aes(x = HR, y = GFR))+
  xlim(0.1, 3.3)+
  theme_bw()+
  ggtitle("Cardiac death")+
  xlab("Hazard Ratio")+ylab("eGFR (mL/min/1.73m^2)")+
  theme(title=element_text(size=20, face="bold"), axis.title=element_text(size=18, face="bold"), axis.text.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=14, face="bold"))+
  geom_errorbar(aes(xmin = lwr.95, xmax = upr.95), width=0.2, size=1)+
  geom_point(size = 3)+
  geom_vline(xintercept = 1, col='red', linewidth = 1.5, linetype = 1)+
  scale_y_discrete(limits = seq_min)
res_cd
ggplot(res_mi, aes(x = HR, y = GFR))+
  xlim(0.1, 3.3)+
  theme_bw()+
  ggtitle("MI")+
  xlab("Hazard Ratio")+ylab("eGFR (mL/min/1.73m^2)")+
  theme(title=element_text(size=20, face="bold"), axis.title=element_text(size=18, face="bold"), axis.text.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=14, face="bold"))+
  geom_errorbar(aes(xmin = lwr.95, xmax = upr.95), width=0.2, size=1)+
  geom_point(size = 3)+
  geom_vline(xintercept = 1, col='red', linewidth = 1.5, linetype = 1)+
  scale_y_discrete(limits = seq_min)
res_mi
ggplot(res_cva, aes(x = HR, y = GFR))+
  xlim(0.1, 3.3)+
  theme_bw()+
  ggtitle("Stroke")+
  xlab("Hazard Ratio")+ylab("eGFR (mL/min/1.73m^2)")+
  theme(title=element_text(size=20, face="bold"), axis.title=element_text(size=18, face="bold"), axis.text.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=14, face="bold"))+
  geom_errorbar(aes(xmin = lwr.95, xmax = upr.95), width=0.2, size=1)+
  geom_point(size = 3)+
  geom_vline(xintercept = 1, col='red', linewidth = 1.5, linetype = 1)+
  scale_y_discrete(limits = seq_min)
res_cva

# ==============================================================================
# interaction term
# ==============================================================================
dd <- datadist(df)
options(datadist='dd')
f_uni <- cph(Surv(acdfu, acd12) ~ PCI_time_2group_early * GFR, data=df)
pred_uni <- Predict(f_uni, GFR, fun=exp)
ggplot(pred_uni) +
  geom_hline(yintercept = 1, col='red') +
  labs(y = "Hazard Ratio") +
  theme_bw()

# GFR RC spline
dd <- datadist(tmp)
options(datadist='dd')
f_uni <- cph(Surv(acdfu, acd12) ~ rcs(GFR), data=tmp)
pred_uni <- Predict(f_uni, GFR, fun=exp)
ggplot(pred_uni) +
  geom_hline(yintercept = 1, col='red') +
  labs(y = "Hazard Ratio") +
  theme_bw()