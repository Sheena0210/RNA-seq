
#Q5:survival analysis----
library(dplyr)
library(survival)
library(survminer)

#(1)----
table(pheno$sex, pheno$smoking_status)
table(pheno$sex, pheno$EGFR)
table(pheno$EGFR, pheno$tumor_recurrence)
table(pheno$sex, pheno$EGFR, pheno$tumor_recurrence)

#把tumor_recurrence轉成recurrence_event:0,1
pheno <- pheno %>%mutate(recurrence_event = ifelse(tumor_recurrence == "Yes", 1, 0))
table(pheno$tumor_recurrence, useNA = "ifany")
table(pheno$tumor_recurrence, pheno$recurrence_event, useNA = "ifany")
head(pheno)
#建立新欄位( rfs_time)納入追縱到死亡但為沒有recurrence  不然直接跑recurrence會刪掉所有沒有recurrence者

pheno <- pheno %>%mutate(recurrence_event = ifelse(tumor_recurrence == "Yes", 1, 0),
                         rfs_time = ifelse(recurrence_event == 1,days_to_tumor_recurrence,days_to_last_contact_or_death))

#km curve
fit_sex <- survfit(Surv(rfs_time, recurrence_event) ~ sex,data = pheno)

#plot
plot(fit_sex,xlab = "Days to tumor recurrence",ylab = "Recurrence-free survival probability",
     main = "Kaplan-Meier Curve by Sex",lty = 1:2)
legend("bottomleft",legend = levels(as.factor(pheno$sex)),lty = 1:2)

#log rank test
survdiff(Surv(rfs_time, recurrence_event) ~ sex,data = pheno)



#(2)男生:recurrence-free survival是否會因smoking status和EGFR不同而不同？----
#patient_pheno_male
#吸菸狀態;smoking2
#建立新欄位:egfr+smoking status
patient_pheno_male<-patient_pheno_male%>%mutate(egfr_smoking_group = paste(smoking2, EGFR, sep = " + "))
table(patient_pheno_male$egfr_smoking_group)
#因為分六組的話人數太少 所以合併EGFR變成兩類:WT,MUTATANT
patient_pheno_male<- pheno %>%filter(sex == "Male") %>%mutate(EGFR_binary = ifelse(EGFR == "WT","WT","EGFR mutant"),
                                                              smoking_egfr_group = paste(smoking2, EGFR_binary, sep = " + ")
)
table(patient_pheno_male$smoking_egfr_group, patient_pheno_male$recurrence_event)

#併EGFR變成兩類:WT,MUTATANT 及有無吸菸 的km curve
fit_male_group <- survfit(Surv(rfs_time, recurrence_event) ~ smoking_egfr_group,data =patient_pheno_male)

plot(fit_male_group,xlab = "Follow-up Time (days)",
     ylab = "Recurrence-Free Survival Probability",
     main = "Kaplan-Meier Curve for Recurrence-Free Survival by Smoking and EGFR Status in Male Patients",
     lty = 1:length(levels(as.factor(patient_pheno_male$smoking_egfr_group))))

legend(
  "bottomleft",
  legend = levels(as.factor(patient_pheno_male$smoking_egfr_group)),
  lty = 1:length(levels(as.factor(patient_pheno_male$smoking_egfr_group)))
)
logrank_male <-survdiff(
  Surv(rfs_time, recurrence_event) ~ smoking_egfr_group,
  data =patient_pheno_male
)



p_logrank <- 1 - pchisq(
  logrank_male$chisq,
  length(logrank_male$n) - 1
)

p_logrank





#男生中 分成兩類去看RFS(Ever+mutant vs. others)
#新增OTHERS
pheno <- pheno %>%
  mutate(
    recurrence_event = ifelse(tumor_recurrence == "Yes", 1, 0),
    rfs_time = ifelse(
      recurrence_event == 1,
      days_to_tumor_recurrence,
      days_to_last_contact_or_death
    )
  )

# 建立 male subgroup 和分析分組
#男生中 分成兩類去看RFS(Ever+mutant vs. others)
#新增 rfs_group_ever_mutant，不覆蓋原本變數

pheno <- pheno %>%
  mutate(
    recurrence_event = ifelse(tumor_recurrence == "Yes", 1, 0),
    rfs_time = ifelse(
      recurrence_event == 1,
      days_to_tumor_recurrence,
      days_to_last_contact_or_death
    )
  )

# 建立 male subgroup 和分析分組
patient_pheno_male <- pheno %>%
  filter(sex == "Male") %>%
  mutate(
    smoking2 = ifelse(
      smoking_status %in% c("current-smoker", "former-smoker"),
      "Ever smoker",
      "Never smoker"
    ),
    
    EGFR_binary = ifelse(
      EGFR == "WT",
      "WT",
      "EGFR mutant"
    ),
    
    # 新增二分類變數，不覆蓋原本四組
    rfs_group_ever_mutant = ifelse(
      smoking2 == "Ever smoker" & EGFR_binary == "EGFR mutant",
      "Ever smoker + EGFR mutant",
      "Others"
    ),
    
    rfs_group_ever_mutant = factor(
      rfs_group_ever_mutant,
      levels = c("Others", "Ever smoker + EGFR mutant")
    )
  )

# 檢查人數和事件數
table(patient_pheno_male$rfs_group_ever_mutant)
table(patient_pheno_male$rfs_group_ever_mutant, patient_pheno_male$recurrence_event)

fit_male_ever_mutant <- survfit(
  Surv(rfs_time, recurrence_event) ~ rfs_group_ever_mutant,
  data = patient_pheno_male
)

plot(
  fit_male_ever_mutant,
  xlab = "Follow-up Time (days)",
  ylab = "Recurrence-Free Survival Probability",
  main = "Kaplan-Meier Curve for RFS in Male Patients",
  lty = 1:2,
  lwd = 2
)

legend(
  "bottomleft",
  legend = levels(patient_pheno_male$rfs_group_ever_mutant),
  lty = 1:2,
  lwd = 2
)

logrank_male_ever_mutant <- survdiff(
  Surv(rfs_time, recurrence_event) ~ rfs_group_ever_mutant,
  data = patient_pheno_male
)

logrank_male_ever_mutant

p_logrank <- 1 - pchisq(
  logrank_male_ever_mutant$chisq,
  length(logrank_male_ever_mutant$n) - 1
)

p_logrank




#(3)女性病人中，EGFR status 是否與 recurrence-free survival 有關----
#patient_pheno_female
pheno <- pheno %>%mutate(recurrence_event = ifelse(tumor_recurrence == "Yes", 1, 0))
patient_pheno_female <- pheno %>%
  filter(sex == "Female") %>%
  mutate(EGFR = factor(EGFR, levels = c("exon19del", "L858R", "WT"))
  )
table(patient_pheno_female$EGFR, patient_pheno_female$recurrence_event)


#km curve

fit_female_egfr <- survfit(Surv(rfs_time, recurrence_event) ~ EGFR,data = patient_pheno_female)
plot(fit_female_egfr,xlab = "Follow-up Time (days)",
     ylab = "Recurrence-Free Survival Probability",
     main = "In female, Kaplan-Meier Curve for Recurrence-Free Survival by EGFR ",
     lty = 1:length(levels(as.factor(patient_pheno_female$EGFR)))
)

legend(
  "bottomleft",
  legend = levels(as.factor(patient_pheno_female$EGFR)),
  lty = 1:length(levels(as.factor(patient_pheno_female$EGFR)))
)
survdiff(Surv(rfs_time, recurrence_event) ~ EGFR,data = patient_pheno_female)


