#Question:
#保留男生n=27
patient_pheno <- pheno_info %>%distinct(ID, .keep_all = TRUE)%>%filter(sex=='Male')
#抓回男生的男生的資料
#male_info是一個人有兩筆資料
male_info <- pheno_info %>%filter(ID %in% patient_pheno$ID)
table(male_info$sample_type)

#pre-proccess:
#smoke status->smoking2:never vs.ever----
#never vs.ever :18:9
#male_patient：一個病人只保留一列
male_patient <- pheno_info %>%distinct(ID, .keep_all = TRUE) %>%filter(sex == "Male")
male_patient <- male_patient %>%mutate(smoking2 = case_when(
      smoking_status == "never-smoker" ~ "Never",
      smoking_status %in% c("current-smoker", "former-smoker") ~ "Ever",
      TRUE ~ NA_character_
    )
  )
table(male_patient$smoking2)
table(male_patient$EGFR)

Q1 <- table(male_patient$EGFR,male_patient$smoking2)
Q1

fisher_Q1 <- fisher.test(Q1)
fisher_Q1


Q1_table <- prop.table(Q1 , margin = 1) * 100
Q1_table 
round(Q1_table , 1)

Q1 <- as.data.frame.matrix(Q1)
Q1_table <- as.data.frame.matrix(round(Q1_table, 1))

Q1_final <- data.frame(
  EGFR = rownames(Q1),
  Never = paste0(Q1$Never, " (", Q1_table$Never, "%)"),
  Ever = paste0(Q1$Ever, " (", Q1_table$Ever, "%)"),
  Total = rowSums(Q1),
  row.names = NULL
)
Q1_final






#Q1:在男生中，EGFR三個類別是否和smoking status有關？(臨床/病人特徵之間的關係)----
#樣本數太少不用卡方>fisher exact test
Q1 <- table(male_info $EGFR,male_info $smoking2)
fisher.test(Q1)


#Q2:All gender中，EGFR三類對Tumor-Normalexpression change的影響----
#EGFR 三類之間的 Tumor-Normal expression difference 是否整體有差異？
#利用Benjamini-Hochberg多重比較校正後沒有任何基因達到adj pvalue < 0.05
#EGFR類別對 genome-wide Tumor-Normal expression difference 沒有達到統計顯著
all(colnames(logCPM_after) == pheno_info$id)
all_patient_ids <- unique(pheno_info$ID)
diff_mat_all <- sapply(all_patient_ids, function(pid) {tumor_sample <- pheno_info$id[pheno_info$ID == pid & pheno_info$sample_type == "Tumor"]
normal_sample <- pheno_info$id[pheno_info$ID == pid & pheno_info$sample_type == "Normal"]
logCPM_after[, tumor_sample] - logCPM_after[, normal_sample]})

#每一位病人的Tumor-Normal expression difference
#一個人有tumor nomal兩列
colnames(diff_mat_all) <- all_patient_ids
dim(diff_mat_all)
diff_mat_all[1:5, 1:5]

#用diff轉成一個病人一列
patient_pheno_all <- pheno_info %>%distinct(ID, .keep_all = TRUE) %>%filter(ID %in% colnames(diff_mat_all)) %>%arrange(match(ID, colnames(diff_mat_all)))
dim(patient_pheno_all)
all(patient_pheno_all$ID == colnames(diff_mat_all)) #檢查順序是否一致

table(patient_pheno_all$EGFR)
patient_pheno_all$EGFR <- factor(patient_pheno_all$EGFR)
levels(patient_pheno_all$EGFR)

#設定ref:WT(dummy)
patient_pheno_all$EGFR <- relevel(patient_pheno_all$EGFR,ref = 'WT')
design_all_egfr <- model.matrix(~ EGFR,data = patient_pheno_all)
colnames(design_all_egfr)

#一次看很多基因用limma overall f test
#limma 則會把所有 genes 的 variance pattern 一起考慮，讓統計結果比較穩
fit_all_egfr <- lmFit(diff_mat_all, design_all_egfr)
fit_all_egfr <- eBayes(fit_all_egfr)

#EGFR 三類之間的 Tumor-Normal expression difference 是否整體有差異？
#利用Benjamini-Hochberg多重比較校正後沒有任何基因達到adj pvalue < 0.05
#EGFR類別對 genome-wide Tumor-Normal expression difference 沒有達到統計顯著
coef_egfr_all <- grep("^EGFR", colnames(design_all_egfr))
coef_egfr_all
colnames(design_all_egfr)[coef_egfr_all]

deg_all_egfr_overall <- topTable(fit_all_egfr,coef = coef_egfr_all,number = Inf,adjust.method = "BH") %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%left_join(gene_info, by = "ensembl_gene_id")
head(deg_all_egfr_overall)
#定義GENE狀態 如果有顯著就是associated
deg_all_egfr_overall <- deg_all_egfr_overall %>%
  mutate(EGFR_status = case_when(adj.P.Val < 0.05 ~ "EGFR-associated",TRUE ~ "Not significant"))

table(deg_all_egfr_overall$EGFR_status)

write.csv(deg_all_egfr_overall,"All_patients_EGFR3_Tumor_Normal_difference_limma.csv",row.names = FALSE)



#全部人：Top gene boxplot
#top gene:"ENSG00000226074"->"PRSS44"
top_gene_id_all <- deg_all_egfr_overall$ensembl_gene_id[1]
top_gene_name_all <- deg_all_egfr_overall$external_gene_name[1]
top_gene_id_all
top_gene_name_all

plot_df_all<- patient_pheno_all %>%
  mutate(gene_diff = as.numeric(diff_mat_all[top_gene_id, ])
  )

ggplot(plot_df_all, aes(x = EGFR, y = gene_diff, fill = EGFR)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  labs(
    title = paste0("Tumor-Normal expression difference of ", top_gene_name_all, " by EGFR category"),
    x = "EGFR category",
    y = "Tumor - Normal logCPM difference"
  )

#f test
deg_all_egfr_overall <- deg_all_egfr_overall %>%
  mutate(
    EGFR_status = case_when(
      adj.P.Val < 0.05 ~ "EGFR-associated",
      TRUE ~ "Not significant"
    )
  )
ggplot(deg_all_egfr_overall, aes(x = F, y = -log10(adj.P.Val), color = EGFR_status)) +
  geom_point(alpha = 0.6, size = 1.2) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "EGFR three-category effect on Tumor-Normal expression difference",
    x = "F-statistic",
    y = "-log10 adjusted p-value"
  )

#heatmap
library(pheatmap)

top_egfr_genes_all <- deg_all_egfr_overall %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  pull(ensembl_gene_id)

top_egfr_genes_all <- top_egfr_genes_all[
  top_egfr_genes_all %in% rownames(diff_mat_all)
]
heat_mat_all_egfr <- diff_mat_all[top_egfr_genes_all, ]
annotation_col_all_egfr <- patient_pheno_all %>%
  dplyr::select(ID, EGFR, sex, stage, smoking_status) %>%
  as.data.frame()

rownames(annotation_col_all_egfr) <- annotation_col_all_egfr$ID
annotation_col_all_egfr$ID <- NULL

annotation_col_all_egfr <- annotation_col_all_egfr[colnames(heat_mat_all_egfr), ]

all(rownames(annotation_col_all_egfr) == colnames(heat_mat_all_egfr))

pheatmap(
  heat_mat_all_egfr,
  scale = "row",
  annotation_col = annotation_col_all_egfr,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Top 50 EGFR-associated Tumor-Normal differences: All patients"
)
#Q3:在男生中，EGFR三類對Tumor-Normalexpression change的影響----

male_patient_ids <- unique(male_info$ID)
diff_mat_male <- sapply(male_patient_ids, function(pid) {tumor_sample <- male_info$id[male_info$ID == pid & male_info$sample_type == "Tumor"]
  normal_sample <- male_info$id[
  male_info$ID == pid & male_info$sample_type == "Normal"]
  logCPM_after[, tumor_sample] - logCPM_after[, normal_sample]})

colnames(diff_mat_male) <- male_patient_ids
dim(diff_mat_male)

#一個人有兩列 要利用diff變成一列
patient_pheno_male <- male_info %>%distinct(ID, .keep_all = TRUE) %>%filter(ID %in% colnames(diff_mat_male)) %>%arrange(match(ID, colnames(diff_mat_male)))

#設定ref:WT(dummy)
table(patient_pheno_male$EGFR)
design_male_egfr <- model.matrix(~ EGFR,data = patient_pheno_male)
colnames(design_male_egfr)

#limma
fit_male_egfr <- lmFit(diff_mat_male,design_male_egfr) #一次對所有gene都fitting model
fit_male_egfr <- eBayes(fit_male_egfr) #empirical Bayes moderation:為了穩定每個基因的變異
coef_egfr_male <- grep("^EGFR", colnames(design_male_egfr))

coef_egfr_male
colnames(design_male_egfr)[coef_egfr_male]

#overall F test:男生中，EGFR三類之間的Tumor-Normal expression difference是否整體有差？
#multiple testing correction->為了校正假陽性 , BH:所有P先由小到大排序 在各自呈上總數/排名數 在看有沒有小於0.05
deg_male_egfr_overall <- topTable(fit_male_egfr,coef = coef_egfr_male,number = Inf,adjust.method = "BH") %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%left_join(gene_info, by = "ensembl_gene_id")


head(deg_male_egfr_overall)


deg_male_egfr_overall <- deg_male_egfr_overall %>%
  mutate(
    EGFR_status = case_when(
      adj.P.Val < 0.05 ~ "EGFR-associated",
      TRUE ~ "Not significant"
    )
  )

table(deg_male_egfr_overall$EGFR_status)


#抓出排序第一名的 gene然畫圖:SCH4
#男生中，SHC4的Tumor - Normal表達差異在EGFR三類之間是否不同
#判斷:y>0Tumor中比Normal高
#top gene:"ENSG00000185634"->"SHC4"


top_gene_id <- deg_male_egfr_overall$ensembl_gene_id[1]
top_gene_id 
top_gene_name <- deg_male_egfr_overall$external_gene_name[1]
top_gene_name
plot_df_male <- patient_pheno_male %>%
  mutate(
    gene_diff = as.numeric(diff_mat_male[top_gene_id, ])
  )

plot_df_male <- plot_df_male %>%
  mutate(
    EGFR = factor(EGFR, levels = c("WT", "exon19del", "L858R"))
  )

ggplot(plot_df_male, aes(x = EGFR, y = gene_diff, fill = EGFR)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  labs(
    title = paste0("In male:
Tumor-Normal expression difference of ", top_gene_name, " by EGFR category"),
    x = "EGFR category",
    y = "Tumor - Normal logCPM difference"
  )


#f test
deg_male_egfr_overall <-deg_male_egfr_overall %>%mutate(EGFR_status = case_when(adj.P.Val < 0.05 ~ "EGFR-associated",TRUE ~ "Not significant"))
ggplot(deg_male_egfr_overall, aes(x = F, y = -log10(adj.P.Val), color = EGFR_status)) +
  geom_point(alpha = 0.6, size = 1.2) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "In male,
EGFR three-category effect on Tumor-Normal expression difference",
    x = "F-statistic",
    y = "-log10 adjusted p-value"
  )


#heatmap
library(pheatmap)
#前50個
top_egfr_genes_male <- deg_male_egfr_overall %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  pull(ensembl_gene_id)

top_egfr_genes_male <- top_egfr_genes_male[
  top_egfr_genes_male %in% rownames(diff_mat_male)
]
heat_mat_male_egfr <- diff_mat_male[top_egfr_genes_male, ]

annotation_col_male_egfr <- patient_pheno_male%>%
  dplyr::select(ID, EGFR, sex, stage, smoking_status) %>%
  as.data.frame()

rownames(annotation_col_male_egfr) <- annotation_col_male_egfr$ID
annotation_col_male_egfr$ID <- NULL

annotation_col_male_egfr <- annotation_col_male_egfr[colnames(heat_mat_male_egfr), ]

all(rownames(annotation_col_male_egfr) == colnames(heat_mat_male_egfr))

pheatmap(
  heat_mat_male_egfr,
  scale = "row",
  annotation_col = annotation_col_male_egfr,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Top 50 EGFR-associated Tumor-Normal differences:in male"
)


#Q4:在女生中，EGFR三類對Tumor-Normalexpression change的影----
#female:66
table(pheno_info$sex)
female_info <- pheno_info %>%filter(sex == "Female")
table(female_info$sample_type)


female_patient_ids <- unique(female_info$ID)

diff_mat_female <- sapply(female_patient_ids, function(pid) {tumor_sample <- female_info$id[female_info$ID == pid & female_info$sample_type == "Tumor"]
  normal_sample <- female_info$id[female_info$ID == pid & female_info$sample_type == "Normal"]
  logCPM_after[, tumor_sample] - logCPM_after[, normal_sample]
})

colnames(diff_mat_female) <- female_patient_ids
dim(diff_mat_female)
diff_mat_female[1:5, 1:5]

patient_pheno_female <- female_info %>%distinct(ID, .keep_all = TRUE) %>%filter(ID %in% colnames(diff_mat_female)) %>%arrange(match(ID, colnames(diff_mat_female)))
all(patient_pheno_female$ID == colnames(diff_mat_female))

table(patient_pheno_female$EGFR)

patient_pheno_female$EGFR <- factor(patient_pheno_female$EGFR)
levels(patient_pheno_female$EGFR)
#wt為ref, 比較:
#exon19del vs WT
#L858R vs WT
patient_pheno_female$EGFR <- relevel(patient_pheno_female$EGFR,ref = "WT")
levels(patient_pheno_female$EGFR)
design_female_egfr <- model.matrix( ~ EGFR,data = patient_pheno_female)
colnames(design_female_egfr)

#limma f test
fit_female_egfr <- lmFit(diff_mat_female,design_female_egfr)
fit_female_egfr <- eBayes(fit_female_egfr)
coef_egfr_female <- grep("^EGFR", colnames(design_female_egfr))
coef_egfr_female
colnames(design_female_egfr)[coef_egfr_female]
deg_female_egfr_overall <- topTable(fit_female_egfr,coef = coef_egfr_female,number = Inf,adjust.method = "BH") %>%tibble::rownames_to_column("ensembl_gene_id") %>%left_join(gene_info, by = "ensembl_gene_id")


deg_female_egfr_overall <- deg_female_egfr_overall %>%mutate(EGFR_status = case_when(adj.P.Val < 0.05 ~ "EGFR-associated",TRUE ~ "Not significant"))

table(deg_female_egfr_overall$EGFR_status)
min(deg_female_egfr_overall$adj.P.Val, na.rm = TRUE)


deg_female_egfr_overall %>%filter(EGFR_status == "EGFR-associated") %>%dplyr::select(ensembl_gene_id,EGFRexon19del,EGFRL858R,
    AveExpr, F,P.Value,adj.P.Val,external_gene_name)



#f test plot
ggplot(deg_female_egfr_overall,aes(x = F, y = -log10(adj.P.Val), color = EGFR_status)) +geom_point(alpha = 0.6, size = 1.2) +theme_bw() +geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "In female:
EGFR three-category effect on Tumor-Normal expression difference",
    x = "F-statistic",
    y = "-log10 adjusted p-value"
  )

#boxplot
#"ENSG00000086696"->"HSD17B2"
top_gene_id_female <- deg_female_egfr_overall$ensembl_gene_id[1]
top_gene_name_female <- deg_female_egfr_overall$external_gene_name[1]
top_gene_id_female

plot_df_female <- patient_pheno_female %>%mutate(gene_diff = as.numeric(diff_mat_female[top_gene_id_female, ]))

ggplot(plot_df_female, aes(x = EGFR, y = gene_diff, fill = EGFR)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  labs(
    title = paste0("In female:
Tumor-Normal expression difference of ", top_gene_name_female, " by EGFR category"),
    x = "EGFR category",
    y = "Tumor - Normal logCPM difference"
  )


#heatmap
library(pheatmap)
top_egfr_genes_female <- deg_female_egfr_overall %>%arrange(adj.P.Val) %>%slice_head(n = 50) %>%pull(ensembl_gene_id)
top_egfr_genes_female <- top_egfr_genes_female[top_egfr_genes_female %in% rownames(diff_mat_female)]
heat_mat_female_egfr <- diff_mat_female[top_egfr_genes_female, ]
#不用smoking
annotation_col_female_egfr <- patient_pheno_female %>%dplyr::select(ID, EGFR, stage) %>%as.data.frame()
rownames(annotation_col_female_egfr) <- annotation_col_female_egfr$ID
annotation_col_female_egfr$ID <- NULL
annotation_col_female_egfr <- annotation_col_female_egfr[colnames(heat_mat_female_egfr), ]

all(rownames(annotation_col_female_egfr) == colnames(heat_mat_female_egfr))

pheatmap(
  heat_mat_female_egfr,
  scale = "row",
  annotation_col = annotation_col_female_egfr,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Top 50 EGFR-associated Tumor-Normal differences,in female"
)

#特別看"HSD17B2"
deg_female_egfr_overall %>%filter(external_gene_name == "HSD17B2") %>%
  dplyr::select(
    ensembl_gene_id,
    external_gene_name,
    EGFRexon19del,
    EGFRL858R,
    AveExpr,
    F,
    P.Value,
    adj.P.Val
  )
deg_female_egfr_overall
target_gene <- "HSD17B2"

target_gene_id <- gene_info %>%
  filter(external_gene_name == target_gene) %>%
  pull(ensembl_gene_id) %>%
  unique()

target_gene_id <- target_gene_id[1]

plot_df_female_hsd <- patient_pheno_female %>%
  mutate(
    gene_diff = as.numeric(diff_mat_female[target_gene_id, ])
  )

ggplot(plot_df_female_hsd, aes(x = EGFR, y = gene_diff, fill = EGFR)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  labs(
    title = "In female: Tumor-Normal expression difference of HSD17B2 by EGFR category",
    x = "EGFR category",
    y = "Tumor - Normal logCPM difference"
  )


#pairwise->看是哪一組有差

# exon19del vs WT
deg_female_exon19_vs_WT <- topTable(
  fit_female_egfr,
  coef = "EGFRexon19del",
  number = Inf,
  adjust.method = "BH"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info, by = "ensembl_gene_id")

deg_female_exon19_vs_WT %>%
  filter(external_gene_name == "HSD17B2")


# L858R vs WT
deg_female_L858R_vs_WT <- topTable(
  fit_female_egfr,
  coef = "EGFRL858R",
  number = Inf,
  adjust.method = "BH"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info, by = "ensembl_gene_id")

deg_female_L858R_vs_WT %>%
  filter(external_gene_name == "HSD17B2")


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
survdiff(
  Surv(rfs_time, recurrence_event) ~ smoking_egfr_group,
  data =patient_pheno_male
)



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


