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
diff_mat_all[1:5, 1:5] #logCPM difference
#logCPM>1:表示差兩倍
# 計算每個 EGFR group 中，每個 gene 的平均 Tumor-Normal difference
egfr_levels <- levels(patient_pheno_all$EGFR)
mean_diff_by_egfr <- sapply(egfr_levels, function(g) {
  rowMeans(diff_mat_all[, patient_pheno_all$EGFR == g, drop = FALSE], na.rm = TRUE)
})

colnames(mean_diff_by_egfr) <- paste0("mean_diff_", make.names(egfr_levels))

mean_diff_df <- as.data.frame(mean_diff_by_egfr) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  mutate(
    max_abs_mean_diff = apply(
      dplyr::select(., starts_with("mean_diff_")),
      1,
      function(x) max(abs(x), na.rm = TRUE)
    )
  )



#用diff轉成一個病人一列
patient_pheno_all <- pheno_info %>%distinct(ID, .keep_all = TRUE) %>%filter(ID %in% colnames(diff_mat_all)) %>%arrange(match(ID, colnames(diff_mat_all)))
dim(patient_pheno_all)
all(patient_pheno_all$ID == colnames(diff_mat_all)) #檢查順序是否一致

table(patient_pheno_all$EGFR)
patient_pheno_all$EGFR <- factor(patient_pheno_all$EGFR)
levels(patient_pheno_all$EGFR)

#設定ref:WT(dummy)
patient_pheno_all$EGFR <- relevel(patient_pheno_all$EGFR, ref = "WT")
design_all_egfr <- model.matrix(~ EGFR, data = patient_pheno_all)
colnames(design_all_egfr)

#一次看很多基因用limma overall f test
#limma 則會把所有 genes 的 variance pattern 一起考慮，讓統計結果比較穩
fit_all_egfr <- lmFit(diff_mat_all, design_all_egfr)
fit_all_egfr <- eBayes(fit_all_egfr)

#EGFR 三類之間的 Tumor-Normal expression difference 是否整體有差異？
#不做fdr(不用BH校正)且P<0.05就顯著 並挑出POST HOC也顯著的才畫圖

coef_egfr_all <- grep("^EGFR", colnames(design_all_egfr))
colnames(design_all_egfr)[coef_egfr_all]

deg_all_egfr_overall <- topTable(fit_all_egfr,coef = coef_egfr_all,
  number = Inf,adjust.method = "none") %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info, by = "ensembl_gene_id")
#合併log2CPM>1:
deg_all_egfr_overall <- deg_all_egfr_overall %>%left_join(mean_diff_df, by = "ensembl_gene_id")
head(deg_all_egfr_overall )

#定義GENE狀態 如果有顯著就是associated
deg_all_egfr_overall <- deg_all_egfr_overall %>%
  mutate(EGFR_status = case_when(P.Value < 0.05 ~ "EGFR-associated",TRUE ~ "Not significant"))
table(deg_all_egfr_overall$EGFR_status)


#post hoc 看三組egfr不同是哪組不同

# 建立沒有 intercept 的 design，方便做 pairwise comparison
design_posthoc <- model.matrix(~ 0 + EGFR, data = patient_pheno_all)
# 讓欄位名稱變乾淨
colnames(design_posthoc) <- make.names(gsub("EGFR", "", colnames(design_posthoc)))
colnames(design_posthoc)
#fit model
fit_posthoc <- lmFit(diff_mat_all, design_posthoc)
contrast_posthoc <- makeContrasts(exon19_vs_WT = exon19del - WT,
  L858R_vs_WT = L858R - WT,
  exon19_vs_point = exon19del - L858R,
  levels = design_posthoc
)
colnames(design_posthoc)
fit_posthoc2 <- contrasts.fit(fit_posthoc, contrast_posthoc)
fit_posthoc2 <- eBayes(fit_posthoc2)

posthoc_exon19_WT <- topTable(
  fit_posthoc2,
  coef = "exon19_vs_WT",
  number = Inf,
  adjust.method = "none"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(
    ensembl_gene_id,
    logFC_exon19_vs_WT = logFC,
    P_exon19_vs_WT = P.Value
  )

posthoc_L858R_WT <- topTable(
  fit_posthoc2,
  coef = "L858R_vs_WT",
  number = Inf,
  adjust.method = "none"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(
    ensembl_gene_id,
    logFC_L858R_vs_WT = logFC,
    P_L858R_vs_WT = P.Value
  )

posthoc_exon19_L858R <- topTable(
  fit_posthoc2,
  coef = "exon19_vs_point",
  number = Inf,
  adjust.method = "none"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(
    ensembl_gene_id,
    logFC_exon19_vs_L858R = logFC,
    P_exon19_vs_L858R = P.Value
  )

deg_all_egfr_final <- deg_all_egfr_overall %>%
  left_join(posthoc_exon19_WT, by = "ensembl_gene_id") %>%
  left_join(posthoc_L858R_WT, by = "ensembl_gene_id") %>%
  left_join(posthoc_exon19_L858R, by = "ensembl_gene_id") %>%
  mutate(
    Ftest_sig = P.Value < 0.05,
    two_fold_sig = max_abs_mean_diff > 1,
    posthoc_sig =
      P_exon19_vs_WT < 0.05 |
      P_L858R_vs_WT < 0.05 |
      P_exon19_vs_L858R < 0.05,
    final_status = case_when(
      Ftest_sig & two_fold_sig & posthoc_sig ~ "Selected candidate genes",
      TRUE ~ "Not selected"
    )
  )


##F-test P.Value < 0.05
#log2CPM > 1
#至少一個 post hoc p-value < 0.05
table(deg_all_egfr_final$final_status)


#全部人：挑選出來的gene boxplot
selected_egfr_genes <- deg_all_egfr_final %>%
  filter(final_status == "Selected candidate genes") %>%
  arrange(P.Value)

selected_egfr_genes %>%
  dplyr::select(
    ensembl_gene_id,
    external_gene_name,
    P.Value,
    F,
    max_abs_mean_diff,
    starts_with("mean_diff_"),
    starts_with("P_")
  ) %>%
  head(20)


top_gene_id_all <- selected_egfr_genes$ensembl_gene_id[1]
top_gene_name_all <- selected_egfr_genes$external_gene_name[1]

top_gene_id_all
top_gene_name_all

plot_df_all<- patient_pheno_all %>%
  mutate(gene_diff = as.numeric(diff_mat_all[top_gene_id_all, ])
  )

ggplot(plot_df_all, aes(x = EGFR, y = gene_diff, fill = EGFR)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  labs(title = paste0(
      "Tumor-Normal expression difference by EGFR category",
        top_gene_name_all,
        " by EGFR category"
      ),
    x = "EGFR category",
    y = "Tumor - Normal log2CPM difference"
  )

#f test
deg_all_egfr_overall <- deg_all_egfr_overall %>%
  mutate(
    EGFR_status = case_when(
      P.Value < 0.05 ~ "EGFR-associated",
      TRUE ~ "Not significant"
    )
  )
ggplot(deg_all_egfr_overall, aes(x = F, y = -log10(P.Value), color = EGFR_status)) +
  geom_point(alpha = 0.6, size = 1.2) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "EGFR effect on Tumor-Normal expression difference",
    x = "F-statistic",
    y = "-log10 p-value"
  )
#共1314個皆有符合條件
sum(deg_all_egfr_final$Ftest_sig, na.rm = TRUE)
sum(deg_all_egfr_final$two_fold_sig, na.rm = TRUE)
sum(deg_all_egfr_final$posthoc_sig, na.rm = TRUE)
sum(deg_all_egfr_final$Ftest_sig & 
      deg_all_egfr_final$two_fold_sig & 
      deg_all_egfr_final$posthoc_sig, 
    na.rm = TRUE)
criteria_count <- data.frame(
  Criteria = c(
    "F-test nominal p < 0.05",
    "At least one EGFR group |mean Tumor-Normal log2CPM difference| > 1",
    "At least one post hoc pairwise p < 0.05",
    "All three criteria"
  ),
  N_genes = c(
    sum(deg_all_egfr_final$Ftest_sig, na.rm = TRUE),
    sum(deg_all_egfr_final$two_fold_sig, na.rm = TRUE),
    sum(deg_all_egfr_final$posthoc_sig, na.rm = TRUE),
    sum(deg_all_egfr_final$Ftest_sig &
          deg_all_egfr_final$two_fold_sig &
          deg_all_egfr_final$posthoc_sig,
        na.rm = TRUE)
  )
)
criteria_count
selected_egfr_genes <- deg_all_egfr_final %>%filter(Ftest_sig == TRUE, two_fold_sig == TRUE, posthoc_sig == TRUE) %>%arrange(P.Value)
nrow(selected_egfr_genes)


#heatmap->只放基因顯著且|log2CPM|>1且至少一個POST HOC顯著
library(pheatmap)

top_egfr_genes_all <- selected_egfr_genes %>%
  arrange(P.Value) %>%
  pull(ensembl_gene_id)

nrow(selected_egfr_genes) #=1314

#呈現前50個
if(length(top_egfr_genes_all) > 50) {
  top_egfr_genes_all <- selected_egfr_genes %>%
    arrange(P.Value) %>%
    slice_head(n = 50) %>%
    pull(ensembl_gene_id)
}



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
  main = "EGFR-associated Tumor-Normal expression differences"
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
#不做fdr(不用BH校正)且P<0.05就顯著 並挑出POST HOC也顯著的才畫圖
deg_male_egfr_overall <- topTable(fit_male_egfr,coef = coef_egfr_male,number = Inf,adjust.method = "none") %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%left_join(gene_info, by = "ensembl_gene_id")


head(deg_male_egfr_overall)

deg_male_egfr_overall <- deg_male_egfr_overall %>%
  mutate(
    EGFR_status = case_when(
      P.Value < 0.05 ~ "EGFR-associated",
      TRUE ~ "Not significant"
    )
  )

table(deg_male_egfr_overall$EGFR_status)

egfr_levels_male <- levels(patient_pheno_male$EGFR)

mean_diff_by_egfr_male <- sapply(egfr_levels_male, function(g) {
  rowMeans(
    diff_mat_male[, patient_pheno_male$EGFR == g, drop = FALSE],
    na.rm = TRUE
  )
})
colnames(mean_diff_by_egfr_male) <- paste0(
  "mean_diff_",
  make.names(egfr_levels_male)
)

mean_diff_male_df <- as.data.frame(mean_diff_by_egfr_male) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  mutate(
    max_abs_mean_diff = apply(
      dplyr::select(., starts_with("mean_diff_")),
      1,
      function(x) max(abs(x), na.rm = TRUE)
    )
  )

deg_male_egfr_overall <- deg_male_egfr_overall %>%
  left_join(mean_diff_male_df, by = "ensembl_gene_id")

#post hoc
design_posthoc_male <- model.matrix(~ 0 + EGFR, data = patient_pheno_male)

colnames(design_posthoc_male) <- make.names(
  gsub("EGFR", "", colnames(design_posthoc_male))
)

colnames(design_posthoc_male)

fit_posthoc_male <- lmFit(diff_mat_male, design_posthoc_male)

contrast_posthoc_male <- makeContrasts(
  exon19_vs_WT = exon19del - WT,
  L858R_vs_WT = L858R - WT,
  exon19_vs_L858R = exon19del - L858R,
  levels = design_posthoc_male
)

posthoc_male_exon19_WT <- topTable(
  fit_posthoc_male2,
  coef = "exon19_vs_WT",
  number = Inf,
  adjust.method = "none"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(
    ensembl_gene_id,
    logFC_exon19_vs_WT = logFC,
    P_exon19_vs_WT = P.Value
  )

posthoc_male_L858R_WT <- topTable(
  fit_posthoc_male2,
  coef = "L858R_vs_WT",
  number = Inf,
  adjust.method = "none"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(
    ensembl_gene_id,
    logFC_L858R_vs_WT = logFC,
    P_L858R_vs_WT = P.Value
  )

posthoc_male_exon19_L858R <- topTable(
  fit_posthoc_male2,
  coef = "exon19_vs_L858R",
  number = Inf,
  adjust.method = "none"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(
    ensembl_gene_id,
    logFC_exon19_vs_L858R = logFC,
    P_exon19_vs_L858R = P.Value
  )

#滿足三個條件: 432
deg_male_egfr_final <- deg_male_egfr_overall %>%
  left_join(posthoc_male_exon19_WT, by = "ensembl_gene_id") %>%
  left_join(posthoc_male_L858R_WT, by = "ensembl_gene_id") %>%
  left_join(posthoc_male_exon19_L858R, by = "ensembl_gene_id") %>%
  mutate(
    Ftest_sig = P.Value < 0.05,
    two_fold_sig = max_abs_mean_diff > 1,
    posthoc_sig =
      P_exon19_vs_WT < 0.05 |
      P_L858R_vs_WT < 0.05 |
      P_exon19_vs_L858R < 0.05,
    final_status = case_when(
      Ftest_sig & two_fold_sig & posthoc_sig ~ "Selected candidate genes",
      TRUE ~ "Not selected"
    )
  )

table(deg_male_egfr_final$final_status)


#挑挑出符合三格條件的基因
selected_egfr_genes_male <- deg_male_egfr_final %>%
  filter(final_status == "Selected candidate genes") %>%
  arrange(P.Value)

nrow(selected_egfr_genes_male)

selected_egfr_genes_male %>%
  dplyr::select(
    ensembl_gene_id,
    external_gene_name,
    P.Value,
    F,
    max_abs_mean_diff,
    starts_with("mean_diff_"),
    starts_with("P_")
  ) %>%
  head(20)

#boxplot
top_gene_id_male <- selected_egfr_genes_male$ensembl_gene_id[1]
top_gene_name_male <- selected_egfr_genes_male$external_gene_name[1]

top_gene_id_male
top_gene_name_male

plot_df_male <- patient_pheno_male %>%
  mutate(
    gene_diff = as.numeric(diff_mat_male[top_gene_id_male, ])
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
    title = paste0(
      "In male: Tumor-Normal expression difference of ",
      top_gene_name_male,
      " by EGFR category"
    ),
    x = "EGFR category",
    y = "Tumor - Normal log2CPM difference"
  )


#f test
deg_all_egfr_overall <- deg_all_egfr_overall %>%
  mutate(
    EGFR_status = case_when(
      P.Value < 0.05 ~ "EGFR-associated",
      TRUE ~ "Not significant"
    )
  )

ggplot(deg_all_egfr_overall, aes(x = F, y = -log10(P.Value), color = EGFR_status)) +
  geom_point(alpha = 0.6, size = 1.2) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(
    values = c(
      "EGFR-associated" = "#F8766D",
      "Not significant" = "#00BFC4"
    )
  ) +
  labs(
    title = "In male, EGFR effect on Tumor-Normal expression difference",
    x = "F-statistic",
    y = "-log10 p-value",
    color = "EGFR_status"
  )

#heatmap
top_egfr_genes_male <- selected_egfr_genes_male %>%
  arrange(P.Value) %>%
  pull(ensembl_gene_id)

# 如果超過50個，為了可讀性呈現前50個
if(length(top_egfr_genes_male) > 50) {
  top_egfr_genes_male <- selected_egfr_genes_male %>%
    arrange(P.Value) %>%
    slice_head(n = 50) %>%
    pull(ensembl_gene_id)
}

top_egfr_genes_male <- top_egfr_genes_male[
  top_egfr_genes_male %in% rownames(diff_mat_male)
]

heat_mat_male_egfr <- diff_mat_male[top_egfr_genes_male, , drop = FALSE]

annotation_col_male_egfr <- patient_pheno_male %>%
  dplyr::select(ID, EGFR, stage, smoking_status) %>%
  mutate(
    EGFR = ifelse(is.na(EGFR) | EGFR == "", "Unknown", as.character(EGFR)),
    stage = ifelse(is.na(stage) | stage == "", "Unknown", as.character(stage)),
    smoking_status = ifelse(is.na(smoking_status) | smoking_status == "", "Unknown", as.character(smoking_status))
  ) %>%
  as.data.frame()

rownames(annotation_col_male_egfr) <- annotation_col_male_egfr$ID
annotation_col_male_egfr$ID <- NULL

annotation_col_male_egfr <- annotation_col_male_egfr[
  colnames(heat_mat_male_egfr),
  ,
  drop = FALSE
]

all(rownames(annotation_col_male_egfr) == colnames(heat_mat_male_egfr))

pheatmap(
  heat_mat_male_egfr,
  scale = "row",
  annotation_col = annotation_col_male_egfr,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "EGFR-associated Tumor-Normal expression differences: Male"
)












#Q4:在女生中，EGFR三類對Tumor-Normalexpression change的影----
#F-test nominal P.Value < 0.05
#至少一組 EGFR subtype 的 |mean Tumor-Normal log2CPM difference| > 1
#至少一組 post hoc pairwise comparison P.Value < 0.05
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
deg_female_egfr_overall <- topTable(
  fit_female_egfr,
  coef = coef_egfr_female,
  number = Inf,
  adjust.method = "none"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info, by = "ensembl_gene_id")

head(deg_female_egfr_overall)
egfr_levels_female <- levels(patient_pheno_female$EGFR)

mean_diff_by_egfr_female <- sapply(egfr_levels_female, function(g) {
  rowMeans(
    diff_mat_female[, patient_pheno_female$EGFR == g, drop = FALSE],
    na.rm = TRUE
  )
})

colnames(mean_diff_by_egfr_female) <- paste0(
  "mean_diff_",
  make.names(egfr_levels_female)
)

mean_diff_female_df <- as.data.frame(mean_diff_by_egfr_female) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  mutate(
    max_abs_mean_diff = apply(
      dplyr::select(., starts_with("mean_diff_")),
      1,
      function(x) max(abs(x), na.rm = TRUE)
    )
  )

deg_female_egfr_overall <- deg_female_egfr_overall %>%
  left_join(mean_diff_female_df, by = "ensembl_gene_id") %>%
  mutate(
    EGFR_status = case_when(
      P.Value < 0.05 ~ "EGFR-associated",
      TRUE ~ "Not significant"
    )
  )

table(deg_female_egfr_overall$EGFR_status)

#post hoc
design_posthoc_female <- model.matrix(~ 0 + EGFR, data = patient_pheno_female)

colnames(design_posthoc_female) <- make.names(
  gsub("EGFR", "", colnames(design_posthoc_female))
)

colnames(design_posthoc_female)

fit_posthoc_female <- lmFit(diff_mat_female, design_posthoc_female)

contrast_posthoc_female <- makeContrasts(
  exon19_vs_WT = exon19del - WT,
  L858R_vs_WT = L858R - WT,
  exon19_vs_L858R = exon19del - L858R,
  levels = design_posthoc_female
)

fit_posthoc_female2 <- contrasts.fit(fit_posthoc_female, contrast_posthoc_female)
fit_posthoc_female2 <- eBayes(fit_posthoc_female2)

colnames(fit_posthoc_female2$coefficients)
posthoc_female_exon19_WT <- topTable(
  fit_posthoc_female2,
  coef = "exon19_vs_WT",
  number = Inf,
  adjust.method = "none"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(
    ensembl_gene_id,
    logFC_exon19_vs_WT = logFC,
    P_exon19_vs_WT = P.Value
  )

posthoc_female_L858R_WT <- topTable(
  fit_posthoc_female2,
  coef = "L858R_vs_WT",
  number = Inf,
  adjust.method = "none"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(
    ensembl_gene_id,
    logFC_L858R_vs_WT = logFC,
    P_L858R_vs_WT = P.Value
  )

posthoc_female_exon19_L858R <- topTable(
  fit_posthoc_female2,
  coef = "exon19_vs_L858R",
  number = Inf,
  adjust.method = "none"
) %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(
    ensembl_gene_id,
    logFC_exon19_vs_L858R = logFC,
    P_exon19_vs_L858R = P.Value
  )


#三個條件皆符合:1154
deg_female_egfr_final <- deg_female_egfr_overall %>%
  left_join(posthoc_female_exon19_WT, by = "ensembl_gene_id") %>%
  left_join(posthoc_female_L858R_WT, by = "ensembl_gene_id") %>%
  left_join(posthoc_female_exon19_L858R, by = "ensembl_gene_id") %>%
  mutate(
    Ftest_sig = P.Value < 0.05,
    two_fold_sig = max_abs_mean_diff > 1,
    posthoc_sig =
      P_exon19_vs_WT < 0.05 |
      P_L858R_vs_WT < 0.05 |
      P_exon19_vs_L858R < 0.05,
    final_status = case_when(
      Ftest_sig & two_fold_sig & posthoc_sig ~ "Selected candidate genes",
      TRUE ~ "Not selected"
    )
  )

table(deg_female_egfr_final$final_status)
selected_egfr_genes_female <- deg_female_egfr_final %>%
  filter(final_status == "Selected candidate genes") %>%
  arrange(P.Value)

nrow(selected_egfr_genes_female)

selected_egfr_genes_female %>%
  dplyr::select(
    ensembl_gene_id,
    external_gene_name,
    P.Value,
    F,
    max_abs_mean_diff,
    starts_with("mean_diff_"),
    starts_with("P_")
  ) %>%
  head(20)

#f test
deg_female_egfr_final <- deg_female_egfr_final %>%
  mutate(
    EGFR_status = case_when(
      P.Value < 0.05 ~ "EGFR-associated",
      TRUE ~ "Not significant"
    )
  )

ggplot(
  deg_female_egfr_final,
  aes(x = F, y = -log10(P.Value), color = EGFR_status)
) +
  geom_point(alpha = 0.6, size = 1.2) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(
    values = c(
      "EGFR-associated" = "#F8766D",
      "Not significant" = "#00BFC4"
    )
  ) +
  labs(
    title = "In female: EGFR effect on Tumor-Normal expression difference",
    x = "F-statistic",
    y = "-log10 p-value",
    color = "EGFR_status"
  )

#boxplot
top_gene_id_female <- selected_egfr_genes_female$ensembl_gene_id[1]
top_gene_name_female <- selected_egfr_genes_female$external_gene_name[1]

top_gene_id_female
top_gene_name_female

plot_df_female <- patient_pheno_female %>%
  mutate(
    gene_diff = as.numeric(diff_mat_female[top_gene_id_female, ])
  ) %>%
  mutate(
    EGFR = factor(EGFR, levels = c("WT", "exon19del", "L858R"))
  )

ggplot(plot_df_female, aes(x = EGFR, y = gene_diff, fill = EGFR)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  labs(
    title = paste0(
      "In female: Tumor-Normal expression difference of ",
      top_gene_name_female,
      " by EGFR category"
    ),
    x = "EGFR category",
    y = "Tumor - Normal log2CPM difference"
  )
#heatmap
top_egfr_genes_female <- selected_egfr_genes_female %>%
  arrange(P.Value) %>%
  pull(ensembl_gene_id)

nrow(selected_egfr_genes_female)

# 如果超過50個，為了可讀性呈現前50個
if(length(top_egfr_genes_female) > 50) {
  top_egfr_genes_female <- selected_egfr_genes_female %>%
    arrange(P.Value) %>%
    slice_head(n = 50) %>%
    pull(ensembl_gene_id)
}

top_egfr_genes_female <- top_egfr_genes_female[
  top_egfr_genes_female %in% rownames(diff_mat_female)
]

heat_mat_female_egfr <- diff_mat_female[top_egfr_genes_female, , drop = FALSE]

annotation_col_female_egfr <- patient_pheno_female %>%
  dplyr::select(ID, EGFR, stage) %>%
  mutate(
    EGFR = ifelse(is.na(EGFR) | EGFR == "", "Unknown", as.character(EGFR)),
    stage = ifelse(is.na(stage) | stage == "", "Unknown", as.character(stage))
  ) %>%
  as.data.frame()

rownames(annotation_col_female_egfr) <- annotation_col_female_egfr$ID
annotation_col_female_egfr$ID <- NULL

annotation_col_female_egfr <- annotation_col_female_egfr[
  colnames(heat_mat_female_egfr),
  ,
  drop = FALSE
]

all(rownames(annotation_col_female_egfr) == colnames(heat_mat_female_egfr))

pheatmap(
  heat_mat_female_egfr,
  scale = "row",
  annotation_col = annotation_col_female_egfr,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "EGFR-associated Tumor-Normal expression differences: Female"
)


