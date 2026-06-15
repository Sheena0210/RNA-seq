# Bioconductor 套件----
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
BiocManager::install(c(
  "edgeR",
  "limma",
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db"
))

library(tidyverse)
library(readxl)
library(edgeR)
library(limma)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(survival)
library(survminer)
library(clusterProfiler)
library(org.Hs.eg.db)

library(readxl)

expression<- read_xlsx( "C:/Users/Sh33na/Desktop/summer_intern/geneExpression_60samples_readcount.xlsx")
pheno<- read_xlsx("C:/Users/Sh33na/Desktop/summer_intern/pheno.xlsx")

dim(expression)
head(expression)
colnames(expression)

dim(pheno)
head(pheno)
colnames(pheno)
#ID   age   sex   tumor_tissue_type   tumor_laterality_site
#tumor_size  stage  EGFR  smoking_status  vital_status  days_to_last_contact_or_death   tumor_recurrence   days_to_tumor_recurrence
 
#Pre-processing----

#將兩個檔案的欄位一致

#data:expression ->  count_exp:每一列是 gene，每一欄是 sample  --
gene_info <- expression%>%dplyr::select(ensembl_gene_id, external_gene_name) #挑出ensembl_gene_id、external_gene_name
count_exp <- expression%>%dplyr::select(-ensembl_gene_id,-external_gene_name) #拿掉ensembl_gene_id、external_gene_name
#放對應的gene name
count_exp <- expression%>%dplyr::select(-ensembl_gene_id, -external_gene_name) %>%as.data.frame()
rownames(count_exp) <- expression$ensembl_gene_id

#data frame->matrix(每一列是 gene，每一欄是 sample)
count_exp<- as.matrix(count_exp)
mode(count_exp) <- "numeric"     #轉成數值
dim(count_exp)
count_exp[1:5, 1:5]

#為了跟pheno matching 要把ID間隔出來->  info:每一列是 sample，記錄--
info<-data.frame(id=colnames(count_exp),stringsAsFactors = FALSE)
info$ID <- sub("-(T|N)$","", info$id) #只要id
info$sample_type <- ifelse(grepl("-T$", info$id),"Tumor","Normal")
info$sample_type <- factor(info$sample_type,levels = c("Normal", "Tumor"))
head(info)


#data:pheno--
#用ID跟expression matching
colnames(count_exp)
pheno$ID

#合併info&pheno:  pheno_info (N=120)--
pheno_info<-info%>%left_join(pheno,by="ID")
head(pheno_info)

#ID單純只有編號 id有T/N
all(pheno_info$id == colnames(count_exp))

#check 每個人都有一個t 一個n
pair_check <- pheno_info %>%count(ID)
table(pair_check$n)



#過濾低表達基因
library(edgeR)
library(limma)
library(dplyr)
library(ggplot2)
library(tidyr)
dge <- DGEList(counts = count_exp,group = pheno_info$sample_type)

dge$samples #每個樣本的基本資訊
dim(dge)

#normalization(TMM)-----
dge <- calcNormFactors(dge, method = "TMM")
dge$samples

#paired design matrix(因為每個人都有T/N->控制每個病人本身差異後，比較Tumor和Normal的gene expression差異
pheno_info$ID <- factor(pheno_info$ID)
pheno_info$sample_type <- factor(pheno_info$sample_type,levels = c("Normal", "Tumor"))
design <- model.matrix(~ ID + sample_type, data =pheno_info)
colnames(design)

#voom
voom <- voom(dge, design, plot = TRUE)

#limma fitting
fit <- lmFit(voom, design)
fit <- eBayes(fit)
deg_limma <- topTable(fit,coef = "sample_typeTumor",number = Inf,adjust.method = "BH")
deg_limma <- deg_limma %>%tibble::rownames_to_column("ensembl_gene_id") %>%left_join(gene_info, by = "ensembl_gene_id")

deg_limma <- deg_limma %>%
  mutate(
    DEG_status = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Up in Tumor",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down in Tumor",
      TRUE ~ "Not significant"
    )
  )

table(deg_limma$DEG_status)

#DEG visualization：Volcano plot
ggplot(deg_limma, aes(x = logFC, y = -log10(adj.P.Val), color = DEG_status)) +
  geom_point(alpha = 0.6, size = 1.2) +
  theme_bw() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano plot: Tumor vs Normal",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  )


#過濾低表達基因後 + TMM normalization 後的 log2 CPM expression matrix
logCPM_after <- cpm(
  dge,
  log = TRUE,
  prior.count = 1,
  normalized.lib.sizes = TRUE
)
dim(logCPM_after)
logCPM_after[1:5, 1:5]

#DEG visualization：Top 50 heatmap
library(pheatmap)
top_genes <- deg_limma %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  pull(ensembl_gene_id)
heat_mat <- logCPM_after[top_genes, ]
annotation_col <- pheno_info %>%
  dplyr::select(id, sample_type, ID, stage, EGFR, smoking_status) %>%
  as.data.frame()
rownames(annotation_col) <- annotation_col$id
annotation_col$id <- NULL
pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Top 50 DEGs: Tumor vs Normal"
)





