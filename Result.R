library(dplyr)
library(tidyr)
library(cluster)
library(factoextra)
library(umap)
library(matrixStats)
library(survival)
library(ggplot2)
library(mixOmics)
library(enrichplot)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(limma)
library(org.Hs.eg.db)
library(msigdbr)
library(tibble)

# Cox Regression Analysis #

rna_mirna_protein_surv_obj <- Surv(time = merge_rna_mirna_protein_filter_BRCA_survival$OS.time, event = merge_rna_mirna_protein_filter_BRCA_survival$OS)

rna_mirna_protein_cox_model <- coxph(rna_mirna_protein_surv_obj ~
                     Factor1 + Factor2 + Factor3 + Factor4 + Factor5 +
                     Factor6 + Factor7 + Factor8 + Factor9 + Factor10 +
                     Factor11 + Factor12 + Factor13 + Factor14 + Factor15,
                   data = merge_rna_mirna_protein_filter_BRCA_survival)

summary(rna_mirna_protein_cox_model)
cox.zph(rna_mirna_protein_cox_model)

filter_merge_rna_mirna_protein_filter_BRCA_survival <- merge_rna_mirna_protein_filter_BRCA_survival %>%
  select(patient, Factor3, Factor6, Factor13)

rownames(filter_merge_rna_mirna_protein_filter_BRCA_survival) <- filter_merge_rna_mirna_protein_filter_BRCA_survival$patient
filter_merge_rna_mirna_protein_filter_BRCA_survival$patient <- NULL

# Cluster Analysis #

## Optimal k Selection
### 1. Gap Statistic
set.seed(123)
gap_stat <- clusGap(filter_merge_rna_mirna_protein_filter_BRCA_survival, FUN = kmeans, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)
### 2. Silhouette Method
fviz_nbclust(filter_merge_rna_mirna_protein_filter_BRCA_survival, kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")

## k-means clustering
rna_mirna_protein_km_res <- kmeans(filter_merge_rna_mirna_protein_filter_BRCA_survival, centers = 3, nstart = 25)

## UMAP: dimension reduction
### 1) MOG
rna_mirna_protein_umap_res <- umap(filter_merge_rna_mirna_protein_filter_BRCA_survival)
rna_mirna_protein_umap_res <- as.data.frame(rna_mirna_protein_umap_res$layout)
colnames(rna_mirna_protein_umap_res) <- c('UMAP1', 'UMAP2')

rna_mirna_protein_umap_res$cluster <- factor(
  rna_mirna_protein_umap_res$cluster,
  levels = c(1, 2, 3),
  labels = c("MOG1", "MOG2", "MOG3")
)

### 2) PAM50 Subtype
merge_rna_mirna_protein_umap_res <- rna_mirna_protein_umap_res
merge_rna_mirna_protein_umap_res$sample <- rownames(merge_rna_mirna_protein_umap_res)
rownames(merge_rna_mirna_protein_umap_res) <- NULL
merge_rna_mirna_protein_umap_res <- merge(merge_rna_mirna_protein_umap_res, rna_mirna_protein_filter_BRCA_TCGA_PAM50_df, by='sample')

## Proportion of PAM50 subtype in the cluster
merge_rna_mirna_protein_umap_df_count <- merge_rna_mirna_protein_umap_df_survival %>% dplyr::count(cluster, BRCA_Subtype_PAM50)
merge_rna_mirna_protein_umap_df_prop <- merge_rna_mirna_protein_umap_df_count %>% 
  group_by(cluster) %>%
  mutate(perc = n / sum(n) * 100)

## Visualization
### 1) MOG
ggplot(rna_mirna_protein_umap_res, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    color = "Cluster"
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "blue")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

### 2) PAM50 Subtype
ggplot(merge_rna_mirna_protein_umap_res, aes(x = UMAP1, y = UMAP2, color = BRCA_Subtype_PAM50)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    color = "PAM50"
  ) +
  scale_color_manual(values = c("tomato", "yellow4", "turquoise3", "purple")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

### 3) Bar plot of proportion
ggplot(merge_rna_mirna_protein_umap_df_prop,
       aes(x = cluster, y = perc, fill = BRCA_Subtype_PAM50)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = n), 
            position = position_stack(vjust = 0.5),
            color = "black", size = 3.5) +
  theme_minimal() +
  labs(x = "", y = "Percentage (%)", fill = "Subtype") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x  = element_text(size = 12, face='bold'),
        plot.title   = element_text(hjust = 0.5))

# Survival analysis #

## Creat Survival object & Kaplan-Meier fit
### 1) MOG
rna_mirna_protein_surv_obj_cluster <- Surv(time = merge_rna_mirna_protein_umap_df_survival$OS.time,
                                               event = merge_rna_mirna_protein_umap_df_survival$OS)
rna_mirna_protein_fit_cluster <- survfit(rna_mirna_protein_surv_obj_cluster ~ cluster,
                                                     data = merge_rna_mirna_protein_umap_df_survival)

### 2) PAM50 Subtype
rna_mirna_protein_surv_obj_subtype <- Surv(time = merge_rna_mirna_protein_umap_df_survival$OS.time,
                                               event = merge_rna_mirna_protein_umap_df_survival$OS)
rna_mirna_protein_fit_subtype <- survfit(rna_mirna_protein_surv_obj_subtype ~ BRCA_Subtype_PAM50,
                                                     data = merge_rna_mirna_protein_umap_df_survival)

## Log-rank test
### 1) MOG
rna_mirna_protein_surv_diff_subtype <- survdiff(rna_mirna_protein_surv_obj_subtype ~ BRCA_Subtype_PAM50,
                                                            data = merge_rna_mirna_protein_umap_df_survival)
rna_mirna_protein_p_val_subtype <- 1 - pchisq(rna_mirna_protein_surv_diff_subtype$chisq, df = length(rna_mirna_protein_surv_diff_subtype$n) - 1)

### 2) PAM50 Subtype
rna_mirna_protein_surv_diff_cluster <- survdiff(rna_mirna_protein_surv_obj_cluster ~ cluster,
                                                            data = merge_rna_mirna_protein_umap_df_survival)
rna_mirna_protein_p_val_cluster <- 1 - pchisq(rna_mirna_protein_surv_diff_cluster$chisq, df = length(rna_mirna_protein_surv_diff_cluster$n) - 1)

## Pairwise log-rank test
### 1) MOG
rna_mirna_protein_pairwise_results_cluster <- pairwise_survdiff(Surv(OS.time, OS) ~ cluster, data = merge_rna_mirna_protein_umap_df_survival)
rna_mirna_protein_pairwise_results_cluster$p.value

### 2) PAM50 Subtype
rna_mirna_protein_pairwise_results_subtype <- pairwise_survdiff(Surv(OS.time, OS) ~ BRCA_Subtype_PAM50, data = merge_rna_mirna_protein_umap_df_survival)
rna_mirna_protein_pairwise_results_subtype$p.value

## Visualization
### 1) MOG
ggsurvplot(
  rna_mirna_protein_fit_subtype,
  data = merge_rna_mirna_protein_umap_df_survival,
  risk.table = TRUE,       
  pval = FALSE,            
  pval.coord = c(1, 0.1),   
  conf.int = FALSE,          
  xlab = "Time (days)",     
  legend.title = "Subtype",
  legend.labs = levels(merge_rna_mirna_protein_umap_df_survival$BRCA_Subtype_PAM50),
  palette = c("tomato", "yellow4", "turquoise3", "purple"),
  risk.table.height = 0.2,
  ggtheme = theme_minimal()
)

### 2) PAM50 Subtype 
ggsurvplot(
  rna_mirna_protein_fit_cluster,
  data = merge_rna_mirna_protein_umap_df_survival,
  risk.table = TRUE, 
  pval = FALSE,   
  pval.coord = c(1, 0.1),
  conf.int = FALSE,  
  xlab = "Time (days)", 
  legend.title = "Cluster",
  legend.labs = levels(merge_rna_mirna_protein_umap_df_survival$cluster),
  palette = c("#1B9E77", "#D95F02", "blue"),
  risk.table.height = 0.2,
  ggtheme = theme_minimal()
)

# Multi-omics signature analysis #
# MOG1-MOG2, MOG1-MOG3 proceed the same way

## Load the data
X <- list(RNAseq = rna_mirna_protein_RNAseq_matrix_top,
          miRNA = rna_mirna_protein_miRNA_matrix,
          Protein = rna_mirna_protein_Protein_matrix)
Y <- rna_mirna_protein_mofa_cluster$cluster
### For X, both miRNA and protein features were utilized. Since RNA-seq has a large number of features, the top 500 genes with the highest variability were filtered using the median absolute deviation.)

## Parameter chioce
### 1) Design matrix
design <- matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) <- 0

### 2) Number of components
diablo.tcga <- block.plsda(X, Y, ncomp = 5, design = design)

set.seed(123)
perf.diablo.tcga = perf(diablo.tcga, validation = 'Mfold', folds = 10, nrepeat = 10)

plot(perf.diablo.tcga)
perf.diablo.tcga$choice.ncomp$WeightedVote
ncomp <- perf.diablo.tcga$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
### ncomp = 3

### 3) Number of variables to select
set.seed(123)
test.keepX <- list(RNAseq = c(5:9, seq(10, 25, 5)),
                   miRNA = c(5:9, seq(10, 20, 2)),
                   Protein = c(seq(5, 25, 5)))

tune.diablo.tcga <- tune.block.splsda(X, Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              dist = "centroids.dist")

list.keepX <- tune.diablo.tcga$choice.keepX


## Final model
diablo.tcga_final <- block.splsda(X, Y, ncomp = ncomp, 
                            keepX = list.keepX, design = design)

## Visualization
circosPlot(diablo.tcga_final, cutoff = 0.8, line = TRUE, comp = c(1, 2, 3),
           color.blocks = c('steelblue', 'chartreuse3', 'tomato'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5,
          cex.labels=2)
network(diablo.tcga_final, blocks = c(1,2,3), cutoff = 0.8, size.node = 0.3)

# GSEA analysis #

## Differential expression analysis
rna_mirna_protein_RNAseq_binary <- ifelse(as.character(filter_rna_mirna_protein_mapping$cluster) == "MOG1", "MOG1", "MOG2")
rna_mirna_protein_RNAseq_group <- factor(rna_mirna_protein_RNAseq_binary, levels=c("MOG2","MOG1"))  # ref = Other
table(rna_mirna_protein_RNAseq_group)

rna_mirna_protein_RNAseq_design <-  model.matrix(~0 + rna_mirna_protein_RNAseq_group)
colnames(rna_mirna_protein_RNAseq_design ) <- levels(rna_mirna_protein_RNAseq_group)

rna_mirna_protein_RNAseq_contrast_matrix  <- makeContrasts(MOG2_vs_MOG1 = MOG2 - MOG1, levels=rna_mirna_protein_RNAseq_design)

rna_mirna_protein_RNAseq_fit <- lmFit(rna_mirna_protein_RNAseq_df, rna_mirna_protein_RNAseq_design)
rna_mirna_protein_RNAseq_fit2 <- contrasts.fit(rna_mirna_protein_RNAseq_fit, rna_mirna_protein_RNAseq_contrast_matrix)
# trend=TRUE is often helpful if mean-variance relationship exists; it makes eBayes more robust
rna_mirna_protein_RNAseq_fit2 <- eBayes(rna_mirna_protein_RNAseq_fit2, trend=TRUE)

rna_mirna_protein_RNAseq_res_all <- topTable(rna_mirna_protein_RNAseq_fit2, coef="MOG2_vs_MOG1", number=Inf, sort.by="P")

rna_mirna_protein_RNAseq_res_all <- rna_mirna_protein_RNAseq_res_all %>% tibble::rownames_to_column("gene")

## Creat GSEA object & MsigDB 
rna_mirna_protein_RNAseq_res_all_stats <- setNames(rna_mirna_protein_RNAseq_res_all$t,
                                                  rna_mirna_protein_RNAseq_res_all$gene)
rna_mirna_protein_RNAseq_res_all_stats <- sort(rna_mirna_protein_RNAseq_res_all_stats, decreasing = TRUE)

pathways_c5 <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>% 
  split(x = .$gene_symbol, f = .$gs_name)
pathways_c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>% 
  split(x = .$gene_symbol, f = .$gs_name)

rna_mirna_protein_fgseaRes_c5 <- fgseaMultilevel(
  pathways = pathways_c5,
  stats = rna_mirna_protein_RNAseq_res_all_stats,
  minSize = 15,
  maxSize = 500,
  eps = 0,
  nPermSimple = 100000
)

rna_mirna_protein_fgseaRes_c2 <- fgseaMultilevel(
  pathways = pathways_c2,
  stats = rna_mirna_protein_RNAseq_res_all_stats,
  minSize = 15,
  maxSize = 500,
  eps = 0,  # p-value 최소값 제한 해제
  nPermSimple = 100000
)

filter_rna_mirna_protein_fgseaRes_c5_MOG2_up <- filter_rna_mirna_protein_fgseaRes_c5 %>% arrange(padj) %>% filter(padj < 0.05) %>% filter(NES > 0) %>% head(10)
filter_rna_mirna_protein_fgseaRes_c5_MOG1_up <- filter_rna_mirna_protein_fgseaRes_c5 %>% arrange(padj) %>% filter(padj < 0.05) %>% filter(NES < 0) %>% head(10) %>% mutate(NES_abs = abs(NES))
filter_rna_mirna_protein_fgseaRes_c2_MOG2_up <- filter_rna_mirna_protein_fgseaRes_c2 %>% arrange(padj) %>% filter(padj < 0.05) %>% filter(NES > 0) %>% head(10)
filter_rna_mirna_protein_fgseaRes_c2_MOG1_up <- filter_rna_mirna_protein_fgseaRes_c2 %>% arrange(padj) %>% filter(padj < 0.05) %>% filter(NES < 0) %>% head(10) %>% mutate(NES_abs = abs(NES))

## Visualization
ggplot(filter_rna_mirna_protein_12_fgseaRes_c5_MOG2_up, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, color = padj)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "GO-BP: MOG2 upregulated",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    color = "p.adjust",
    size = "Gene Count"
  ) +
  theme_minimal(base_size = 12) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),  # y축 텍스트 설정
    strip.text = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 20, 10, 80),  # 여백 추가 (top, right, bottom, left)
  ) +
  theme(plot.title = element_text(hjust = 0.5))
### The rest is the same.
