####################################
#
# created on 2024-12-16
####################################

# load libraries
library(tidyverse)
library(rlang)
library(openxlsx)
library(readxl)
library(tableone)
library(dplyr)
library(foreign)
library(e1071)
library(pROC)
library(ModelMetrics)
library(caret)
library(gmodels)
library(janitor)
library(magrittr)
library(gridExtra)
library(robustbase)
library(StatMatch)
library(RColorBrewer)
library(colourvalues)
library(ggfortify)
library(ggforce)
library(ggsignif)
library(readr)
library(vroom)
library(data.table)
library(psych)
library(survival)
library(survminer)
library(patchwork)
library(rms)
library(glmnet)
library(ROCR)
library(ggplot2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(sva)
library(DESeq2)
library(edgeR)
library(factoextra)
library(cluster)
library(tximport)
library(jsonlite)
library(txdbmaker)
library(GenomicFeatures)
library(Biostrings)
library(tidyr)
library(limma)
library(rstatix)

# Read metadata files
GSE63042 <- read.csv("GSE63042_PEP_raw.csv", check.names = FALSE)
GSE185263 <- read.csv("GSE185263_PEP_raw.csv", check.names = FALSE)
GSE131411 <- read.csv("GSE131411_metadata.csv", check.names = FALSE)
GSE222393 <- read.csv("GSE222393_metadata.csv", check.names = FALSE)

# tximport
tx.exp <- tximport(salmon_files, type = "salmon", txOut = TRUE, tx2gene = tx2gene, ignoreTxVersion = TRUE)
# gene
gene.exp <- summarizeToGene(tx.exp, tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

# preprocess
## data integration
dge <- DGEList(counts = count1, group = coldata1$group)
keep <- filterByExpr(dge)
dge_filtered <- dge[keep, , keep.lib.sizes=FALSE]
# ComBat-seq
batch_corrected_counts <- ComBat_seq(
    counts = dge_filtered$counts,
    batch = coldata1$paper,
    group = coldata1$group
)

dge_corrected <- DGEList(counts = batch_corrected_counts, 
                        group = coldata1$group)
dge_corrected <- calcNormFactors(dge_corrected)

design_corrected <- model.matrix(~ group, data = coldata1)

# voom
v <- voom(dge_corrected, design_corrected, plot = TRUE)
v_batch <- v

dge_before <- calcNormFactors(dge_filtered)
design_before <- model.matrix(~ group + paper, data = coldata1)
v_before <- voom(dge_before, design_before, plot = FALSE)

corrected_counts_df <- as.data.frame(v_batch$E)

### clustering
corrected_counts_df %>% 
  as.matrix() -> corrected_counts_df_mat
TITLE="km_output"
SEED = 123
res_km_euc <- ConsensusClusterPlus::ConsensusClusterPlus(d = corrected_counts_df_mat, 
                                                             maxK = 7,
                                                             reps = 100, 
                                                             pItem=0.8, 
                                                             pFeature=1,         
                                                             clusterAlg="km",   
                                                             distance="euclidean",
                                                             seed = SEED,
                                                             title=TITLE, plot ="png")
# Get output of clusters for each pt
km_euc_df <- corrected_counts_df %>% 
                    t() %>%
                    as.data.frame() %>%
                    rownames(.) %>% 
                    as.data.frame() %>%
                    dplyr::select("public_id" = 1) %>% 
      mutate(km_cluster_k3 = res_km_euc[[3]][["consensusClass"]]) 

corrected_counts_df %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "public_id") %>% 
  inner_join(km_euc_df, by = "public_id") -> merged_batch_km_df

# To assess
icl_km_euc_norm = ConsensusClusterPlus::calcICL(res_km_euc,
                                    title=TITLE, plot="png")
icl_km_euc_norm[["clusterConsensus"]]
icl_km_euc_norm[["itemConsensus"]][1:5,]
km_euc_cluster_consensus_df <- icl_km_euc_norm[["clusterConsensus"]] %>% 
                        as.data.frame() %>% 
                        dplyr::filter(k %in%  c("2", "3", "4", "5")) %>% 
                        dplyr::select(k = k, Cluster = cluster, "Cluster consensus value" = clusterConsensus)

km_euc_cluster_consensus_df_3 <- km_euc_cluster_consensus_df %>% 
                          dplyr::filter(k %in%  c("3")) %>% 
                          dplyr::select(Cluster, "Cluster consensus value") %>% 
                          mutate(Cluster = as.factor(Cluster))

### Mortality
meta_batch_km_df$death_status <- ifelse(meta_batch_km_df$group == "Death", 1, 0)

lowest_mortality_cluster <- meta_batch_km_df %>%
  group_by(km_cluster_k3) %>%
  summarise(mortality_rate = mean(death_status)) %>%
  filter(mortality_rate == min(mortality_rate)) %>%
  pull(km_cluster_k3)

# GLM
glm(death_status ~ relevel(factor(km_cluster_k3), ref = lowest_mortality_cluster),
    family = binomial(link = "logit"),
    data = meta_batch_km_df) %>% 
  broom::tidy(exp = TRUE, conf.int = TRUE) %>% 
  dplyr::filter(!grepl("Intercept", term)) %>% 
  dplyr::select(term, estimate, conf.low, conf.high, p.value, std.error, statistic)

### Visualization
# PCA
pca_result <- prcomp(t(corrected_counts_df_mat), scale. = TRUE)

# DEG
### 1cluster vs others
norm_batch_seq_df <- meta_batch_km_df %>% 
                      column_to_rownames(var = "public_id") %>% 
                      dplyr::select(-c(group, 
                                      km_cluster_k3)) %>% 
                      t() %>% 
                      as.data.frame()

km_cluster_k3 <- meta_batch_km_df$km_cluster_k3

km_cluster_k3_1 <- ifelse(km_cluster_k3 == 1, "Cluster1", "Others")

design <- model.matrix(~0+km_cluster_k3_1)
colnames(design) <- levels(factor(km_cluster_k3_1))

fit <- lmFit(norm_batch_seq_df, design)

contrast.matrix <- makeContrasts(Cluster1 - Others, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 

# DEG
km_cluster_k3_1_res <- topTable(fit2, coef=1, n=Inf)

ensembl_ids <- rownames(km_cluster_k3_1_res)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

#### GSEA
# HALLMARK
m_df <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(m_df$gene_symbol, m_df$gs_name)
fgsea(pathways = pathways, 
                  stats = ranks,
                  minSize = 15,
                  maxSize = 500,
                  nPerm = 1000)
# GO
gseGO(geneList = res_km_cluster_k3_1_gsea1, 
                    OrgDb        = org.Hs.eg.db,
                    ont          = "BP",
                    minGSSize    = 100,
                    maxGSSize = 500, 
                    pAdjustMethod = "fdr", 
                    verbose = TRUE, 
                    by = "fgsea",
                    pvalueCutoff  = 0.2
                    )

# Classifier LASSO
# Choose Top200 gene
selected_genes_cluster1 <- km_cluster_k3_1_res3 %>%
    arrange(adj.P.Val) %>%
    head(200) %>%
    pull(gene_symbol)
selected_genes_cluster2 <- km_cluster_k3_2_res3 %>%
    arrange(adj.P.Val) %>%
    head(200) %>%
    pull(gene_symbol)
selected_genes_cluster3 <- km_cluster_k3_3_res3 %>%
    arrange(adj.P.Val) %>%
    head(200) %>%
    pull(gene_symbol)

all_selected_genes <- unique(c(selected_genes_cluster1, 
                             selected_genes_cluster2, 
                             selected_genes_cluster3))
# LASSO
select_genes_multinomial_lasso <- function(expression_data, cluster_labels, target_n_genes) {

  x <- as.matrix(t(expression_data))
  y <- factor(cluster_labels)

  x <- scale(x)

  set.seed(42)
  cv_fit <- cv.glmnet(
    x = x,
    y = y,
    family = "multinomial",
    alpha = 1,
    nfolds = 10,
    lambda = exp(seq(log(1e-4), log(10), length.out = 100)) 
  )
  best_lambda <- cv_fit$lambda.min
  coef_matrix <- coef(cv_fit$glmnet.fit, s = best_lambda)
  importance_matrix <- matrix(0, nrow = length(rownames(expression_data)), ncol = 1)
  rownames(importance_matrix) <- rownames(expression_data)

  for(i in seq_along(coef_matrix)) {
    class_coef <- abs(as.matrix(coef_matrix[[i]]))[-1, , drop = FALSE]  
    importance_matrix <- importance_matrix + class_coef
  }

  importance_matrix <- importance_matrix / length(coef_matrix)

  importance_df <- data.frame(
    gene = rownames(importance_matrix),
    importance = as.vector(importance_matrix),
    stringsAsFactors = FALSE
  )

  importance_df <- importance_df[order(-importance_df$importance), ]
  importance_over0_df <- importance_df[importance_df$importance > 0, ]
  selected_genes <- importance_df$gene[1:target_n_genes]
  selected_importance <- importance_df$importance[1:target_n_genes]
  names(selected_importance) <- selected_genes
  return(list(
    genes = selected_genes,
    lambda = best_lambda,
    importance = selected_importance,
    all_importance = importance_over0_df
  ))
}

selected_genes_multinomial <- select_genes_multinomial_lasso(
  expression_data = train_expression_symbol,
  cluster_labels = train_coldata$cluster,
  target_n_genes = 14
)

evaluate_single_run <- function(expression_data, selected_genes, endotype_labels) {
  expression_subset <- as.matrix(expression_data[selected_genes, ])
  model_data <- data.frame(t(expression_subset))
  model_data$endotype <- factor(paste0("Cluster_", endotype_labels))
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    savePredictions = TRUE
  ) 
  final_model <- train(
    endotype ~ .,
    data = model_data,
    method = "glmnet",
    trControl = ctrl,
    metric = "Accuracy"
  )  
  cv_predictions <- final_model$pred %>%
    group_by(rowIndex) %>%
    slice_tail(n = 1) %>%
    ungroup()
  auc_results <- lapply(levels(model_data$endotype), function(endo) {
    binary_labels <- ifelse(cv_predictions$obs == endo, 1, 0)
    prob_positive <- cv_predictions[[endo]]
    roc_obj <- pROC::roc(binary_labels, prob_positive, quiet = TRUE)
    data.frame(
      endotype = endo,
      auc = as.numeric(pROC::auc(roc_obj)),
      ci_lower = as.numeric(pROC::ci(roc_obj)[1]),
      ci_upper = as.numeric(pROC::ci(roc_obj)[3])
    )
  }) %>% bind_rows()  
  conf_mat <- confusionMatrix(cv_predictions$pred, cv_predictions$obs)  
  list(
    model = final_model,
    performance = final_model$results,
    auc_by_endotype = auc_results,
    mean_auc = mean(auc_results$auc),
    confusion_matrix = conf_mat,
    cv_predictions = cv_predictions,
    model_data = model_data 
  )
}
evaluate_endotype_model <- function(expression_data, selected_genes, endotype_labels, seed = 123) {
  set.seed(seed)
  result <- evaluate_single_run(
    expression_data = expression_data,
    selected_genes = selected_genes,
    endotype_labels = endotype_labels
  )
  # AUC
  cat("\nAUC by Endotype:\n")
  print(result$auc_by_endotype %>%
    mutate(across(where(is.numeric), ~round(. * 100, 1))))
  
  cat(sprintf("\nMean AUC: %.1f%%\n", result$mean_auc * 100))
  
  cat("\nConfusion Matrix:\n")
  print(result$confusion_matrix$table)
  
  cat("\nModel Accuracy:\n")
  print(sprintf("%.2f%%", result$confusion_matrix$overall["Accuracy"] * 100))
  
  return(result)
}

model_eval <- evaluate_endotype_model(
  expression_data = train_expression_symbol,
  selected_genes = selected_genes_multinomial$genes,
  endotype_labels = train_coldata$cluster
)

# microarray GSE236713
RG.bg <- backgroundCorrect(RG, method="normexp")
RG.bg$E <- RG.bg$E + min(RG.bg$E[RG.bg$E > 0])  

RG.log <- log2(RG.bg$E)

normalized.data <- normalizeBetweenArrays(RG.log, method="quantile")  

## filtering
is.control <- RG$genes$ControlType == 1 | RG$genes$ControlType == -1
normalized.filtered <- normalized.data[!is.control,]

ave.expr <- rowMeans(normalized.filtered)

expr.threshold <- quantile(ave.expr, 0.1)  
min.samples <- ncol(normalized.filtered) * 0.2  

is.expressed <- ave.expr > expr.threshold &  
               rowSums(normalized.filtered > expr.threshold) >= min.samples  

normalized.filtered <- normalized.filtered[is.expressed,]
RG.filtered <- RG.bg[!is.control,][is.expressed,]  

# gene classifier
selected_genes <- selected_genes_multinomial$genes
matched_genes <- intersect(selected_genes, rownames(expression_df_final))

selected_df <- expression_df_final[matched_genes, ]

selected_df_mat <- as.matrix(selected_df)

# Clustering
TITLE <- "km_selected_genes_microarray"
SEED <- 123
res_km_euc_selected <- ConsensusClusterPlus(d = selected_df_mat,
                                           maxK = 7,
                                           reps = 100,
                                           pItem = 0.87,
                                           pFeature = 1,
                                           clusterAlg = "pam",
                                           distance = "euclidean",
                                           seed = SEED,
                                           title = TITLE,
                                           plot = "png")

km_selected_df <- selected_df %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "geo_accession") %>%
  dplyr::select(geo_accession) %>%
  mutate(
    km_selected_cluster_k2 = res_km_euc_selected[[2]][["consensusClass"]],
    km_selected_cluster_k3 = res_km_euc_selected[[3]][["consensusClass"]],
    km_selected_cluster_k4 = res_km_euc_selected[[4]][["consensusClass"]],
    km_selected_cluster_k5 = res_km_euc_selected[[5]][["consensusClass"]],
    km_selected_cluster_k6 = res_km_euc_selected[[6]][["consensusClass"]]
  )

merged_targets_km_df <- targets %>%
  dplyr::select(geo_accession, Group) %>%
  inner_join(km_selected_df, by = "geo_accession")

# ICL
icl_km_euc_selected <- ConsensusClusterPlus::calcICL(res_km_euc_selected,
                                                    title = "km_selected_genes_microarray",
                                                    plot = "png")

km_selected_cluster_consensus_df <- icl_km_euc_selected[["clusterConsensus"]] %>% 
  as.data.frame() %>% 
  dplyr::filter(k %in% c("2", "3", "4", "5")) %>% 
  dplyr::select(k = k, Cluster = cluster, "Cluster consensus value" = clusterConsensus)

km_selected_cluster_consensus_df_2 <- km_selected_cluster_consensus_df %>% 
  dplyr::filter(k %in% c("2")) %>% 
  dplyr::select(Cluster, "Cluster consensus value") %>% 
  mutate(Cluster = as.factor(Cluster))
km_selected_cluster_consensus_df_3 <- km_selected_cluster_consensus_df %>% 
  dplyr::filter(k %in% c("3")) %>% 
  dplyr::select(Cluster, "Cluster consensus value") %>% 
  mutate(Cluster = as.factor(Cluster))
km_selected_cluster_consensus_df_4 <- km_selected_cluster_consensus_df %>% 
  dplyr::filter(k %in% c("4"))
km_selected_cluster_consensus_df_5 <- km_selected_cluster_consensus_df %>% 
  dplyr::filter(k %in% c("5"))

# mortality
calculate_mortality_rate <- function(df, group_col, outcome_col) {
  df %>%
    group_by(!!sym(group_col)) %>%
    summarise(
      Total = n(),
      Deaths = sum(!!sym(outcome_col) == "Death"),
      MortalityRate = Deaths / Total,
      MortalityPercentage = sprintf("%.2f%%", MortalityRate * 100)
    ) %>%
    arrange(!!sym(group_col))
}

mortality_rates_k3 <- calculate_mortality_rate(merged_targets_km_df, "km_selected_cluster_k3", "Group")

## PCA
pca_result <- prcomp(t(selected_df_mat), scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$cluster <- merged_targets_km_df$km_selected_cluster_k3
pca_df <- pca_df %>%
  left_join(mortality_rates_k3, by = c("cluster" = "km_selected_cluster_k3"))

cluster_centroids <- pca_df %>%
  group_by(cluster) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

variance <- pca_result$sdev^2
variance_prop <- variance / sum(variance) * 100
variance_prop_pc1 <- round(variance_prop[1], 1)
variance_prop_pc2 <- round(variance_prop[2], 1)

ggplot(pca_df, aes(x = PC1, y = PC2, color = factor(cluster))) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_line(color = "grey95"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(
    x = paste0("PC1 (", variance_prop_pc1, "%)"),
    y = paste0("PC2 (", variance_prop_pc2, "%)"),
    color = "Cluster"
  ) +
  scale_color_manual(values = c("red", "blue", "#006400"))

## GSE236713
### DEG
norm_microarray_df <- expression_df_final %>% 
  as.data.frame()

km_selected_cluster_k3 <- merged_targets_km_df$km_selected_cluster_k3

km_selected_cluster_k3_1 <- ifelse(km_selected_cluster_k3 == 1, "Cluster1", "Others")

design <- model.matrix(~0+km_selected_cluster_k3_1)
colnames(design) <- levels(factor(km_selected_cluster_k3_1))

fit <- lmFit(norm_microarray_df, design)

contrast.matrix <- makeContrasts(Cluster1 - Others, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

km_selected_cluster_k3_1_res <- topTable(fit2, coef=1, n=Inf)

km_selected_cluster_k3_1_res <- km_selected_cluster_k3_1_res %>%
  as.data.frame() %>%
  rownames_to_column(var = "ProbeName") %>%
  dplyr::arrange(adj.P.Val, P.Value)


### GSEA

ranks <- km_selected_cluster_k3_1_res$logFC
names(ranks) <- km_selected_cluster_k3_1_res$ProbeName
ranks <- sort(ranks, decreasing = TRUE)

m_df <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(m_df$gene_symbol, m_df$gs_name)

set.seed(123)
fgsea(pathways = pathways,
                                   stats = ranks,
                                   minSize = 15,
                                   maxSize = 500,
                                   nPerm = 1000)
### GO
mapIds(org.Hs.eg.db, 
       keys = km_selected_cluster_k3_1_res$ProbeName,
       column = "ENTREZID", 
       keytype = "SYMBOL") %>% 
  as.data.frame() %>%
  dplyr::select("ENTREZID" = 1) %>% 
  rownames_to_column(var = "ProbeName") %>% 
  inner_join(km_selected_cluster_k3_1_res, by = "ProbeName") %>% 
  dplyr::select("ENTREZID", "logFC") %>% 
  mutate_at(vars(1), as.character) %>% 
  distinct(ENTREZID, .keep_all = T) -> gsea_cluster1_data

gsea_cluster1_ranked <- as.numeric(gsea_cluster1_data[,2])
names(gsea_cluster1_ranked) = as.character(gsea_cluster1_data[, 1])
gsea_cluster1_ranked = sort(gsea_cluster1_ranked, decreasing = TRUE)

gseGO(geneList = gsea_cluster1_ranked, 
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         minGSSize = 100,
                         maxGSSize = 500, 
                         pAdjustMethod = "fdr", 
                         verbose = TRUE, 
                         by = "fgsea",
                         pvalueCutoff = 0.2,
                         scoreType = "pos")

## Clinical data
calculate_odds_ratios <- function(merged_data, categorical_cols) {
  odds_ratios <- list()
  for(col in categorical_cols) {
    pairs <- combn(unique(merged_data$Predicted_Endotype), 2)
    for(i in 1:ncol(pairs)) {
      cluster1 <- pairs[1,i]
      cluster2 <- pairs[2,i]
      for(category in unique(merged_data[[col]][!is.na(merged_data[[col]])])) {
        contingency <- table(
          merged_data[[col]][merged_data$Predicted_Endotype %in% c(cluster1, cluster2)] == category,
          merged_data$Predicted_Endotype[merged_data$Predicted_Endotype %in% c(cluster1, cluster2)]
        )
        if(nrow(contingency) == 2 && ncol(contingency) == 2) {
          or <- fisher.test(contingency)
          odds_ratios[[paste(col, category, "Cluster", cluster1, "vs", cluster2)]] <- list(
            odds_ratio = or$estimate,
            ci_lower = or$conf.int[1],
            ci_upper = or$conf.int[2],
            p_value = or$p.value
          )
        }
      }
    }
  }

  or_df <- do.call(rbind, lapply(names(odds_ratios), function(name) {
    data.frame(
      comparison = name,
      odds_ratio = odds_ratios[[name]]$odds_ratio,
      ci_lower = odds_ratios[[name]]$ci_lower,
      ci_upper = odds_ratios[[name]]$ci_upper,
      p_value = odds_ratios[[name]]$p_value
    )
  }))
  return(or_df)
}

clinical_analysis <- function(targets, res_km_euc_selected, k = 3) {
  cluster_assignments <- data.frame(
    geo_accession = rownames(t(selected_df)),
    Predicted_Endotype = res_km_euc_selected[[k]][["consensusClass"]]
  ) 
  merged_data <- targets %>%
    inner_join(cluster_assignments, by = "geo_accession")  
  numeric_cols <- c(
    "Age", 
    "APACHEIIScoreonadmission",
    "SOFA Overall Score Day 1", 
    "SOFA Overall Score Day 2", 
    "SOFA Overall Score Day 5", 
    "SOFA Overall Score Day Discharge",
    "White Blood Cells D1",
    "White Blood Cells D2",
    "White Blood Cells D5",
    "White Blood Cells Discharge", 
    "Neutrophils D1",
    "Neutrophils D2",
    "Neutrophils D5",
    "Neutrophils Discharge",
    "Lymphocytes D1",
    "Lymphocytes D2",
    "Lymphocytes D5",
    "Lymphocytes Discharge",
    "Basophils D1",
    "Basophils D2", 
    "Basophils D5",
    "Basophils Discharge", 
    "Platelets D1",
    "Platelets D2",
    "Platelets D5",
    "Platelets Discharge",
    "CRP D1",
    "CRP D2",
    "CRP D5...35",
    "CRP D5...41"
  )  
  categorical_cols <- c(
    "Group", 
    "Disease sub-type", 
    "Sex",
    "Ethnicity",
    "Outofhospitalcardiacarrest"
  )

  total_n_summary <- merged_data %>%
    group_by(Predicted_Endotype) %>%
    summarise(Total_N = n())

  numeric_summary <- merged_data %>%
    group_by(Predicted_Endotype) %>%
    summarise(across(all_of(numeric_cols), function(x) {
      n_missing <- sum(is.na(x))
      total_n <- n()
      valid_data <- x[!is.na(x)]     
      if(length(valid_data) > 0) {
        sprintf("%0.1f (%0.1f-%0.1f) [Missing: %d (%0.1f%%)]",
                median(valid_data),
                quantile(valid_data, 0.25),
                quantile(valid_data, 0.75),
                n_missing,
                (n_missing/total_n) * 100)
      } else {
        "No valid data"
      }
    }))

  categorical_summary <- merged_data %>%
    group_by(Predicted_Endotype) %>%
    summarise(across(all_of(categorical_cols), function(x) {
      total_n <- n()
      n_missing <- sum(is.na(x))
      valid_data <- table(x[!is.na(x)])      
      if(length(valid_data) > 0) {
        valid_props <- prop.table(valid_data) * 100
        category_stats <- paste(names(valid_data), 
                              sprintf("%d (%0.1f%%)", 
                                     valid_data, 
                                     valid_props),
                              collapse = ", ")       
        sprintf("%s [Missing: %d (%0.1f%%)]",
                category_stats,
                n_missing,
                (n_missing/total_n) * 100)
      } else {
        "No valid data"
      }
    }))

  final_summary <- total_n_summary %>%
    left_join(numeric_summary, by = "Predicted_Endotype") %>%
    left_join(categorical_summary, by = "Predicted_Endotype") %>%
    rename(Cluster = Predicted_Endotype)

  final_summary <- final_summary %>%
    select(Cluster, Total_N, everything())

  kw_tests <- sapply(numeric_cols, function(col) {
    tryCatch({
      test_data <- merged_data[!is.na(merged_data[[col]]), ]
      if (nrow(test_data) > 0) {
        test <- kruskal.test(as.formula(sprintf("`%s` ~ Predicted_Endotype", col)), data = test_data)
        test$p.value
      } else {
        NA
      }
    }, error = function(e) NA)
  })
  chi_tests <- sapply(categorical_cols, function(col) {
    tryCatch({
      test_data <- merged_data[!is.na(merged_data[[col]]), ]
      if (nrow(test_data) > 0) {
        tbl <- table(test_data[[col]], test_data$Predicted_Endotype)
        if (all(dim(tbl) > 1)) {
          test <- chisq.test(tbl, simulate.p.value = TRUE)
          test$p.value
        } else {
          NA
        }
      } else {
        NA
      }
    }, error = function(e) NA)
  })


  p_values <- c(kw_tests, chi_tests)
  names(p_values) <- c(numeric_cols, categorical_cols)

  return(list(
    summary_table = final_summary,
    p_values = p_values,
    odds_ratios = odds_ratios_results,
    merged_data = merged_data 
  ))
}
