# GSE5281_pipeline_fixed.R
# Robust pipeline for GSE5281 (EC or any brain region) - Ishaan Bhardwaj
# Set root project path & create folders first (change if needed)
proj_dir <- "C:/Users/ishaa/project/GSE5281_EC"
dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- file.path(proj_dir, "data"); dir.create(data_dir, showWarnings = FALSE)
out_dir  <- file.path(proj_dir, "output"); dir.create(out_dir, showWarnings = FALSE)
save_dir <- file.path(proj_dir, "saved_objects"); dir.create(save_dir, showWarnings = FALSE)
setwd(proj_dir)

# --------------------------
# 0. Package install / load
# --------------------------
cran_pkgs <- c("ggplot2","pheatmap","randomForest","caret","pROC","e1071","Boruta","dplyr","readr")
bioc_pkgs  <- c("GEOquery","Biobase","limma","hgu133plus2.db","org.Hs.eg.db","clusterProfiler","EnhancedVolcano","AnnotationDbi")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

# load
suppressPackageStartupMessages({
  lapply(c(bioc_pkgs, cran_pkgs), library, character.only = TRUE)
})

# set GEOquery options
options(download.file.method.GEOquery = "auto")
options(GEOquery.inmemory.gpl = FALSE)

# --------------------------
# 1. Load GSE/Series matrix
# --------------------------
series_matrix_path <- "C:/Users/ishaa/Downloads/GSE5281_series_matrix.txt.gz"
if (!file.exists(series_matrix_path)) stop("Series matrix not found at: ", series_matrix_path)

gset <- getGEO(filename = series_matrix_path, GSEMatrix = TRUE)

# getGEO sometimes returns a list or an ExpressionSet object
if (is(gset, "ExpressionSet")) {
  eset <- gset
} else if (is.list(gset) && length(gset) >= 1 && is(gset[[1]], "ExpressionSet")) {
  eset <- gset[[1]]
} else {
  stop("getGEO did not return an ExpressionSet. Inspect 'gset' manually.")
}

# quick checks
message("Annotation: ", annotation(eset))
message("Dimensions (genes x samples): ", paste(dim(exprs(eset)), collapse = " x "))

# --------------------------
# 2. Build robust phenotype table and region_base
# --------------------------
pheno <- pData(eset)   # data.frame

# region from title (first token), but many samples use names like EC_AFFECTED_3 - create region_base
pheno$region <- toupper(sub(" .*", "", pheno$title))
# collapse "_AFFECTED..." and suffixes to a base region
pheno$region_base <- toupper(gsub("(_AFFECTED.*|_AFFECTED_.*)$", "", pheno$region, ignore.case = TRUE))
# if the above didn't change anything, also try splitting at first underscore
pheno$region_base[is.na(pheno$region_base) | pheno$region_base == ""] <- toupper(sapply(strsplit(pheno$region, "_"), "[", 1))

# Robustly assign group labels
pheno$group <- NA_character_
# check title patterns
pheno$group[grepl("\\bcontrol\\b|\\bnormal\\b|non[- ]?demented|UNAFFECTED|UNAFF|CONTROL", pheno$title, ignore.case = TRUE)] <- "Control"
pheno$group[grepl("Alzheimer|Alz|AD|Affected|AFFECTED|dementia|AFFECT", pheno$title, ignore.case = TRUE)] <- "AD"

# fallback: check characteristic columns if exist
char_cols <- grep("^characteristics", colnames(pheno), ignore.case = TRUE, value = TRUE)
if (length(char_cols) > 0) {
  for (cc in char_cols) {
    pheno$group[is.na(pheno$group) & grepl("control|normal|non[- ]?demented|UNAFFECTED|UNAFF|CONTROL", pheno[[cc]], ignore.case = TRUE)] <- "Control"
    pheno$group[is.na(pheno$group) & grepl("Alzheimer|Alz|AD|Affected|AFFECTED", pheno[[cc]], ignore.case = TRUE)] <- "AD"
  }
}

# if still NA but region name contains AFFECTED -> mark AD, else if region name contains UNAFFECTED or CONTROL -> Control
pheno$group[is.na(pheno$group) & grepl("AFFECTED|AFF", pheno$region, ignore.case = TRUE)] <- "AD"
pheno$group[is.na(pheno$group) & grepl("UNAFFECTED|UNAFF|CONTROL", pheno$region, ignore.case = TRUE)] <- "Control"

# write back to eset
pData(eset) <- pheno

# --------------------------
# 3. Check which base regions have both groups
# --------------------------
tab <- table(pheno$region_base, pheno$group, useNA = "no")
good_regions <- rownames(tab)[apply(tab, 1, function(x) sum(x > 0) >= 2)]
message("Regions with both AD & Control (good choices): ", paste(good_regions, collapse = ", "))
message("Full table (region_base x group):")
print(tab)

# STOP here for manual choice if you want:
# choose one of the 'good_regions' below, e.g. "EC", "HIP", "MTG", "PC", "SFG", "VCX"
# If you want to continue automatically and EC is available, we choose EC; otherwise pick the first good one.
if ("EC" %in% good_regions) {
  chosen_region_base <- "EC"
} else if (length(good_regions) > 0) {
  chosen_region_base <- good_regions[1]
} else {
  stop("No region has both AD and Control. Inspect pheno and grouping logic.")
}
message("Selected region_base for analysis: ", chosen_region_base)

# --------------------------
# 4. Subset ExpressionSet for chosen region_base
# --------------------------
sel_idx <- which(pheno$region_base == chosen_region_base & !is.na(pheno$group))
if (length(sel_idx) < 4) stop("Too few samples in chosen region for meaningful analysis: ", length(sel_idx))

eset_sub <- eset[, sel_idx]
pheno_sub <- pData(eset_sub)
expr_sub <- exprs(eset_sub)   # genes x samples
message("Samples selected: ", ncol(expr_sub), " (", paste(table(pheno_sub$group), collapse = ", "), ")")

# --------------------------
# 5. Preprocess: log2 if needed and filtering
# --------------------------
if (max(expr_sub, na.rm = TRUE) > 100) {
  expr_sub <- log2(expr_sub + 1)
  message("Applied log2(x+1) transform")
}

# filter: remove lowest 20% by mean (keep top 80% by mean)
keep <- rowMeans(expr_sub, na.rm = TRUE) > quantile(rowMeans(expr_sub, na.rm = TRUE), 0.2)
exprf <- expr_sub[keep, ]
message("Filtered probes: kept ", nrow(exprf), " of ", nrow(expr_sub))

# --------------------------
# 6. Differential expression (limma) - AD vs Control
# --------------------------
group <- factor(pheno_sub$group)
if (length(levels(group)) < 2) stop("Only one group present - cannot run limma.")

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(exprf, design)
contrast.matrix <- makeContrasts(AD_vs_Control = AD - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, number = Inf, adjust.method = "BH")
res$probe <- rownames(res)

# map probe -> symbol (first mapping)
res$Symbol <- mapIds(hgu133plus2.db, keys = res$probe, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")

# save results
write.csv(res, file = file.path(out_dir, paste0("GSE5281_", chosen_region_base, "_limma_all.csv")), row.names = FALSE)
message("Limma results saved: ", file.path(out_dir, paste0("GSE5281_", chosen_region_base, "_limma_all.csv")))

# --------------------------
# 7. Plots: Volcano, PCA, Heatmap (top 50)
# --------------------------
# Volcano (EnhancedVolcano)
try({
  EnhancedVolcano(res,
                  lab = ifelse(is.na(res$Symbol), res$probe, res$Symbol),
                  x = 'logFC', y = 'P.Value',
                  pCutoff = 0.05, FCcutoff = 1,
                  title = paste0('GSE5281 ', chosen_region_base, ' Volcano'))
  ggsave(file.path(out_dir, paste0("GSE5281_", chosen_region_base, "_volcano.png")), width = 7, height = 7)
  message("Volcano saved.")
}, silent = TRUE)

# PCA (top 500 variable genes)
topVarGenes <- head(order(apply(exprf, 1, var), decreasing = TRUE), 500)
pca <- prcomp(t(exprf[topVarGenes, , drop = FALSE]), scale. = TRUE)
pcaDf <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = pheno_sub$group, sample = colnames(exprf))
png(file.path(out_dir, paste0("GSE5281_", chosen_region_base, "_PCA.png")), width = 900, height = 700)
print(ggplot(pcaDf, aes(PC1, PC2, color = group)) + geom_point(size=4) + theme_minimal() + labs(title=paste0("GSE5281 ", chosen_region_base," PCA (top 500 var genes)")))
dev.off()
message("PCA saved.")

# Heatmap top 50 DEGs (if available)
sig <- subset(res, adj.P.Val < 0.05 & abs(logFC) >= 1)
if (nrow(sig) >= 5) {
  top50 <- head(sig[order(sig$adj.P.Val), "probe"], 50)
  mat <- exprf[top50, , drop = FALSE]
  rowlabs <- ifelse(is.na(mapIds(hgu133plus2.db, keys=top50, column="SYMBOL", keytype="PROBEID", multiVals="first")),
                    top50,
                    mapIds(hgu133plus2.db, keys=top50, column="SYMBOL", keytype="PROBEID", multiVals="first"))
  rowlabs <- make.unique(rowlabs)
  rownames(mat) <- rowlabs
  mat_z <- t(scale(t(mat)))
  annotation_col <- data.frame(Group = factor(pheno_sub$group)); rownames(annotation_col) <- colnames(mat_z)
  pheatmap(mat_z, annotation_col = annotation_col, show_rownames = TRUE, filename = file.path(out_dir, paste0("GSE5281_", chosen_region_base, "_top50_heatmap.png")))
  message("Heatmap saved.")
} else {
  message("Too few significant DEGs to make top50 heatmap.")
}

# --------------------------
# 8. Prepare ML: features & data
# --------------------------
# choose top genes for ML: top 200 by adj.P.Val (from significant set); fallback use top variable genes if sig too few
if (nrow(sig) >= 10) {
  top_genes <- head(sig[order(sig$adj.P.Val), "probe"], 200)
} else {
  message("Not enough significant DEGs; using top variable genes for ML features.")
  top_genes <- head(order(apply(exprf,1,var), decreasing = TRUE), 200)
  top_genes <- names(top_genes)
}

X <- t(exprf[top_genes, , drop = FALSE])  # samples x genes
y <- factor(pheno_sub$group)

saveRDS(list(X = X, y = y, chosen_region = chosen_region_base), file = file.path(save_dir, "ml_data.rds"))
message("ML data saved: ", file.path(save_dir, "ml_data.rds"))

# --------------------------
# 9. Train simple models (RF and SVM) with caret (if enough samples)
# --------------------------
if (nrow(X) >= 6 && length(unique(y)) >= 2 && min(table(y)) >= 3) {
  set.seed(123)
  train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
  X_train <- X[train_idx, , drop=FALSE]; y_train <- y[train_idx]
  X_test  <- X[-train_idx, , drop=FALSE]; y_test  <- y[-train_idx]
  
  fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5,
                             classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = "final")
  
  # Random Forest
  set.seed(123)
  rf_model <- train(x = X_train, y = y_train, method = "rf", metric = "ROC",
                    trControl = fitControl, preProcess = c("center", "scale"), tuneLength = 3, ntree = 500)
  saveRDS(rf_model, file.path(save_dir, "rf_model.rds"))
  rf_pred_class <- predict(rf_model, X_test)
  rf_pred_prob  <- predict(rf_model, X_test, type = "prob")[,1]
  conf_rf <- confusionMatrix(rf_pred_class, y_test)
  write.csv(as.data.frame(conf_rf$table), file.path(out_dir, "rf_confusion_table.csv"), row.names = FALSE)
  message("RF trained and saved.")
  
  # SVM radial
  set.seed(123)
  svm_model <- train(x = X_train, y = y_train, method = "svmRadial", metric = "ROC",
                     trControl = fitControl, preProcess = c("center", "scale"), tuneLength = 6)
  saveRDS(svm_model, file.path(save_dir, "svm_model.rds"))
  message("SVM trained and saved.")
  
  # save feature importance (RF)
  rf_imp <- varImp(rf_model, scale = TRUE)
  rf_imp_df <- data.frame(Gene = rownames(rf_imp$importance), Importance = rf_imp$importance$Overall)
  rf_imp_df <- dplyr::arrange(rf_imp_df, desc(Importance))
  write.csv(rf_imp_df, file.path(out_dir, "rf_feature_importance.csv"), row.names = FALSE)
} else {
  message("Skipping ML training: insufficient samples per class (need at least ~3 per class in train/test).")
}

# --------------------------
# 10. Save key R objects / session info
# --------------------------
saveRDS(list(eset = eset, eset_sub = eset_sub, res = res), file = file.path(save_dir, paste0("GSE5281_", chosen_region_base, "_objects.rds")))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("✅ Pipeline complete. Outputs: ", out_dir, "  |  saved objects: ", save_dir)


