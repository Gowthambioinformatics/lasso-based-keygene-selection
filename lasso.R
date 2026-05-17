# Load required libraries
library(limma)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(caret)
library(randomForest)
library(e1071)
library(pROC)
library(xgboost)
library(glmnet)
library(keras)


rpkm_file <- "...Own your path..."
rpkm_data <- read.csv(rpkm_file, sep="\t", row.names=1)


il6_columns <- grep("_il6_", colnames(rpkm_data), value=TRUE)
pbs_columns <- grep("_pbs_", colnames(rpkm_data), value=TRUE)
selected_data <- rpkm_data[, c(il6_columns, pbs_columns)]

dge <- DGEList(counts = selected_data)
dge <- calcNormFactors(dge)
v <- voom(dge)


group <- factor(c(rep("IL6", length(il6_columns)), rep("PBS", length(pbs_columns))))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)


fit <- lmFit(v, design)

contrast_matrix <- makeContrasts(IL6_vs_PBS = IL6 - PBS, levels=design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)


results <- topTable(fit2, adjust="fdr", number=nrow(rpkm_data))


significant_genes <- subset(results, adj.P.Val < 0.05)


write.csv(significant_genes, "...Own your path...")

n_sig_genes <- sum(significant_genes$adj.P.Val < 0.05)
cat("Number of significant genes:", n_sig_genes, "\n")

significant_gene_ids <- rownames(significant_genes)
expr_data <- rpkm_data[significant_gene_ids, ]


normalized_expr_data <- t(scale(t(expr_data)))


short_colnames <- sapply(strsplit(colnames(normalized_expr_data), "_"), function(x) x[1])
colnames(normalized_expr_data) <- short_colnames


heatmap <- pheatmap(normalized_expr_data, scale="row", clustering_distance_rows="correlation",
                    clustering_distance_cols="euclidean", display_numbers=TRUE,
                    fontsize_row=6, fontsize_col=8, angle_col=45)


heatmap_file <- "...Own your path..."
png(heatmap_file, width=800, height=800, res=150)
print(heatmap)
dev.off()

cat("Heatmap of gene expression saved as", heatmap_file, "\n")


deg_results_file <- "...Own your path..."
deg_results <- read.csv(deg_results_file, row.names = 1)


significant_genes <- rownames(subset(deg_results, adj.P.Val < 0.05))
selected_data <- rpkm_data[significant_genes, ]


selected_data_t <- t(selected_data)


labels <- factor(c(rep("IL6", length(il6_columns)), rep("PBS", length(pbs_columns))))


print(table(labels))


ml_data <- data.frame(selected_data_t)
ml_data$Group <- labels


preProcValues <- preProcess(ml_data[, -ncol(ml_data)], method = c("center", "scale"))
ml_data_preprocessed <- predict(preProcValues, ml_data)


pca <- prcomp(ml_data_preprocessed[, -ncol(ml_data_preprocessed)], center = TRUE, scale. = TRUE)


num_pcs <- min(ncol(ml_data_preprocessed) - 1, nrow(pca$x))


pca_data <- data.frame(pca$x[, 1:num_pcs])


pca_data$Group <- ml_data_preprocessed$Group

lasso_model <- cv.glmnet(as.matrix(ml_data_preprocessed[, -ncol(ml_data_preprocessed)]), ml_data_preprocessed$Group, family="binomial", alpha=1)
lasso_coefficients <- as.matrix(coef(lasso_model, s = "lambda.min"))
selected_genes_lasso <- rownames(lasso_coefficients)[lasso_coefficients != 0]
selected_genes_lasso <- selected_genes_lasso[selected_genes_lasso != "(Intercept)"]


cat("Selected genes by LASSO:", selected_genes_lasso, "\n")


lasso_selected_data <- ml_data_preprocessed[, c(selected_genes_lasso, "Group")]


train_control <- trainControl(method = "LOOCV", 
                              classProbs = TRUE, 
                              summaryFunction = twoClassSummary)


set.seed(123)
models <- list()


models$LogisticRegression <- train(Group ~ ., data = lasso_selected_data, 
                                   method = "glm", 
                                   trControl = train_control, 
                                   metric = "ROC")


models$RandomForest <- train(Group ~ ., data = lasso_selected_data, 
                             method = "rf", 
                             trControl = train_control, 
                             metric = "ROC")


models$SVM <- train(Group ~ ., data = lasso_selected_data, 
                    method = "svmRadial", 
                    trControl = train_control, 
                    metric = "ROC")


models$XGBoost <- train(Group ~ ., data = lasso_selected_data, 
                        method = "xgbTree", 
                        trControl = train_control, 
                        metric = "ROC")


results <- lapply(models, function(model) {
  pred <- predict(model, lasso_selected_data)
  prob <- predict(model, lasso_selected_data, type = "prob")
  conf_matrix <- confusionMatrix(pred, lasso_selected_data$Group)
  roc_curve <- roc(lasso_selected_data$Group, prob[, "IL6"])
  list(conf_matrix = conf_matrix, roc_curve = roc_curve)
})

for (model_name in names(results)) {
  cat("\nModel:", model_name, "\n")
  print(results[[model_name]]$conf_matrix)
  cat("AUC:", auc(results[[model_name]]$roc_curve), "\n")
}


importance <- lapply(models, varImp, scale = FALSE)
for (model_name in names(importance)) {
  cat("\nVariable Importance for Model:", model_name, "\n")
  print(importance[[model_name]])
  plot(importance[[model_name]])
}


best_model_name <- names(models)[which.max(sapply(results, function(x) auc(x$roc_curve)))]
best_model <- models[[best_model_name]]
saveRDS(best_model, "best_ml_model.rds")

cat("Best model saved as best_ml_model.rds\n")
