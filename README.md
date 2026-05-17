# RNA-Seq Differential Gene Expression and Machine Learning Analysis Pipeline

## Introduction

This repository contains a comprehensive R-based bioinformatics and machine learning workflow for RNA-Seq transcriptomics analysis. The pipeline integrates differential gene expression analysis, visualization, feature selection, and predictive machine learning modeling using multiple statistical and artificial intelligence approaches.

The workflow focuses on identifying significantly differentially expressed genes between IL6 and PBS experimental conditions and uses these genes for downstream machine learning classification.

The pipeline combines:

- Differential gene expression analysis
- RNA-Seq normalization
- Heatmap visualization
- Feature selection using LASSO
- Machine learning classification
- Model evaluation using ROC-AUC
- Variable importance analysis

---

# Objectives

The primary objectives of this project are:

- Import and preprocess RNA-Seq expression data
- Perform normalization and differential expression analysis
- Identify significant differentially expressed genes (DEGs)
- Visualize gene expression patterns
- Apply feature selection methods
- Build predictive machine learning models
- Evaluate classification performance
- Identify important biomarker genes

---

# Libraries Used

The following R libraries are required:

```r
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
