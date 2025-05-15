Overview

This repository contains R code for the analysis of gene expression data to identify molecular endotypes in sepsis patients. The analysis combines data from multiple public datasets, performs batch correction, unsupervised clustering, and develops a molecular classifier to predict patient outcomes.

Features

    Data integration across multiple gene expression datasets (RNA-seq and microarray)
    Batch effect correction using ComBat-seq
    Unsupervised consensus clustering to identify patient endotypes
    Differential gene expression analysis between endotypes
    Gene set enrichment analysis (GSEA) to identify biological pathways
    Development of a LASSO-based classifier for endotype prediction
    Validation in an independent microarray dataset
    Clinical phenotype association with identified endotypes

Workflow

    Data Import and Preprocessing
        Read RNA-seq counts from Salmon files using tximport
        Summarize transcript-level to gene-level expression
        Filter low-expression genes
        Batch correction with ComBat-seq
        Normalization using voom transformation
    Unsupervised Clustering
        Consensus clustering using ConsensusClusterPlus
        Determination of optimal cluster number
        Assessment of cluster stability
    Mortality Analysis
        Comparison of mortality rates between endotypes
        Logistic regression analysis
    Differential Expression Analysis
        Identification of differentially expressed genes between endotypes
        Visualization with PCA
    Pathway Analysis
        Gene set enrichment analysis (GSEA) using MSigDB Hallmark pathways
        Gene Ontology enrichment analysis
    Classifier Development
        Selection of top differentially expressed genes
        LASSO-based regularization for feature selection
        Multinomial model development
        Cross-validation for performance assessment
    External Validation
        Preprocessing of microarray data (GSE236713)
        Application of the gene classifier
        Cluster analysis in validation cohort
        Comparison of mortality outcomes
    Clinical Association Analysis
        Calculation of odds ratios for categorical variables
        Comparison of clinical features between endotypes
        Statistical testing (Kruskal-Wallis, Chi-square)

Usage

The script is designed to be run in R environment with all required packages installed. Dataset-specific paths and parameters may need to be adjusted before execution.
