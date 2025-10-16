# 🧬 Endo-ML: Machine Learning for Endometriosis Biomarker Discovery

A multimodal bioinformatics and machine learning project for identifying biomarkers and predicting endometriosis status using RNA-Seq and DNA methylation data.

---

## 🚀 Overview & Motivation

Endometriosis is an under-researched gynecological disorder affecting an estimated 5–10 % of women.  
Current diagnostics are invasive and often delayed.  

This project applies **bioinformatics and ML techniques** to support **non-invasive diagnosis** and **molecular understanding** of endometriosis.

This study combines: 

- Differential gene expression analysis (RNA-Seq)
- Differential methylation analysis (e.g., MBD-Seq)
- Multi-omics integration
- Supervised ML models for disease prediction (**Two machine learning architectures**)
- Feature selection and biological interpretation (GO / pathway enrichment)

This project demonstrates **omics data processing, differential analysis, functional enrichment, feature engineering, model development (ensemble ML modeling), and biological reasoning**.


---

## 🧠 Workflow Summary

1. **RNA-Seq analysis (DESeq2)**  
   - Filtering low-count genes, normalization (median of ratios), log2 fold-change shrinkage (apeglm)  
   - 100 significant DEGs (padj < 0.05): 83 downregulated, 17 upregulated  
   - Top genes: `COL1A1`, `TENM2`, `IGF2`, `KLF2P1`

2. **Methylation analysis (edgeR)**  
   - Normalization via TMM, GLM + LRT testing  
   - 27 significant DMRs: 11 hypo-, 16 hyper-methylated  
   - Top regions: `GGA1`, `TOM1L1`, `RPL5P10`, `NT5DC1P1`

3. **Functional Enrichment (R / Bioconductor)**  
   - Tools: `clusterProfiler`, `ReactomePA`, `GOstats`, `g:Profiler`  
   - Significant Reactome pathways: collagen formation, IGF signaling  
   - Key gene hits: `IGFBP3`, `COL1A1`, `GGA1`

4. **Machine Learning (Python / scikit-learn)**  
   - **Model 1:** Two logistic regressions (RNA + Methylation) → Meta ensemble  
     - Combined probability averaging  
     - Best performance: Accuracy = 0.92, AUC = 0.94  
     - Key predictive genes: `KLF2P1`, `SCAF1`, `IGF2`, `CDCA2`
   - **Model 2:** Stacked classifier (Logistic + Decision Tree + SVM)  
     - Combined DEGs + DMRs  
     - Accuracy = 0.83, AUC = 0.83  



---
## 🧰 Tech Stack


| Category | Tools & Libraries |
|-----------|------------------|
| Languages | **Python 3.11**, **R 4.3** |
| RNA-Seq | DESeq2, biomaRt, org.Hs.eg.db |
| Methylation | edgeR, Bioconductor |
| Enrichment | clusterProfiler, ReactomePA, GOstats, g:Profiler |
| Machine Learning | scikit-learn (Logistic Regression, Decision Tree, SVM, StackingClassifier) |
| Visualization | matplotlib, seaborn, ggplot2 |
| Utilities | pandas, numpy, MinMaxScaler, GridSearchCV |


---

## 📁 Project Structure




--- 
**Additionally, a readme and a report with a detailed description of the project included.**
