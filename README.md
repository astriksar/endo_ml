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

## 📁 Repository Content

Diff_Analysis_mbdseq.R        → MBD-seq differential methylation analysis  
GoTermEnrichment.Rmd          → GO enrichment pipeline  
GoTermEnrichment.html         → Rendered enrichment report  
ML_model1.ipynb / .html       → Machine learning model 1  
ML_model2.ipynb / .html       → Machine learning model 2  
MethylationAnnotation.R       → Methylation site annotation  
rnaseq.Rmd / .html            → RNA-seq analysis pipeline  
Report_Group3.pdf             → Final report  
ReadMe.docx.pdf               → Documentation  
README.md                     → Project overview (this file)


## 📈 Model Performance & Results

The objective of this project was to build a **high-performance machine learning classifier** capable of predicting endometriosis status using integrated **multi-omics data** (RNA-Seq + Methylation).  
The results clearly show that **ensemble learning** and **meta-model integration** significantly outperform baseline models trained on individual data types.

---

### 🔹 1. Ensemble Learning vs. Baseline Models

A **Stacked Classifier** was developed by combining Logistic Regression, Decision Tree, and SVM base learners.  
It achieved the **highest and most balanced performance** across all metrics.

| **Classifier** | **Test Accuracy** | **F1-score** | **AUC** |
| :-------------- | :---------------: | :-----------: | :-----: |
| Logistic Regression | 0.75 | 0.77 | 0.75 |
| Decision Tree | 0.42 | 0.22 | 0.42 |
| Support Vector Machine | 0.75 | 0.73 | 0.75 |
| **Stacked Classifier** | **0.83** | **0.83** | **0.83** |

**Conclusion:**  
The ensemble approach boosts predictive performance and improves stability, increasing the F1-score and AUC to **0.83**.

---

### 🔹 2. Multi-Omics Meta-Model Integration

To further enhance prediction, a **Meta-Model** was trained by combining features derived from both **RNA expression** and **DNA methylation** datasets.  
This model **substantially outperformed** any single-omics model.

| **Classifier** | **Test Accuracy** | **F1-score** | **AUC** |
| :-------------- | :---------------: | :-----------: | :-----: |
| RNA Logistic Regression | 0.62 | 0.67 | 0.69 |
| Methylation Logistic Regression | 0.77 | 0.67 | 0.74 |
| **Meta-Model (RNA + Methylation)** | **0.92** | **0.91** | **0.94** |

**Conclusion:**  
Integrating multi-omics features delivers **state-of-the-art performance**, achieving an **F1-score of 0.91** and **AUC of 0.94** — demonstrating the effectiveness of **data fusion** in endometriosis biomarker prediction.

---

📊 *Overall Insight:*  
Multi-omics integration and ensemble learning dramatically improve classification accuracy and robustness, validating this approach for complex biological datasets.


**Highlighted features:**  
- *IGF2* (RNA downregulated, literature-supported biomarker)  
- *CDCA2* (hypermethylated region)  
- *IGFBP3* (methylation linked to lesion development)  
- *KLF2P1*, *SCAF1* (novel candidates)

---

## 🔬 Limitations & Future Work

- Small overlapping sample size between RNA and MBD datasets  
- Limited external validation due to data availability  
- Future directions:
  - Larger cohort validation  
  - Bootstrapping for variance reduction  
  - Deep learning integration (autoencoders / multimodal fusion)  
  - Expanded omics integration (proteomics, metabolomics)  
  - Web API or Streamlit interface for inference  

---

## ⚙️ Usage

```bash
# Clone repository
git clone https://github.com/astriksar/endo_ml.git
cd endo_ml
```

--- 
