# Alzheimers_GSE5281_Analysis
Comprehensive gene expression and pathway analysis of the GSE5281 Alzheimer’s Disease dataset using R (limma, clusterProfiler, and machine learning models).
# 🧠 Gene Expression Data Analysis from Public Microarray Datasets for Alzheimer’s Disease

**Author:** Ishaan Bhardwaj  
**Guide:** Dr. Rajendra Kumar  
**Institution:** Department of Bioinformatics, Jamia Millia Islamia  
**Year:** 2025  
**Dataset:** [GSE5281 - NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281)

---

## 📘 Overview
This repository contains all scripts, results, and outputs from a gene expression analysis of Alzheimer’s Disease using the public microarray dataset **GSE5281**.  
The workflow integrates **differential expression**, **pathway enrichment**, and **machine learning classification** to uncover key molecular mechanisms of AD.

---

## ⚙️ Workflow Summary
1. **Data retrieval:** GEOquery (GSE5281)
2. **Preprocessing:** Log₂ transformation and probe filtering  
3. **DEG identification:** limma package (AD vs Control)
4. **Enrichment analysis:** GO (BP) and KEGG using clusterProfiler  
5. **Visualization:** Volcano, PCA, Heatmap
6. **Machine Learning:** Random Forest & SVM (caret, e1071)

---

## 📊 Key Findings
- **302 strong DEGs** (adj. p < 0.05, |log₂FC| ≥ 1)  
- **6,148 significant probes** before filtering  
- **Top DEGs:** APP, PSEN1, MAPT, SORL1, APOE  
- **Top enriched GO terms:** Synaptic signaling, Neuron projection morphogenesis  
- **Key KEGG pathways:** Alzheimer’s disease, MAPK signaling, Oxidative phosphorylation  
- **ML Accuracy:** Random Forest (≈85%), SVM (≈82%)
📁 Repository Structure
├── GSE5281_pipeline.R # Full analysis script
├── enrichment_analysis.R # GO/KEGG enrichment
├── output/
│ ├── GSE5281_EC_limma_all.csv
│ ├── GSE5281_EC_GO_BP.csv
│ ├── GSE5281_EC_KEGG.csv
│ ├── GSE5281_EC_volcano.png
│ ├── GSE5281_EC_PCA.png
│ └── GSE5281_EC_top50_heatmap.png
├── saved_objects/
│ ├── rf_model.rds
│ ├── svm_model.rds
│ └── ml_data.rds
└── ppt/
└── Ishaan_Bhardwaj_Alzheimer_GSE5281_Presentation.pptx

---

## 🧩 Tools and Packages
- R (v4.5.1)
- GEOquery, Biobase, limma
- hgu133plus2.db, org.Hs.eg.db
- clusterProfiler, EnhancedVolcano, pheatmap
- caret, randomForest, e1071, ggplot2

---

## 🧬 Biological Insight
The results highlight dysregulation in **synaptic**, **oxidative phosphorylation**, and **mitochondrial** pathways — core mechanisms implicated in Alzheimer’s progression.

---

## 📜 License
Released under the **MIT License** – you are free to reuse and modify with attribution.

---

## 🎓 Citation
> Ishaan Bhardwaj (2025). *Gene Expression Data Analysis from Public Microarray Datasets for Alzheimer’s Disease.*  
> MSc Bioinformatics, Jamia Millia Islamia. 

---

## 💬 Contact
📧 ishaanbhardwaj007@gmail.com
---

## 📁 Repository Structure
