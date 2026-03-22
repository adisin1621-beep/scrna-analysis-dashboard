# 🔬 Single-Cell RNA-seq Analysis Dashboard

An interactive web-based platform for **single-cell RNA sequencing (scRNA-seq) analysis**, built using **Scanpy, Streamlit, and Python**.

---

## 🔗 Live Demo

https://scrna-analysis-dashboard-rnjqqts7f2jbtgyncnidpp.streamlit.app/

---

## 🚀 Features

### 📊 Data Handling

* Load demo datasets (PBMC)
* Upload custom `.h5ad` datasets
* Load datasets directly via URL

### 🔬 Analysis Pipeline

* Normalization & log transformation
* Highly variable gene selection
* PCA & neighborhood graph
* Leiden clustering
* UMAP visualization

### 🧬 Marker Gene Analysis

* Marker gene identification
* Heatmap visualization
* Dotplot visualization
* Export marker gene tables

### 📈 Differential Expression

* Cluster vs cluster comparison
* Volcano plot with:

  * log2 fold change
  * p-value thresholds
  * gene labeling
* Export DE results
* Download significant genes

### 🧪 Cell Type Annotation

* Marker-based annotation assistant

### 🧭 Trajectory Analysis

* Diffusion pseudotime (DPT)

---

## 🖥️ Run Locally

```bash
pip install -r requirements.txt
streamlit run Dashboard/app.py
```

---

## 📁 Project Structure

```
Scrna-immune-cell-analysis
│
├── Dashboard/
│   └── app.py
├── dataset_loader.py
├── analysis_pipeline.py
├── requirements.txt
├── README.md
```

---

## 🧠 Tech Stack

* Python
* Scanpy
* Streamlit
* Plotly
* Matplotlib

---

## 👨‍💻 Author

Adi
Bioinformatics / AI-ML Enthusiast

---

## ⭐ Future Improvements

* GEO dataset integration
* Automatic cell type annotation
* Pathway enrichment visualization
* UI improvements

---
