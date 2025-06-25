# Cellula-Bulk-RNA-seq-Clustering-Clinical-Prediction-Engine
This project provides a modular engine for rapid clustering and dimensionality reduction (e.g., PCA + KMeans) of bulk RNA-seq data, with integrated clinical outcome prediction (Cox model, Kaplan-Meier curves).  The default demo uses the S100 protein family (e.g., S100A4/A8/A9/A11).
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

> A streamlined RNA-seq analysis pipeline with integrated survival modeling and statistical validation.

## 🎯 Overview

This project provides a **modular engine** for rapid clustering and dimensionality reduction (PCA + K-means) of bulk RNA-seq data, with integrated clinical outcome prediction (Cox model, Kaplan-Meier curves).

**Default focus:** S100 protein family (S100A4/A8/A9/A11) - key biomarkers in cancer progression  
**Flexible usage:** Easily specify any gene set for analysis

| Input | Output |
|-------|--------|
| 📊 RNA-seq matrix (patients × genes) | 🎯 Patient clusters |
| 🏥 Clinical data (time, event) | 📈 PCA visualization |
| | 📉 Survival curves & statistics |

## 🚀 Quick Start

```bash
pip install numpy pandas scikit-learn lifelines matplotlib
```
pip install numpy pandas scikit-learn matplotlib lifelines umap-learn **vo3 version*

```python
from genomics_cluster_engine import ClusterEngine, ClinicalEngine

# Initialize
cluster_engine = ClusterEngine(n_components=2, n_clusters=2)
clinical_engine = ClinicalEngine()

# Analysis pipeline
X_reduced, labels = cluster_engine.fit_transform(expression_data)
cluster_engine.plot_clusters(X_reduced, labels)

cox_model = clinical_engine.fit_predict(clinical_data, labels)
clinical_engine.plot_survival(clinical_data, labels)

# Validation
validation = cluster_engine.validation_summary(X_reduced, expression_data, clinical_data, labels)
print(f"Result confidence: {validation['overall_assessment']['recommendation']}")
```

## 🧬 Core Features

### **ClusterEngine** - Dimensionality Reduction & Clustering
- PCA dimensionality reduction for high-dimensional gene expression
- K-means clustering to identify patient subgroups  
- Integrated statistical validation (log-rank test, ROC-AUC, stability)

### **ClinicalEngine** - Survival Analysis
- Cox proportional hazards modeling
- Kaplan-Meier survival curves with statistical testing
- Multi-group comparison and visualization

### **Validation Framework**
- **Statistical significance**: Log-rank test for survival differences
- **Biomarker performance**: ROC-AUC analysis for diagnostic potential
- **Clustering stability**: Adjusted Rand Index across multiple runs

## 📊 Data Format

### Expression Data
```
        S100A4  S100A8  S100A9  S100A11
Patient_001  5.2    3.1     4.8     2.9
Patient_002  7.1    2.3     6.2     4.1
```

### Clinical Data  
```
        time  event  age  stage
Patient_001  24.5    1    65    II
Patient_002  18.2    0    58   III
```
- `time`: Follow-up time (months/days)
- `event`: Event indicator (1=event, 0=censored)

## 🔬 S100 Protein Family

Pre-configured analysis for key cancer biomarkers:
- **S100A4**: Metastasis-associated protein
- **S100A8**: Neutrophil cytosolic factor  
- **S100A9**: Calgranulin B
- **S100A11**: Calgizzarin

Critical roles in cancer progression, inflammation, and tumor microenvironment.

## 📈 Validation Output Example

```python
{
    'statistical_validation': {'p_value': 0.023, 'significant': True},
    'biomarker_performance': {'S100A4': {'auc': 0.82, 'diagnostic_value': 'Good'}},
    'cluster_stability': {'mean_stability': 0.75, 'stable': True},
    'overall_assessment': {'recommendation': 'High confidence results'}
}
```

## 🎨 Advanced Usage

```python
# Custom gene panel
custom_genes = ['BRCA1', 'BRCA2', 'TP53']
validation = cluster_engine.biomarker_performance(data, labels, target_genes=custom_genes)

# Multi-cluster analysis  
cluster_engine = ClusterEngine(n_components=3, n_clusters=4)

# Batch processing
for dataset in datasets:
    results = cluster_engine.validation_summary(...)
```

## 🏥 Applications

- **Cancer subtype identification**
- **Prognostic biomarker discovery** 
- **Clinical trial stratification**
- **Drug response prediction**
- **Personalized medicine research**

## 📚 Documentation

- [Installation Guide](docs/installation.md)
- [API Reference](docs/api.md)
- [Tutorial Notebooks](examples/)
- [Contributing Guidelines](CONTRIBUTING.md)

## 🤝 Contributing

Contributions welcome! Please read our [contributing guidelines](CONTRIBUTING.md) and submit pull requests.

## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## 📎 Citation

```bibtex
@software{clinical_genomics_engine,
  title={Clinical Genomics Clustering Engine},
  author={[Your Name]},
  year={2025},
  url={https://github.com/[username]/[repo-name]},
  doi={10.5281/zenodo.XXXXXXX}
}
```

---

**⚠️ Note**: This tool is designed for research purposes. Clinical applications require appropriate validation and regulatory approval.
