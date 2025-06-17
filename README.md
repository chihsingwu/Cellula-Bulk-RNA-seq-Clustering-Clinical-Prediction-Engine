# Cellula-Bulk-RNA-seq-Clustering-Clinical-Prediction-Engine
This project provides a modular engine for rapid clustering and dimensionality reduction (e.g., PCA + KMeans) of bulk RNA-seq data, with integrated clinical outcome prediction (Cox model, Kaplan-Meier curves).  The default demo uses the S100 protein family (e.g., S100A4/A8/A9/A11).
Overview
This project provides a modular engine for rapid clustering and dimensionality reduction (e.g., PCA + KMeans) of bulk RNA-seq data, with integrated clinical outcome prediction (Cox model, Kaplan-Meier curves).
The default demo uses the S100 protein family (e.g., S100A4/A8/A9/A11), which are key biomarkers in cancer progression and tumor microenvironment studies. Users can easily specify any set of genes for analysis.
Input: bulk RNA-seq matrix (patients × genes) and clinical data (time, event, etc.)
Output: cluster assignments, 2D plots, and survival analysis by cluster
Features
Core Analysis Engines
1. ClusterEngine - Dimensionality Reduction & Clustering
PCA dimensionality reduction for high-dimensional gene expression data
K-means clustering to identify patient subgroups
Interactive visualization with explained variance ratios
Flexible gene selection (default: S100 protein family)
pythoncluster_engine = ClusterEngine(n_components=2, n_clusters=2)
X_reduced, cluster_labels = cluster_engine.fit_transform(expression_matrix)
cluster_engine.plot_clusters(X_reduced, cluster_labels)
2. ClinicalEngine - Survival Analysis & Prognosis
Cox proportional hazards modeling for risk assessment
Kaplan-Meier survival curves with statistical testing
Clinical outcome prediction based on molecular clusters
Multi-group comparison with automatic statistical testing

pythonclinical_engine = ClinicalEngine()
cox_model = clinical_engine.fit_predict(clinical_data, cluster_labels)
clinical_engine.plot_survival(clinical_data, cluster_labels)
3. ValidationEngine - Statistical Validation & Quality Control
Statistical significance testing (Log-rank test for survival differences)
Biomarker performance evaluation (ROC-AUC analysis)
Clustering stability assessment (Adjusted Rand Index across multiple runs)
Comprehensive validation reporting with confidence assessment
python# Integrated validation pipeline
validation_results = cluster_engine.validation_summary(
    X_reduced, expression_data, clinical_data, cluster_labels
)

# Example output:
{
    'statistical_validation': {'p_value': 0.023, 'significant': True},
    'biomarker_performance': {'S100A4': {'auc': 0.82, 'diagnostic_value': 'Good'}},
    'cluster_stability': {'mean_stability': 0.75, 'stable': True},
    'overall_assessment': {'recommendation': 'High confidence results'}
}
Installation
bashpip install numpy pandas scikit-learn lifelines matplotlib
Quick Start
pythonfrom genomics_cluster_engine import ClusterEngine, ClinicalEngine

# Initialize engines
cluster_engine = ClusterEngine(n_components=2, n_clusters=2)
clinical_engine = ClinicalEngine()

# Load your data
# X: Gene expression matrix (samples × genes)
# clinical_df: DataFrame with ['time', 'event'] columns

# 1. Clustering analysis
X_reduced, cluster_labels = cluster_engine.fit_transform(X)
cluster_engine.plot_clusters(X_reduced, cluster_labels)

# 2. Survival analysis
cox_model = clinical_engine.fit_predict(clinical_df, cluster_labels)
clinical_engine.plot_survival(clinical_df, cluster_labels)

# 3. Validation & quality assessment
validation = cluster_engine.validation_summary(
    X_reduced, X_df, clinical_df, cluster_labels
)
print(f"Overall assessment: {validation['overall_assessment']['recommendation']}")
Data Format Requirements
Gene Expression Data (X)
      S100A4  S100A8  S100A9  S100A11  ...
Sample_001   5.2     3.1     4.8      2.9   ...
Sample_002   7.1     2.3     6.2      4.1   ...
Sample_003   3.8     5.5     2.1      3.6   ...
Clinical Data (clinical_df)
      time    event   age   stage  ...
Sample_001  24.5    1      65    II   ...
Sample_002  18.2    0      58    III  ...
Sample_003  36.1    1      72    I    ...

time: Follow-up time (months/days)
event: Event indicator (1=death/progression, 0=censored)

S100 Protein Family Default Analysis
The engine is pre-configured for S100 protein family analysis, focusing on:

S100A4: Metastasis-associated protein
S100A8: Neutrophil cytosolic factor
S100A9: Calgranulin B
S100A11: Calgizzarin

These proteins are crucial in:
Cancer progression and metastasis
Inflammatory responses
Tumor microenvironment modulation
Drug resistance mechanisms

Validation Framework
Statistical Validation
Log-rank test: Assesses survival differences between clusters
P-value interpretation: Statistical significance threshold (p < 0.05)

Biomarker Performance
ROC-AUC analysis: Evaluates diagnostic potential of individual genes
Performance categories: Excellent (>0.9), Good (>0.8), Fair (>0.7), Poor (≤0.7)

Clustering Stability
Adjusted Rand Index: Measures consistency across multiple K-means runs
Stability threshold: Stable clustering (ARI > 0.7)

Advanced Usage
Custom Gene Sets
python# Analyze custom gene panel
custom_genes = ['BRCA1', 'BRCA2', 'TP53', 'MYC']
validation = cluster_engine.biomarker_performance(
    expression_data, cluster_labels, target_genes=custom_genes
)
Multi-cluster Analysis
python# Analyze multiple patient subgroups
cluster_engine = ClusterEngine(n_components=3, n_clusters=4)
X_reduced, labels = cluster_engine.fit_transform(X)
Batch Analysis
python# Process multiple datasets
for dataset_name, data in datasets.items():
    X_reduced, labels = cluster_engine.fit_transform(data['expression'])
    validation = cluster_engine.validation_summary(
        X_reduced, data['expression'], data['clinical'], labels
    )
    print(f"{dataset_name}: {validation['overall_assessment']['recommendation']}")
Output Interpretation
High Confidence Results

Statistically significant survival differences (p < 0.05)
Stable clustering (ARI > 0.7)
Good biomarker performance (AUC > 0.7)

Moderate Confidence Results
Either statistical significance OR stable clustering
Some biomarkers show diagnostic potential
Low Confidence Results
No statistical significance
Unstable clustering
Poor biomarker performance
Recommendation: Increase sample size or try different parameters

Applications
Cancer subtype identification
Prognostic biomarker discovery
Drug response prediction
Clinical trial stratification
Personalized medicine research

Contributing
Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.
License
This project is licensed under the MIT License.
Citation
If you use this engine in your research, please cite:
Enhanced Clinical Genomics Clustering Engine
GitHub: [Your Repository URL]

Note: This engine is designed for research purposes. Clinical applications should undergo appropriate validation and regulatory approval.
