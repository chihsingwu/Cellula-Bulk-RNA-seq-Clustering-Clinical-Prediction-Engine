# Enhanced Clinical Genomics Clustering Engine
# A streamlined RNA-seq analysis pipeline with integrated statistical validation
# for biomarker discovery and survival analysis

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import roc_auc_score, adjusted_rand_score
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

class ClusterEngine:
    """
    Advanced clustering engine for RNA-seq data analysis with PCA dimensionality reduction
    and K-means clustering, specifically designed for biomarker discovery.
    """
    
    def __init__(self, n_components=2, n_clusters=2):
        """
        Initialize the clustering engine.
        
        Parameters:
        -----------
        n_components : int, default=2
            Number of principal components for PCA reduction
        n_clusters : int, default=2
            Number of clusters for K-means algorithm
        """
        self.n_components = n_components
        self.n_clusters = n_clusters
        self.pca = None
        self.kmeans = None
        self.cluster_labels_ = None

    def fit_transform(self, X):
        """
        Perform PCA dimensionality reduction followed by K-means clustering.
        
        Parameters:
        -----------
        X : array-like, shape (n_samples, n_features)
            Gene expression matrix (samples × genes)
            
        Returns:
        --------
        X_reduced : array-like, shape (n_samples, n_components)
            PCA-transformed data
        cluster_labels : array-like, shape (n_samples,)
            Cluster assignments for each sample
        """
        # PCA dimensionality reduction
        self.pca = PCA(n_components=self.n_components)
        X_reduced = self.pca.fit_transform(X)
        
        # K-means clustering
        self.kmeans = KMeans(n_clusters=self.n_clusters, random_state=42)
        cluster_labels = self.kmeans.fit_predict(X_reduced)
        self.cluster_labels_ = cluster_labels
        
        return X_reduced, cluster_labels

    def plot_clusters(self, X_reduced, labels):
        """
        Visualize clustering results in PCA space.
        
        Parameters:
        -----------
        X_reduced : array-like, shape (n_samples, n_components)
            PCA-transformed data
        labels : array-like, shape (n_samples,)
            Cluster labels
        """
        plt.figure(figsize=(8, 6))
        colors = plt.cm.Set1(np.linspace(0, 1, self.n_clusters))
        
        for i, color in zip(np.unique(labels), colors):
            mask = labels == i
            plt.scatter(X_reduced[mask, 0], X_reduced[mask, 1], 
                       c=[color], label=f'Cluster {i}', alpha=0.7, s=50)
        
        plt.xlabel(f"PC1 ({self.pca.explained_variance_ratio_[0]:.1%} variance)")
        plt.ylabel(f"PC2 ({self.pca.explained_variance_ratio_[1]:.1%} variance)")
        plt.title("PCA + K-means Clustering Results")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    # === Validation Functions ===
    
    def survival_significance(self, clinical_data, cluster_labels):
        """
        Perform log-rank test to assess survival differences between clusters.
        
        Parameters:
        -----------
        clinical_data : DataFrame
            Must contain 'time' and 'event' columns
        cluster_labels : array-like
            Cluster assignments
            
        Returns:
        --------
        dict : Statistical test results
        """
        data = clinical_data.copy()
        data['cluster'] = cluster_labels
        
        groups = []
        for cluster_id in np.unique(cluster_labels):
            cluster_data = data[data['cluster'] == cluster_id]
            groups.append((cluster_data['time'], cluster_data['event']))
        
        if len(groups) >= 2:
            result = logrank_test(*groups[0], *groups[1])
            return {
                'p_value': result.p_value,
                'test_statistic': result.test_statistic,
                'significant': result.p_value < 0.05,
                'interpretation': 'Significant survival difference' if result.p_value < 0.05 else 'No significant difference'
            }
        else:
            return {'error': 'Need at least 2 clusters for comparison'}

    def biomarker_performance(self, expression_data, cluster_labels, target_genes=None):
        """
        Evaluate diagnostic performance of target genes using ROC-AUC analysis.
        
        Parameters:
        -----------
        expression_data : DataFrame
            Gene expression data with genes as columns
        cluster_labels : array-like
            Binary cluster labels (0/1)
        target_genes : list, optional
            Specific genes to evaluate. Defaults to S100 family genes.
            
        Returns:
        --------
        dict : AUC scores and diagnostic potential for each gene
        """
        if target_genes is None:
            target_genes = ['S100A4', 'S100A8', 'S100A9', 'S100A11']
        
        results = {}
        available_genes = [gene for gene in target_genes if gene in expression_data.columns]
        
        for gene in available_genes:
            try:
                auc = roc_auc_score(cluster_labels, expression_data[gene])
                results[gene] = {
                    'auc': round(auc, 3),
                    'diagnostic_value': 'Excellent' if auc > 0.9 else 
                                      'Good' if auc > 0.8 else 
                                      'Fair' if auc > 0.7 else 'Poor'
                }
            except Exception as e:
                results[gene] = {'error': str(e)}
        
        return results

    def cluster_stability(self, X, n_runs=10):
        """
        Assess clustering stability through repeated K-means runs.
        
        Parameters:
        -----------
        X : array-like
            PCA-transformed data
        n_runs : int, default=10
            Number of stability test runs
            
        Returns:
        --------
        dict : Stability metrics
        """
        if self.cluster_labels_ is None:
            return {'error': 'Must run fit_transform first'}
        
        stability_scores = []
        base_labels = self.cluster_labels_
        
        for _ in range(n_runs):
            kmeans_temp = KMeans(n_clusters=self.n_clusters, random_state=None)
            temp_labels = kmeans_temp.fit_predict(X)
            stability_scores.append(adjusted_rand_score(base_labels, temp_labels))
        
        mean_stability = np.mean(stability_scores)
        return {
            'mean_stability': round(mean_stability, 3),
            'std_stability': round(np.std(stability_scores), 3),
            'stable': mean_stability > 0.7,
            'interpretation': 'Highly stable' if mean_stability > 0.8 else
                           'Moderately stable' if mean_stability > 0.6 else 'Unstable'
        }

    def validation_summary(self, X, expression_data, clinical_data, cluster_labels):
        """
        Generate comprehensive validation report for clustering results.
        
        Parameters:
        -----------
        X : array-like
            PCA-transformed data
        expression_data : DataFrame
            Original gene expression data
        clinical_data : DataFrame
            Clinical data with survival information
        cluster_labels : array-like
            Cluster assignments
            
        Returns:
        --------
        dict : Complete validation results
        """
        validation_results = {
            'statistical_validation': self.survival_significance(clinical_data, cluster_labels),
            'biomarker_performance': self.biomarker_performance(expression_data, cluster_labels),
            'cluster_stability': self.cluster_stability(X)
        }
        
        # Overall assessment
        stat_sig = validation_results['statistical_validation'].get('significant', False)
        stability = validation_results['cluster_stability'].get('stable', False)
        good_biomarkers = sum(1 for gene_result in validation_results['biomarker_performance'].values() 
                             if isinstance(gene_result, dict) and gene_result.get('auc', 0) > 0.7)
        
        validation_results['overall_assessment'] = {
            'statistically_significant': stat_sig,
            'stable_clustering': stability,
            'good_biomarkers_found': good_biomarkers,
            'recommendation': 'High confidence results' if (stat_sig and stability and good_biomarkers > 0) else
                            'Moderate confidence results' if (stat_sig or stability) else
                            'Results require further validation'
        }
        
        return validation_results


class ClinicalEngine:
    """
    Clinical survival analysis engine for prognostic modeling and visualization.
    """
    
    def __init__(self):
        """Initialize clinical analysis engine."""
        self.model = None

    def fit_predict(self, clinical_data, cluster_labels):
        """
        Fit Cox proportional hazards model for survival analysis.
        
        Parameters:
        -----------
        clinical_data : DataFrame
            Must contain 'time' and 'event' columns
        cluster_labels : array-like
            Cluster assignments as predictor variable
            
        Returns:
        --------
        CoxPHFitter : Fitted Cox regression model
        """
        data = clinical_data.copy()
        data['cluster'] = cluster_labels
        
        self.model = CoxPHFitter()
        self.model.fit(data, duration_col='time', event_col='event', formula="cluster")
        
        return self.model

    def plot_survival(self, clinical_data, cluster_labels):
        """
        Generate Kaplan-Meier survival curves for each cluster.
        
        Parameters:
        -----------
        clinical_data : DataFrame
            Clinical data with survival information
        cluster_labels : array-like
            Cluster assignments
        """
        plt.figure(figsize=(10, 6))
        
        data = clinical_data.copy()
        data['cluster'] = cluster_labels
        
        colors = plt.cm.Set1(np.linspace(0, 1, len(np.unique(cluster_labels))))
        
        for cluster_id, color in zip(np.unique(cluster_labels), colors):
            mask = data['cluster'] == cluster_id
            cluster_data = data[mask]
            
            kmf = KaplanMeierFitter()
            kmf.fit(cluster_data['time'], 
                   event_observed=cluster_data['event'], 
                   label=f'Cluster {cluster_id} (n={sum(mask)})')
            
            kmf.plot_survival_function(color=color, linewidth=2)
        
        plt.title("Kaplan-Meier Survival Curves by Cluster")
        plt.xlabel("Time")
        plt.ylabel("Survival Probability")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()


# Main analysis pipeline
if __name__ == "__main__":
    """
    Example usage of the Enhanced Clinical Genomics Clustering Engine
    
    Required data format:
    - X: Gene expression matrix (samples × genes), e.g., S100A4/A8/A9/A11 expression
    - clinical_df: DataFrame with columns ['time', 'event'] + other clinical variables
    """
    
    # Initialize engines
    cluster_engine = ClusterEngine(n_components=2, n_clusters=2)
    clinical_engine = ClinicalEngine()
    
    # Example workflow (uncomment when you have real data)
    """
    # 1. Clustering and dimensionality reduction
    X_reduced, cluster_labels = cluster_engine.fit_transform(X)
    cluster_engine.plot_clusters(X_reduced, cluster_labels)
    
    # 2. Clinical survival analysis
    cox_model = clinical_engine.fit_predict(clinical_df, cluster_labels)
    clinical_engine.plot_survival(clinical_df, cluster_labels)
    
    # 3. Comprehensive validation
    validation_results = cluster_engine.validation_summary(
        X_reduced, X_df, clinical_df, cluster_labels)
    
    # Print results
    print("=== Cox Regression Results ===")
    print(cox_model.summary)
    
    print("\n=== Validation Summary ===")
    for category, results in validation_results.items():
        print(f"{category}: {results}")
    """
    
    print("Enhanced Clinical Genomics Clustering Engine initialized successfully!")
    print("Ready for RNA-seq analysis with integrated statistical validation.")