# ==============================================================================
#  A Self-Contained Script for Simulation, Analysis, and Visualization
# ==============================================================================

import numpy as np
import pandas as pd
import umap
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# ==============================================================================
#  1. Data Simulation Function
# ==============================================================================
def simulate_data(n_samples=300, n_features=1000, n_clusters=3):
    """
    Simulates gene expression and corresponding clinical data with distinct clusters
    and different survival outcomes.
    """
    print(f"--> Generating expression data with {n_samples} samples, {n_features} features, and {n_clusters} clusters...")
    
    # Generate clustered data points
    X, y_true = make_blobs(n_samples=n_samples, centers=n_clusters, 
                           n_features=n_features, random_state=42, cluster_std=2.0)
    
    patient_ids = [f'Patient_{i+1:03d}' for i in range(n_samples)]
    expression_df = pd.DataFrame(X, index=patient_ids, columns=[f'Gene_{i+1}' for i in range(n_features)])
    
    print("--> Generating clinical data with different survival profiles...")
    times = []
    events = []
    
    # Define survival characteristics for each true cluster
    # Cluster 0: Poor prognosis, Cluster 1: Good prognosis, Cluster 2: Intermediate
    survival_params = {0: 2.0, 1: 10.0, 2: 5.0} 
    
    for label in y_true:
        # Simulate survival time from a Weibull distribution
        time = np.random.weibull(a=survival_params[label]) * 20 + 5
        times.append(time)
        
        # Simulate censoring (approx. 20% of data will be censored)
        events.append(np.random.choice([0, 1], p=[0.2, 0.8]))
        
    clinical_df = pd.DataFrame({'time': times, 'event': events}, index=patient_ids)
    
    return expression_df, clinical_df

# ==============================================================================
#  2. Analysis and Plotting Functions
# ==============================================================================
def run_analysis_and_plot(expression_data, clinical_data, n_clusters):
    """
    Runs the entire analysis pipeline and generates plots.
    """
    # --- UMAP Reduction and K-means Clustering ---
    print("\n--> Performing UMAP dimensionality reduction...")
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=30,  # A key parameter for UMAP
        min_dist=0.0,    # Pack points tightly
        random_state=42
    )
    X_reduced = reducer.fit_transform(expression_data)

    print("--> Performing K-means clustering on UMAP results...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init='auto')
    predicted_labels = kmeans.fit_predict(X_reduced)
    
    # --- Plotting UMAP Clusters ---
    print("--> Generating UMAP cluster plot...")
    plt.figure(figsize=(10, 8))
    unique_labels = sorted(np.unique(predicted_labels))
    colors = plt.cm.viridis(np.linspace(0, 1, len(unique_labels)))
    
    for i, color in zip(unique_labels, colors):
        mask = predicted_labels == i
        plt.scatter(X_reduced[mask, 0], X_reduced[mask, 1], color=color, 
                    label=f'Predicted Cluster {i}', alpha=0.8, s=60, 
                    edgecolors='k', linewidth=0.5)
    
    plt.xlabel("UMAP Dimension 1", fontsize=12)
    plt.ylabel("UMAP Dimension 2", fontsize=12)
    plt.title("UMAP + K-means Clustering of Simulated Data", fontsize=16)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig("umap_cluster_plot.png")
    plt.close()
    print("    ... saved to umap_cluster_plot.png")

    # --- Plotting Survival Curves ---
    print("--> Generating Kaplan-Meier survival plot...")
    plt.figure(figsize=(10, 8))
    data_for_survival = clinical_data.copy()
    data_for_survival['cluster'] = predicted_labels
    
    for cluster_id, color in zip(unique_labels, colors):
        mask = data_for_survival['cluster'] == cluster_id
        kmf = KaplanMeierFitter()
        kmf.fit(data_for_survival.loc[mask, 'time'], 
                event_observed=data_for_survival.loc[mask, 'event'], 
                label=f'Cluster {cluster_id} (n={sum(mask)})')
        kmf.plot_survival_function(color=color, linewidth=3)
    
    results = logrank_test(data_for_survival['time'], data_for_survival['cluster'], data_for_survival['event'])
    p_val = results.p_value
    
    plt.title(f"Kaplan-Meier Survival Curves by Predicted Cluster\n(Overall Log-rank test p-value: {p_val:.2e})", fontsize=16)
    plt.xlabel("Time (Simulated Months)", fontsize=12)
    plt.ylabel("Survival Probability", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.ylim(0, 1.05)
    plt.tight_layout()
    plt.savefig("survival_plot.png")
    plt.close()
    print("    ... saved to survival_plot.png")

# ==============================================================================
#  3. Main Execution
# ==============================================================================

print("🚀 Starting Full Simulation and Visualization Script...")

# Simulate the data
expression_df, clinical_df = simulate_data(n_samples=300, n_features=1000, n_clusters=3)
print("✅ Data simulation complete.")

# Run the entire analysis and plotting pipeline
run_analysis_and_plot(expression_df, clinical_df, n_clusters=3)

print("\n🎉 Pipeline finished successfully!")
print("Please check your folder for 'umap_cluster_plot.png' and 'survival_plot.png'.")