```python
import numpy as np
import pandas as pd
from sklearn.datasets import make_blobs
import warnings
warnings.filterwarnings('ignore')

def simulate_data(n_samples=300, n_features=1000, n_clusters=3):
    """
    Simulates a gene expression dataset with distinct clusters and
    corresponding clinical data where clusters have different survival outcomes.
    """
    # 1. Simulate Expression Data with clear clusters
    print(f"Generating expression data with {n_samples} samples, {n_features} features, and {n_clusters} clusters...")
    X, y_true = make_blobs(n_samples=n_samples, centers=n_clusters, 
                           n_features=n_features, random_state=42, cluster_std=2.0)
    
    gene_names = [f'Gene_{i+1}' for i in range(n_features)]
    patient_ids = [f'Patient_{i+1:03d}' for i in range(n_samples)]
    X_df = pd.DataFrame(X, index=patient_ids, columns=gene_names)

    # 2. Simulate Clinical Data linked to the true clusters
    print("Generating corresponding clinical data with different survival profiles...")
    times = []
    events = []
    
    # Cluster 0: Poor prognosis, Cluster 1: Good prognosis, Cluster 2: Intermediate
    survival_params = {0: 2.0, 1: 10.0, 2: 5.0} 
    
    for label in y_true:
        # Simulate survival time from Weibull distribution, add baseline time
        time = np.random.weibull(a=survival_params[label]) * 20 + 5
        times.append(time)
        
        # Simulate censoring event (80% non-censored)
        event = np.random.choice([0, 1], p=[0.2, 0.8])
        events.append(event)
        
    clinical_df = pd.DataFrame({'time': times, 'event': events}, index=patient_ids)
    
    # Add true cluster labels for later validation if needed
    clinical_df['true_cluster'] = y_true

    print("Simulation complete.")
    return X_df, clinical_df

# ==============================================================================
#  Main Execution Script to Generate and Save Data
# ==============================================================================
print("🚀 STEP 1: Starting Data Simulation and File Generation...")

# --- Simulation ---
N_CLUSTERS = 3
expression_data, clinical_data = simulate_data(n_samples=300, n_features=1000, n_clusters=N_CLUSTERS)

# --- Save to CSV files ---
expression_data_filename = 'simulated_expression_data.csv'
clinical_data_filename = 'simulated_clinical_data.csv'

expression_data.to_csv(expression_data_filename)
clinical_data.to_csv(clinical_data_filename)

print(f"\n🎉 Data generation successful!")
print(f"--> Expression data saved to: {expression_data_filename}")
print(f"--> Clinical data saved to: {clinical_data_filename}")
print("\nPlease confirm the files have been created before proceeding to the next step.")
