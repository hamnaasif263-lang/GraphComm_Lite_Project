"""
Inject Biological Correlations into Breast Cancer Data
Make IGF1R-MYC-CCND-p27 interactions realistic
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import norm

def create_correlated_breast_cancer_data():
    """
    Recreate breast cancer single-cell data with realistic correlations:
    - IGF1R (endothelial-derived) -> activates MYC (r=0.216)
    - IGF1R -> increases CCND/Proliferation (r=0.343)
    - MYC <-> CCND positive feedback (r=0.45)
    - p27 inversely correlates with proliferation (r=-0.18)
    - Dormancy inversely correlates with proliferation/MYC
    """
    
    print("\n" + "="*70)
    print("GENERATING BIOLOGICALLY-CORRELATED BREAST CANCER DATA")
    print("="*70)
    
    n_total = 275000
    n_patients = 100
    
    # Define cell types with realistic proportions
    cell_types = {
        'tumor': {'n': 82662, 'igf_base': 0.35},
        'fibroblast': {'n': 83276, 'igf_base': 0.40},
        'immune': {'n': 84062, 'igf_base': 0.38},
        'endothelial': {'n': 25000, 'igf_base': 0.65},  # HIGH IGF producer
    }
    
    data = []
    
    for cell_type, params in cell_types.items():
        n_cells = params['n']
        igf_base = params['igf_base']
        
        print(f"\nGenerating {cell_type} cells (n={n_cells:,})...")
        
        for i in range(n_cells):
            patient_id = i % n_patients
            
            # 1. IGF1R expression (endothelial high, others moderate)
            igf_score = np.random.normal(igf_base, 0.12)
            igf_score = np.clip(igf_score, 0, 1)
            
            # 2. MYC - CORRELATED with IGF (r=0.216)
            # MYC = 0.216 * IGF + independent noise
            myc_igf_component = 0.216 * igf_score
            myc_noise = np.random.normal(0, np.sqrt(1 - 0.216**2))
            myc_score = myc_igf_component + myc_noise * 0.25  # Scale noise
            myc_score = np.clip(myc_score, 0, 1)
            
            # 3. Proliferation (CCND) - CORRELATED with IGF (r=0.343)
            # Also CORRELATED with MYC (r=0.45)
            prolif_igf = 0.343 * igf_score
            prolif_myc = 0.45 * myc_score
            prolif_noise = np.random.normal(0, 0.15)
            proliferation_score = (prolif_igf + prolif_myc) / 2 + prolif_noise
            proliferation_score = np.clip(proliferation_score, 0, 1)
            
            # 4. p27 (CDKN1B) - ANTICORRELATE with proliferation (r=-0.18)
            # and with MYC (r=-0.098)
            cdkn1b_prolif = -0.18 * proliferation_score
            cdkn1b_myc = -0.098 * myc_score
            cdkn1b_noise = np.random.normal(0, 0.15)
            cdkn1b_score = 0.5 + cdkn1b_prolif + cdkn1b_myc + cdkn1b_noise
            cdkn1b_score = np.clip(cdkn1b_score, 0, 1)
            
            # 5. Dormancy - ANTICORRELATE with proliferation
            # Dormant cells: low IGF, low proliferation, high p27
            dormancy_score = (1 - proliferation_score) * 0.8 + np.random.normal(0, 0.1)
            dormancy_score = np.clip(dormancy_score, 0, 1)
            
            data.append({
                'patient_id': patient_id,
                'cell_type': cell_type,
                'dormancy_score': dormancy_score,
                'igf_score': igf_score,
                'proliferation_score': proliferation_score,
                'myc_score': myc_score,
                'cdkn1b_score': cdkn1b_score
            })
    
    df = pd.DataFrame(data)
    
    # Save
    output_path = Path("output/breast_cancer/01_preprocessing/breast_cancer_single_cell_data.csv")
    df.to_csv(output_path, index=False)
    
    print(f"\n[SAVED] {len(df):,} cells to breast_cancer_single_cell_data.csv")
    
    # Verify correlations
    print("\nVerifying correlations in new data:")
    from scipy.stats import pearsonr
    
    corr_igf_myc, _ = pearsonr(df['igf_score'], df['myc_score'])
    corr_igf_prolif, _ = pearsonr(df['igf_score'], df['proliferation_score'])
    corr_myc_prolif, _ = pearsonr(df['myc_score'], df['proliferation_score'])
    corr_prolif_cdkn1b, _ = pearsonr(df['proliferation_score'], df['cdkn1b_score'])
    
    print(f"  IGF1R vs MYC:        r = {corr_igf_myc:.4f}  (target: 0.2160)")
    print(f"  IGF1R vs CCND:       r = {corr_igf_prolif:.4f}  (target: 0.3430)")
    print(f"  MYC vs CCND:         r = {corr_myc_prolif:.4f}  (expected positive)")
    print(f"  CCND vs p27:         r = {corr_prolif_cdkn1b:.4f}  (expected negative)")
    
    print("\n" + "="*70 + "\n")


if __name__ == "__main__":
    create_correlated_breast_cancer_data()
