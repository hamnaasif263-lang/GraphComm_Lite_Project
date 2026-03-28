"""
UNIFIED PIPELINE RUNNER
Predicts window of risk for all three cancer types and organizes outputs
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
from stage_2_clinical_risk_modeling import ClinicalRiskModeler


def create_output_structure(base_dir: str = "output"):
    """
    Create the complete output folder structure:
    output/
    ├── breast_cancer/
    │   ├── 01_preprocessing/
    │   ├── 02_window_of_risk/
    │   └── 03_drug_repurposing/
    ├── lung_cancer/
    │   ├── 01_preprocessing/
    │   ├── 02_window_of_risk/
    │   └── 03_drug_repurposing/
    └── prostate_cancer/
        ├── 01_preprocessing/
        ├── 02_window_of_risk/
        └── 03_drug_repurposing/
    """
    cancer_types = ["breast_cancer", "lung_cancer", "prostate_cancer"]
    stages = ["01_preprocessing", "02_window_of_risk", "03_drug_repurposing"]
    
    paths = {}
    
    for cancer in cancer_types:
        cancer_dir = os.path.join(base_dir, cancer)
        os.makedirs(cancer_dir, exist_ok=True)
        
        paths[cancer] = {}
        for stage in stages:
            stage_dir = os.path.join(cancer_dir, stage)
            os.makedirs(stage_dir, exist_ok=True)
            paths[cancer][stage] = stage_dir
    
    print(f"✓ Output folder structure created at: {base_dir}/")
    return paths


def generate_synthetic_data(cancer_type: str, n_patients: int = 50, n_cells: int = 5000):
    """
    Generate synthetic single-cell and patient metadata.
    In production, this would load real GEO data.
    """
    np.random.seed(hash(cancer_type) % 2**32)
    
    # Single-cell data
    single_cell_data = pd.DataFrame({
        'patient_id': np.repeat(range(n_patients), n_cells),
        'cell_type': np.random.choice(['tumor', 'immune', 'fibroblast'], n_patients * n_cells),
        'dormancy_score': np.random.beta(2, 5, n_patients * n_cells),
        'igf_score': np.random.beta(3, 5, n_patients * n_cells),
        'proliferation_score': np.random.beta(3, 5, n_patients * n_cells),
        'myc_score': np.random.beta(2, 6, n_patients * n_cells),
        'cdkn1b_score': np.random.beta(4, 3, n_patients * n_cells),
    })
    
    # Cancer-specific adjustments
    if cancer_type == "breast_cancer":
        # Breast cancer: higher IGF sensitivity
        single_cell_data['igf_score'] *= 1.3
        relapse_prob = 0.35
    elif cancer_type == "lung_cancer":
        # Lung cancer: more heterogeneous
        single_cell_data['dormancy_score'] *= 1.2
        relapse_prob = 0.45
    else:  # prostate_cancer
        # Prostate cancer: slower progression
        single_cell_data['dormancy_score'] *= 0.8
        relapse_prob = 0.25
    
    # Patient metadata
    patient_metadata = pd.DataFrame({
        'patient_id': range(n_patients),
        'early_relapse': np.random.binomial(1, relapse_prob, n_patients),
        'months_to_relapse': np.random.exponential(18, n_patients),
        'event_observed': np.random.binomial(1, 0.7, n_patients),
        'cancer_type': cancer_type
    })
    
    return single_cell_data, patient_metadata


def run_stage2_pipeline(cancer_type: str, output_paths: dict):
    """
    Execute Stage 2 (Clinical Risk Modeling) for a single cancer type.
    """
    print(f"\n{'='*70}")
    print(f"🔬 STAGE 2: CLINICAL RISK MODELING - {cancer_type.upper().replace('_', ' ')}")
    print(f"{'='*70}")
    
    # Step 1: Generate/Load data
    print("\n1️⃣ Loading patient data...")
    single_cell_data, patient_metadata = generate_synthetic_data(cancer_type)
    print(f"   ✓ Loaded {len(patient_metadata)} patients, {len(single_cell_data)} cells")
    
    # Save raw data
    preprocessing_dir = output_paths[cancer_type]["01_preprocessing"]
    single_cell_data.to_csv(
        os.path.join(preprocessing_dir, f"{cancer_type}_single_cell_data.csv"),
        index=False
    )
    patient_metadata.to_csv(
        os.path.join(preprocessing_dir, f"{cancer_type}_patient_metadata.csv"),
        index=False
    )
    print(f"   ✓ Saved preprocessing data to {preprocessing_dir}/")
    
    # Step 2: Initialize modeler
    print("\n2️⃣ Aggregating patient-level features...")
    modeler = ClinicalRiskModeler(cancer_type=cancer_type)
    patient_features = modeler.aggregate_patient_features(single_cell_data, patient_metadata)
    print(f"   ✓ Aggregated features for {len(patient_features)} patients")
    
    # Save aggregated features
    patient_features.to_csv(
        os.path.join(preprocessing_dir, f"{cancer_type}_aggregated_features.csv"),
        index=False
    )
    
    # Step 3: Fit sequential regressors
    print("\n3️⃣ Fitting sequential logistic regressors...")
    modeler.fit_sequential_regressors()
    
    # Step 4: Kaplan-Meier analysis
    print("\n4️⃣ Performing Kaplan-Meier survival analysis...")
    kmf_high, kmf_low, results = modeler.fit_kaplan_meier()
    print(f"   ✓ Log-rank test p-value: {results.p_value:.4f}")
    
    # Step 5: Calculate window of risk
    print("\n5️⃣ Calculating window of risk (6mo, 12mo, 24mo)...")
    risk_df = modeler.calculate_window_of_risk(windows_months=[6, 12, 24])
    print(f"   ✓ Calculated risk scores for {len(risk_df)} patients")
    
    # Save window of risk
    window_risk_dir = output_paths[cancer_type]["02_window_of_risk"]
    risk_df.to_csv(
        os.path.join(window_risk_dir, f"{cancer_type}_window_of_risk.csv"),
        index=False
    )
    print(f"   ✓ Saved window of risk to {window_risk_dir}/")
    
    # Print risk summary
    print("\n6️⃣ Risk Summary:")
    modeler.print_risk_summary()
    
    # Save risk summary statistics
    summary_stats = []
    for window in [6, 12, 24]:
        col = f'risk_{window}mo'
        summary_stats.append({
            'cancer_type': cancer_type,
            'window_months': window,
            'mean_risk': risk_df[col].mean(),
            'median_risk': risk_df[col].median(),
            'std_risk': risk_df[col].std(),
            'high_risk_count': (risk_df[col] > 0.5).sum(),
            'total_patients': len(risk_df)
        })
    
    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(
        os.path.join(window_risk_dir, f"{cancer_type}_risk_summary_stats.csv"),
        index=False
    )
    
    # Step 7: Generate detailed patient report
    print("\n7️⃣ Generating patient-level risk report...")
    patient_report = risk_df.copy()
    patient_report['cancer_type'] = cancer_type
    
    # Categorize risk levels
    def categorize_risk(prob):
        if prob < 0.15:
            return 'Low Risk'
        elif prob < 0.35:
            return 'Intermediate Risk'
        else:
            return 'High Risk'
    
    patient_report['risk_6mo_category'] = patient_report['risk_6mo'].apply(categorize_risk)
    patient_report['risk_12mo_category'] = patient_report['risk_12mo'].apply(categorize_risk)
    patient_report['risk_24mo_category'] = patient_report['risk_24mo'].apply(categorize_risk)
    
    patient_report.to_csv(
        os.path.join(window_risk_dir, f"{cancer_type}_patient_risk_report.csv"),
        index=False
    )
    print(f"   ✓ Saved patient risk report to {window_risk_dir}/")
    
    return modeler, risk_df, patient_report, summary_df


def run_full_pipeline():
    """
    Execute the complete pipeline for all three cancer types.
    """
    print("\n" + "="*70)
    print("🚀 GRAPHCOMM-LITE: UNIFIED THREE-CANCER PIPELINE")
    print("="*70)
    
    # Step 1: Create output structure
    print("\n📁 Creating output folder structure...")
    output_paths = create_output_structure(base_dir="output")
    
    # Step 2: Run Stage 2 for each cancer type
    cancer_types = ["breast_cancer", "lung_cancer", "prostate_cancer"]
    all_results = {}
    
    for cancer_type in cancer_types:
        try:
            modeler, risk_df, patient_report, summary_df = run_stage2_pipeline(
                cancer_type,
                output_paths
            )
            all_results[cancer_type] = {
                'modeler': modeler,
                'risk_df': risk_df,
                'patient_report': patient_report,
                'summary_df': summary_df
            }
        except Exception as e:
            print(f"\n❌ Error processing {cancer_type}: {str(e)}")
            continue
    
    # Step 3: Generate comparative analysis
    print("\n\n" + "="*70)
    print("📊 COMPARATIVE ANALYSIS ACROSS CANCER TYPES")
    print("="*70)
    
    comparative_summary = []
    for cancer_type in all_results.keys():
        summary_df = all_results[cancer_type]['summary_df']
        for _, row in summary_df.iterrows():
            comparative_summary.append(row)
    
    if comparative_summary:
        comparative_df = pd.DataFrame(comparative_summary)
        
        # Save comparative summary
        comparative_df.to_csv(
            os.path.join("output", "COMPARATIVE_RISK_SUMMARY.csv"),
            index=False
        )
        
        print("\n📈 Cross-Cancer Risk Comparison (6-month window):")
        print("-" * 70)
        for cancer_type in cancer_types:
            data = comparative_df[
                (comparative_df['cancer_type'] == cancer_type) & 
                (comparative_df['window_months'] == 6)
            ]
            if len(data) > 0:
                row = data.iloc[0]
                print(f"\n{cancer_type.upper().replace('_', ' ')}:")
                print(f"   Mean Risk (6mo):      {row['mean_risk']:.1%}")
                print(f"   Median Risk (6mo):    {row['median_risk']:.1%}")
                print(f"   High Risk Patients:   {row['high_risk_count']:.0f}/{row['total_patients']:.0f}")
    
    # Step 4: Generate master patient registry
    print("\n\n📋 Generating master patient registry...")
    master_registry = []
    
    for cancer_type in all_results.keys():
        patient_report = all_results[cancer_type]['patient_report']
        master_registry.append(patient_report)
    
    if master_registry:
        master_df = pd.concat(master_registry, ignore_index=True)
        master_df.to_csv(
            os.path.join("output", "MASTER_PATIENT_REGISTRY.csv"),
            index=False
        )
        print(f"✓ Master registry created with {len(master_df)} patients across 3 cancer types")
    
    # Step 5: Summary statistics
    print("\n\n" + "="*70)
    print("✅ PIPELINE EXECUTION COMPLETE")
    print("="*70)
    
    print("\n📁 Output Structure:")
    print("""
    output/
    ├── breast_cancer/
    │   ├── 01_preprocessing/
    │   │   ├── breast_cancer_single_cell_data.csv
    │   │   ├── breast_cancer_patient_metadata.csv
    │   │   └── breast_cancer_aggregated_features.csv
    │   ├── 02_window_of_risk/
    │   │   ├── breast_cancer_window_of_risk.csv
    │   │   ├── breast_cancer_risk_summary_stats.csv
    │   │   └── breast_cancer_patient_risk_report.csv
    │   └── 03_drug_repurposing/
    │       └── [TO BE FILLED IN STAGE 3]
    │
    ├── lung_cancer/
    │   ├── 01_preprocessing/
    │   ├── 02_window_of_risk/
    │   └── 03_drug_repurposing/
    │
    ├── prostate_cancer/
    │   ├── 01_preprocessing/
    │   ├── 02_window_of_risk/
    │   └── 03_drug_repurposing/
    │
    ├── COMPARATIVE_RISK_SUMMARY.csv
    └── MASTER_PATIENT_REGISTRY.csv
    """)
    
    print("\n📊 Key Output Files:")
    print("   • COMPARATIVE_RISK_SUMMARY.csv - Risk metrics across all cancer types")
    print("   • MASTER_PATIENT_REGISTRY.csv - All patients with risk stratification")
    print("\n✓ All data organized by cancer type and pipeline stage")
    print("✓ Ready for Stage 3: IGF-Specific Drug Repurposing")
    
    return all_results, output_paths


if __name__ == "__main__":
    all_results, output_paths = run_full_pipeline()
    print("\n" + "="*70)
    print("🎯 STAGE 2 COMPLETE - All patients stratified by relapse risk")
    print("="*70)
