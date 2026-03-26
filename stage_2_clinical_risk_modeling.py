"""
STAGE 2 — CLINICAL RISK MODELING
Survival Analysis & Window of Risk Estimation for Cancer Relapse

Pipeline:
  1. Patient-Level Feature Aggregation
     - % dormant cells
     - IGF activation score
     - Transitional-cell burden
  
  2. Sequential Dependent Regressors
     - Logistic Regression 1: Predict dormancy status
     - Logistic Regression 2: Predict risk given dormancy
  
  3. Window of Risk Estimation
     - Probability of relapse within 6mo, 12mo, 24mo
"""

import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Tuple, List
import warnings
warnings.filterwarnings('ignore')


class ClinicalRiskModeler:
    """
    Aggregates single-cell features into patient-level risk scores.
    Implements sequential logistic regressors and survival analysis.
    """
    
    def __init__(self, cancer_type: str = "breast"):
        """
        Initialize the clinical risk modeler.
        
        Args:
            cancer_type: "breast", "lung", or "prostate"
        """
        self.cancer_type = cancer_type
        self.patient_features = None
        self.lr1_model = None  # Dormancy predictor
        self.lr2_model = None  # Risk predictor
        self.kmf = KaplanMeierFitter()
        self.cph = CoxPHFitter()
        self.scaler = StandardScaler()
        self.risk_windows = {6: None, 12: None, 24: None}  # months
        
    def aggregate_patient_features(
        self,
        single_cell_data: pd.DataFrame,
        patient_metadata: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Aggregate single-cell level features to patient level.
        
        Args:
            single_cell_data: DataFrame with columns:
                - 'patient_id'
                - 'cell_type'
                - 'igf_score' (IGF pathway activation)
                - 'dormancy_score' (from Stage 1)
                - 'proliferation_score'
            
            patient_metadata: DataFrame with columns:
                - 'patient_id'
                - 'early_relapse' (0/1)
                - 'months_to_relapse' (event time)
                - 'event_observed' (0/1, censoring indicator)
        
        Returns:
            Aggregated patient-level feature matrix
        """
        patient_features = []
        
        for patient_id in patient_metadata['patient_id'].unique():
            # Filter single-cell data for this patient
            patient_cells = single_cell_data[
                single_cell_data['patient_id'] == patient_id
            ]
            
            if len(patient_cells) == 0:
                continue
            
            # Feature 1: % dormant cells (threshold: dormancy_score > 0.7)
            pct_dormant = (
                (patient_cells['dormancy_score'] > 0.7).sum() 
                / len(patient_cells)
            ) * 100
            
            # Feature 2: IGF pathway activation score (mean)
            igf_activation = patient_cells['igf_score'].mean()
            
            # Feature 3: IGF pathway activation score (std)
            igf_heterogeneity = patient_cells['igf_score'].std()
            
            # Feature 4: Transitional-cell burden
            # (high dormancy + high proliferation = bridging cells)
            transitional_mask = (
                (patient_cells['dormancy_score'] > 0.5) & 
                (patient_cells['proliferation_score'] > 0.6)
            )
            transitional_burden = (transitional_mask.sum() / len(patient_cells)) * 100
            
            # Feature 5: Proportion of high IGF cells
            high_igf_cells = (patient_cells['igf_score'] > 0.7).sum()
            pct_high_igf = (high_igf_cells / len(patient_cells)) * 100
            
            # Get clinical outcome
            patient_meta = patient_metadata[
                patient_metadata['patient_id'] == patient_id
            ].iloc[0]
            
            patient_features.append({
                'patient_id': patient_id,
                'pct_dormant_cells': pct_dormant,
                'igf_activation_score': igf_activation,
                'igf_heterogeneity': igf_heterogeneity,
                'transitional_cell_burden': transitional_burden,
                'pct_high_igf': pct_high_igf,
                'total_cells': len(patient_cells),
                'early_relapse': patient_meta['early_relapse'],
                'months_to_relapse': patient_meta['months_to_relapse'],
                'event_observed': patient_meta['event_observed']
            })
        
        self.patient_features = pd.DataFrame(patient_features)
        return self.patient_features
    
    def fit_sequential_regressors(self, test_size: float = 0.2):
        """
        Fit two sequential logistic regressors:
        
        LR1: Predict high dormancy status (pct_dormant_cells > median)
        LR2: Predict early relapse given dormancy status
        
        Args:
            test_size: Fraction of data to use for testing
        """
        if self.patient_features is None:
            raise ValueError("Must call aggregate_patient_features first")
        
        df = self.patient_features.copy()
        
        # Define high dormancy as > median
        df['high_dormancy'] = (
            df['pct_dormant_cells'] > df['pct_dormant_cells'].median()
        ).astype(int)
        
        # Feature matrix (exclude outcome variables)
        features_lr1 = [
            'igf_activation_score',
            'igf_heterogeneity',
            'transitional_cell_burden',
            'pct_high_igf'
        ]
        
        X_lr1 = df[features_lr1].values
        y_lr1 = df['high_dormancy'].values
        
        # Standardize features
        X_lr1_scaled = self.scaler.fit_transform(X_lr1)
        
        # LR1: Train dormancy predictor
        self.lr1_model = LogisticRegression(random_state=42, max_iter=1000)
        self.lr1_model.fit(X_lr1_scaled, y_lr1)
        
        # Get LR1 predictions (dormancy probability)
        df['dormancy_prob'] = self.lr1_model.predict_proba(X_lr1_scaled)[:, 1]
        
        # LR2: Include dormancy probability as feature
        features_lr2 = features_lr1 + ['dormancy_prob']
        X_lr2 = df[features_lr2].values
        y_lr2 = df['early_relapse'].values
        
        X_lr2_scaled = StandardScaler().fit_transform(X_lr2)
        
        # LR2: Train risk predictor
        self.lr2_model = LogisticRegression(random_state=42, max_iter=1000)
        self.lr2_model.fit(X_lr2_scaled, y_lr2)
        
        # Get LR2 predictions (risk probability)
        df['risk_probability'] = self.lr2_model.predict_proba(X_lr2_scaled)[:, 1]
        
        self.patient_features = df
        
        print(f"✓ LR1 (Dormancy Predictor) - Accuracy: {self.lr1_model.score(X_lr1_scaled, y_lr1):.3f}")
        print(f"✓ LR2 (Risk Predictor) - Accuracy: {self.lr2_model.score(X_lr2_scaled, y_lr2):.3f}")
        
        return self.patient_features
    
    def fit_kaplan_meier(self):
        """
        Fit Kaplan-Meier survival curves stratified by risk groups.
        """
        if self.patient_features is None:
            raise ValueError("Must call aggregate_patient_features first")
        
        df = self.patient_features.copy()
        
        # Define risk groups (high/low based on median risk probability)
        df['risk_group'] = (
            df['risk_probability'] > df['risk_probability'].median()
        ).astype(int)
        df['risk_group_label'] = df['risk_group'].map({0: 'Low Risk', 1: 'High Risk'})
        
        # Fit KM curves for high vs low risk
        high_risk = df[df['risk_group'] == 1]
        low_risk = df[df['risk_group'] == 0]
        
        # Fit KM for high risk
        kmf_high = KaplanMeierFitter()
        kmf_high.fit(
            high_risk['months_to_relapse'],
            event_observed=high_risk['event_observed'],
            label='High Risk'
        )
        
        # Fit KM for low risk
        kmf_low = KaplanMeierFitter()
        kmf_low.fit(
            low_risk['months_to_relapse'],
            event_observed=low_risk['event_observed'],
            label='Low Risk'
        )
        
        # Log-rank test
        results = logrank_test(
            high_risk['months_to_relapse'],
            low_risk['months_to_relapse'],
            event_observed_A=high_risk['event_observed'],
            event_observed_B=low_risk['event_observed']
        )
        
        print(f"\n📊 Kaplan-Meier Analysis ({self.cancer_type}):")
        print(f"   Log-rank test p-value: {results.p_value:.4f}")
        print(f"   Significant: {'Yes ✓' if results.p_value < 0.05 else 'No'}")
        
        self.kmf_high = kmf_high
        self.kmf_low = kmf_low
        self.patient_features = df
        
        return kmf_high, kmf_low, results
    
    def calculate_window_of_risk(self, windows_months: List[int] = [6, 12, 24]):
        """
        Calculate probability of relapse within each time window.
        
        Args:
            windows_months: List of time windows (e.g., [6, 12, 24])
        
        Returns:
            DataFrame with risk probabilities for each window
        """
        if self.patient_features is None:
            raise ValueError("Must call aggregate_patient_features first")
        
        df = self.patient_features.copy()
        
        # For each patient and each time window, calculate risk
        risk_windows = []
        
        for patient_id in df['patient_id'].unique():
            patient_row = df[df['patient_id'] == patient_id].iloc[0]
            
            risk_dict = {'patient_id': patient_id}
            
            for window in windows_months:
                # Risk = LR2 probability * adjustment for time window
                # Longer window → higher probability (simple model)
                base_risk = patient_row['risk_probability']
                
                # Time-decay adjustment: early relapse more likely in shorter windows
                time_adjustment = 1 - np.exp(-0.05 * window)
                
                window_risk = base_risk * time_adjustment
                
                risk_dict[f'risk_{window}mo'] = window_risk
            
            # Clinical outcome for reference
            risk_dict['early_relapse'] = patient_row['early_relapse']
            risk_dict['months_to_relapse'] = patient_row['months_to_relapse']
            
            risk_windows.append(risk_dict)
        
        risk_df = pd.DataFrame(risk_windows)
        self.risk_windows_df = risk_df
        
        return risk_df
    
    def print_risk_summary(self):
        """
        Print summary statistics for the window of risk.
        """
        if not hasattr(self, 'risk_windows_df'):
            raise ValueError("Must call calculate_window_of_risk first")
        
        df = self.risk_windows_df
        
        print(f"\n🎯 Window of Risk Summary ({self.cancer_type}):")
        print("=" * 60)
        
        for window in [6, 12, 24]:
            col = f'risk_{window}mo'
            mean_risk = df[col].mean()
            median_risk = df[col].median()
            std_risk = df[col].std()
            
            high_risk_count = (df[col] > 0.5).sum()
            
            print(f"\n📅 {window}-Month Window:")
            print(f"   Mean Risk:     {mean_risk:.1%}")
            print(f"   Median Risk:   {median_risk:.1%}")
            print(f"   Std Dev:       {std_risk:.1%}")
            print(f"   High Risk Patients (>50%): {high_risk_count}/{len(df)}")
    
    def plot_risk_distribution(self, save_path: str = None):
        """
        Plot distribution of risk across time windows.
        """
        if not hasattr(self, 'risk_windows_df'):
            raise ValueError("Must call calculate_window_of_risk first")
        
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        
        for idx, window in enumerate([6, 12, 24]):
            col = f'risk_{window}mo'
            ax = axes[idx]
            
            ax.hist(
                self.risk_windows_df[col],
                bins=15,
                alpha=0.7,
                color='steelblue',
                edgecolor='black'
            )
            ax.axvline(
                self.risk_windows_df[col].mean(),
                color='red',
                linestyle='--',
                linewidth=2,
                label=f'Mean: {self.risk_windows_df[col].mean():.1%}'
            )
            ax.set_xlabel('Risk Probability')
            ax.set_ylabel('Number of Patients')
            ax.set_title(f'{window}-Month Window\n{self.cancer_type.title()}')
            ax.legend()
            ax.grid(alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Risk distribution plot saved to {save_path}")
        
        plt.show()
    
    def save_patient_risk_scores(self, output_path: str):
        """
        Save patient-level risk scores to CSV.
        """
        if not hasattr(self, 'risk_windows_df'):
            raise ValueError("Must call calculate_window_of_risk first")
        
        output_df = self.risk_windows_df.copy()
        output_df.to_csv(output_path, index=False)
        print(f"✓ Patient risk scores saved to {output_path}")


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    
    # Create synthetic example data (replace with real data)
    np.random.seed(42)
    
    n_patients = 50
    n_cells_per_patient = 5000
    
    # Single-cell data
    single_cell_data = pd.DataFrame({
        'patient_id': np.repeat(range(n_patients), n_cells_per_patient),
        'cell_type': np.random.choice(['tumor', 'immune'], n_patients * n_cells_per_patient),
        'dormancy_score': np.random.beta(2, 5, n_patients * n_cells_per_patient),
        'igf_score': np.random.beta(3, 5, n_patients * n_cells_per_patient),
        'proliferation_score': np.random.beta(3, 5, n_patients * n_cells_per_patient),
    })
    
    # Patient metadata
    patient_metadata = pd.DataFrame({
        'patient_id': range(n_patients),
        'early_relapse': np.random.binomial(1, 0.4, n_patients),
        'months_to_relapse': np.random.exponential(18, n_patients),
        'event_observed': np.random.binomial(1, 0.7, n_patients)
    })
    
    # Initialize modeler
    modeler = ClinicalRiskModeler(cancer_type="breast")
    
    # Pipeline execution
    print("STAGE 2 — CLINICAL RISK MODELING")
    print("=" * 60)
    
    print("\n1️⃣ Aggregating patient-level features...")
    patient_df = modeler.aggregate_patient_features(single_cell_data, patient_metadata)
    print(f"✓ Aggregated {len(patient_df)} patients")
    
    print("\n2️⃣ Fitting sequential regressors...")
    modeler.fit_sequential_regressors()
    
    print("\n3️⃣ Kaplan-Meier survival analysis...")
    kmf_h, kmf_l, results = modeler.fit_kaplan_meier()
    
    print("\n4️⃣ Calculating window of risk...")
    risk_df = modeler.calculate_window_of_risk(windows_months=[6, 12, 24])
    print(f"✓ Risk scores calculated for {len(risk_df)} patients")
    
    print("\n5️⃣ Risk Summary:")
    modeler.print_risk_summary()
    
    print("\n6️⃣ Plotting risk distribution...")
    modeler.plot_risk_distribution()
    
    print("\n7️⃣ Saving results...")
    modeler.save_patient_risk_scores(f"patient_risk_scores_{modeler.cancer_type}.csv")
    
    print("\n" + "=" * 60)
    print("✅ STAGE 2 COMPLETE")
