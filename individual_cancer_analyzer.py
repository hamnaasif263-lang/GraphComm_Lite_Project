"""
INDIVIDUAL CANCER TYPE ANALYZER
Separate analysis for each cancer: preprocessing → IGF pathway → dormancy → window of risk → FDA drugs
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from igf_pathway_analyzer import IGFPathwayAnalyzer
from fda_drug_database import FDADrugDatabase
from stage_2_clinical_risk_modeling import ClinicalRiskModeler

# Set publication-quality style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

class IndividualCancerAnalyzer:
    """
    Complete analysis pipeline for individual cancer types
    """
    
    def __init__(self, cancer_type: str, output_base: str = "output"):
        """
        Initialize analyzer for specific cancer type.
        
        Args:
            cancer_type: "breast_cancer", "lung_cancer", or "prostate_cancer"
            output_base: Base output directory
        """
        self.cancer_type = cancer_type
        self.output_base = output_base
        self.output_dir = os.path.join(output_base, cancer_type)
        
        # Subdirectories for each analysis stage
        self.dirs = {
            'preprocessing': os.path.join(self.output_dir, '01_preprocessing'),
            'igf_pathway': os.path.join(self.output_dir, '01_igf_pathway_analysis'),
            'dormancy': os.path.join(self.output_dir, '01_dormancy_signatures'),
            'risk': os.path.join(self.output_dir, '02_window_of_risk'),
            'drugs': os.path.join(self.output_dir, '03_drug_repurposing'),
            'figures': os.path.join(self.output_dir, 'figures'),
        }
        
        # Create directories
        for dir_path in self.dirs.values():
            os.makedirs(dir_path, exist_ok=True)
        
        # Initialize components
        self.igf_analyzer = IGFPathwayAnalyzer(cancer_type=cancer_type)
        self.risk_modeler = ClinicalRiskModeler(cancer_type=cancer_type)
        self.drug_db = FDADrugDatabase()
        
        # Data storage
        self.single_cell_data = None
        self.patient_metadata = None
        self.gene_expression = None
        self.pathway_scores = None
        self.dormancy_signature = None
        self.patient_features = None
        self.risk_scores = None
        self.drug_recommendations = None
        
    def generate_data(self, n_patients: int = 50, n_cells: int = 5000):
        """Generate realistic single-cell and patient data for cancer type."""
        np.random.seed(hash(self.cancer_type) % 2**32)
        
        # Single-cell data with required columns for risk modeling
        self.single_cell_data = pd.DataFrame({
            'patient_id': np.repeat(range(n_patients), n_cells),
            'cell_index': range(n_patients * n_cells),
            'cell_type': np.random.choice(['tumor', 'immune', 'fibroblast'], n_patients * n_cells),
            'dormancy_score': np.random.beta(2, 5, n_patients * n_cells),
            'igf_score': np.random.beta(3, 5, n_patients * n_cells),
            'proliferation_score': np.random.beta(3, 5, n_patients * n_cells),
        })
        
        # Patient metadata
        if self.cancer_type == "breast_cancer":
            relapse_rate = 0.35
        elif self.cancer_type == "lung_cancer":
            relapse_rate = 0.45
        else:  # prostate_cancer
            relapse_rate = 0.25
        
        self.patient_metadata = pd.DataFrame({
            'patient_id': range(n_patients),
            'early_relapse': np.random.binomial(1, relapse_rate, n_patients),
            'months_to_relapse': np.random.exponential(18, n_patients),
            'event_observed': np.random.binomial(1, 0.7, n_patients),
        })
        
        print(f"✓ Generated {n_patients} patients with {n_cells} cells each")
    
    def run_igf_pathway_analysis(self):
        """Analyze IGF pathway expression and activation."""
        print("\n2️⃣ IGF PATHWAY ANALYSIS")
        print("-" * 60)
        
        # Generate cancer-specific IGF expression
        self.gene_expression = self.igf_analyzer.generate_cancer_specific_expression(
            self.single_cell_data
        )
        print("   ✓ Generated cancer-specific IGF pathway expression")
        
        # Calculate pathway scores
        self.pathway_scores = self.igf_analyzer.calculate_pathway_scores()
        print("   ✓ Calculated pathway activation scores")
        
        # Get summary statistics
        summary_df = self.igf_analyzer.get_summary_statistics()
        summary_df.to_csv(
            os.path.join(self.dirs['igf_pathway'], f"{self.cancer_type}_igf_summary.csv"),
            index=False
        )
        
        print("\n   📊 IGF Pathway Summary:")
        for col in summary_df.columns[1:]:
            val = summary_df[col].values[0]
            print(f"      {col}: {val:.3f}")
        
        return summary_df
    
    def identify_dormancy_signatures(self):
        """Identify and characterize dormancy signatures."""
        print("\n3️⃣ DORMANCY SIGNATURE IDENTIFICATION")
        print("-" * 60)
        
        self.dormancy_signature = self.igf_analyzer.identify_dormancy_signature()
        
        sig = self.dormancy_signature
        print(f"   Dormant Cells:      {sig['dormant_pct']:.1f}%")
        print(f"   Proliferative Cells: {sig['proliferative_pct']:.1f}%")
        print(f"   Transitional Cells: {sig['transitional_pct']:.1f}%")
        print(f"   IGF Activation:     {sig['mean_igf_activation']:.3f}")
        print(f"   Dormancy Score:     {sig['mean_dormancy']:.3f}")
        
        # Save dormancy signature
        sig_df = pd.DataFrame([sig])
        sig_df.to_csv(
            os.path.join(self.dirs['dormancy'], f"{self.cancer_type}_dormancy_signature.csv"),
            index=False
        )
        
        return sig
    
    def run_window_of_risk_analysis(self):
        """Calculate patient-specific window of risk."""
        print("\n4️⃣ WINDOW OF RISK ANALYSIS")
        print("-" * 60)
        
        # Aggregate patient features
        self.patient_features = self.risk_modeler.aggregate_patient_features(
            self.single_cell_data,
            self.patient_metadata
        )
        print(f"   ✓ Aggregated features for {len(self.patient_features)} patients")
        
        # Fit regressors
        self.risk_modeler.fit_sequential_regressors()
        
        # Kaplan-Meier analysis
        kmf_h, kmf_l, results = self.risk_modeler.fit_kaplan_meier()
        
        # Calculate window of risk
        self.risk_scores = self.risk_modeler.calculate_window_of_risk(
            windows_months=[6, 12, 24]
        )
        
        # Save risk scores
        self.risk_scores.to_csv(
            os.path.join(self.dirs['risk'], f"{self.cancer_type}_risk_scores.csv"),
            index=False
        )
        
        print(f"   ✓ Risk scores calculated")
        
        return self.risk_scores
    
    def suggest_fda_drugs(self):
        """Get FDA drug recommendations based on IGF activation."""
        print("\n5️⃣ FDA DRUG RECOMMENDATIONS")
        print("-" * 60)
        
        if self.dormancy_signature is None:
            raise ValueError("Must run dormancy analysis first")
        
        igf_level = self.dormancy_signature['mean_igf_activation']
        
        self.drug_recommendations = self.drug_db.recommend_drugs(
            cancer_type=self.cancer_type,
            igf_activation_level=max(0, min(1, igf_level)),  # Clip to 0-1
            target_dormancy=True
        )
        
        print(f"\n   Top 10 Recommended Drugs (IGF Activation: {igf_level:.2f}):")
        print("   " + "-" * 76)
        
        for idx, row in self.drug_recommendations.head(10).iterrows():
            print(f"   {idx+1}. {row['drug_name']:<30} | Target: {row['target']:<18} | Score: {row['score']:.1f}")
        
        # Save recommendations
        self.drug_recommendations.to_csv(
            os.path.join(self.dirs['drugs'], f"{self.cancer_type}_drug_recommendations.csv"),
            index=False
        )
        
        return self.drug_recommendations
    
    def generate_igf_expression_heatmap(self):
        """Generate heatmap of IGF pathway expression by cancer type."""
        if self.gene_expression is None:
            return
        
        # Select representative IGF pathway genes
        igf_genes = ['IGF1', 'IGF2', 'IGF1R', 'AKT1', 'MTOR', 'MYC', 'CDKN1B']
        available_genes = [g for g in igf_genes if g in self.gene_expression.columns]
        
        if len(available_genes) < 3:
            return
        
        # Sample cells for visualization
        sample_size = min(500, len(self.gene_expression))
        sample_indices = np.random.choice(len(self.gene_expression), sample_size, replace=False)
        expr_sample = self.gene_expression.iloc[sample_indices][available_genes]
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Normalize by gene
        expr_normalized = (expr_sample - expr_sample.mean()) / expr_sample.std()
        
        sns.heatmap(
            expr_normalized.T,
            cmap='RdYlBu_r',
            ax=ax,
            cbar_kws={'label': 'Normalized Expression'},
            xticklabels=False,
            yticklabels=True
        )
        
        ax.set_title(f'{self.cancer_type.upper().replace("_", " ")}\nIGF Pathway Gene Expression Heatmap',
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Cells (n=500 sampled)', fontsize=12)
        ax.set_ylabel('Genes', fontsize=12)
        
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.dirs['figures'], f"{self.cancer_type}_igf_expression_heatmap.png"),
            dpi=300, bbox_inches='tight'
        )
        plt.close()
        
        print(f"   ✓ Saved: IGF expression heatmap")
    
    def generate_pathway_activation_distribution(self):
        """Generate distribution plots of pathway activation scores."""
        if self.pathway_scores is None:
            return
        
        pathway_cols = [c for c in self.pathway_scores.columns if c.endswith('_score')]
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        colors = plt.cm.Set2(np.linspace(0, 1, len(pathway_cols)))
        
        for idx, col in enumerate(pathway_cols[:6]):
            ax = axes[idx]
            
            data = self.pathway_scores[col].dropna()
            
            ax.hist(data, bins=30, alpha=0.7, color=colors[idx], edgecolor='black')
            ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {data.mean():.2f}')
            ax.set_xlabel('Activation Score', fontsize=11)
            ax.set_ylabel('Number of Cells', fontsize=11)
            ax.set_title(col.replace('_score', '').upper(), fontsize=12, fontweight='bold')
            ax.legend()
            ax.grid(alpha=0.3)
        
        fig.suptitle(f'{self.cancer_type.upper().replace("_", " ")}\nPathway Activation Score Distribution',
                    fontsize=14, fontweight='bold', y=1.00)
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.dirs['figures'], f"{self.cancer_type}_pathway_activation_distribution.png"),
            dpi=300, bbox_inches='tight'
        )
        plt.close()
        
        print(f"   ✓ Saved: Pathway activation distribution")
    
    def generate_dormancy_stacked_bar(self):
        """Generate stacked bar chart of cell populations."""
        if self.dormancy_signature is None:
            return
        
        sig = self.dormancy_signature
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Stacked bar chart
        categories = ['Dormant', 'Transitional', 'Proliferative']
        values = [sig['dormant_pct'], sig['transitional_pct'], sig['proliferative_pct']]
        colors = ['#FF6B6B', '#FFA500', '#4ECDC4']
        
        ax1.bar(['Cell Population'], values, color=colors, edgecolor='black', linewidth=2)
        ax1.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
        ax1.set_title(f'{self.cancer_type.upper().replace("_", " ")}\nCell Population Distribution',
                     fontsize=12, fontweight='bold')
        ax1.set_ylim([0, 100])
        
        for i, (cat, val) in enumerate(zip(categories, values)):
            ax1.text(0, sum(values[:i]) + val/2, f'{cat}\n{val:.1f}%',
                    ha='center', va='center', fontsize=11, fontweight='bold', color='black')
        
        ax1.legend(categories, loc='upper left', bbox_to_anchor=(0, 1), frameon=True)
        
        # Key metrics
        ax2.axis('off')
        
        metrics_text = f"""
        DORMANCY SIGNATURE METRICS
        
        Dormant Cell %:        {sig['dormant_pct']:.1f}%
        Transitional Cell %:   {sig['transitional_pct']:.1f}%
        Proliferative Cell %:  {sig['proliferative_pct']:.1f}%
        
        Mean IGF Activation:   {sig['mean_igf_activation']:.3f}
        Mean Dormancy Score:   {sig['mean_dormancy']:.3f}
        Mean Proliferation:    {sig['mean_proliferation']:.3f}
        
        Cancer Type:           {sig['cancer_type'].replace('_', ' ').title()}
        """
        
        ax2.text(0.1, 0.5, metrics_text, fontsize=12, family='monospace',
                verticalalignment='center', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.dirs['figures'], f"{self.cancer_type}_dormancy_distribution.png"),
            dpi=300, bbox_inches='tight'
        )
        plt.close()
        
        print(f"   ✓ Saved: Dormancy distribution")
    
    def generate_window_of_risk_plots(self):
        """Generate window of risk visualization."""
        if self.risk_scores is None:
            return
        
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        for idx, window in enumerate([6, 12, 24]):
            col = f'risk_{window}mo'
            ax = axes[idx]
            
            data = self.risk_scores[col]
            
            ax.hist(data, bins=20, alpha=0.7, color='steelblue', edgecolor='black')
            ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2, 
                      label=f'Mean: {data.mean():.1%}')
            ax.axvline(data.median(), color='green', linestyle='--', linewidth=2,
                      label=f'Median: {data.median():.1%}')
            
            ax.set_xlabel('Risk Probability', fontsize=11)
            ax.set_ylabel('Number of Patients', fontsize=11)
            ax.set_title(f'{window}-Month Window', fontsize=12, fontweight='bold')
            ax.set_xlim([0, 1])
            ax.legend()
            ax.grid(alpha=0.3)
        
        fig.suptitle(f'{self.cancer_type.upper().replace("_", " ")}\nWindow of Relapse Risk',
                    fontsize=14, fontweight='bold', y=1.00)
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.dirs['figures'], f"{self.cancer_type}_window_of_risk.png"),
            dpi=300, bbox_inches='tight'
        )
        plt.close()
        
        print(f"   ✓ Saved: Window of risk")
    
    def generate_top_drugs_visualization(self):
        """Generate top FDA drug recommendations visualization."""
        if self.drug_recommendations is None:
            return
        
        top_drugs = self.drug_recommendations.head(10).copy()
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        y_pos = np.arange(len(top_drugs))
        colors_map = {'FDA APPROVED': '#2ECC71', 'Clinical Trial': '#F39C12', 'Experimental': '#E74C3C'}
        colors = [colors_map.get(status, '#95A5A6') for status in top_drugs['fda_status']]
        
        bars = ax.barh(y_pos, top_drugs['score'], color=colors, edgecolor='black', linewidth=1.5)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_drugs['drug_name'], fontsize=11)
        ax.invert_yaxis()
        ax.set_xlabel('Recommendation Score', fontsize=12, fontweight='bold')
        ax.set_title(f'{self.cancer_type.upper().replace("_", " ")}\nTop 10 FDA-Approved Drugs for IGF Pathway Inhibition',
                    fontsize=13, fontweight='bold', pad=20)
        
        # Add value labels
        for i, (bar, score) in enumerate(zip(bars, top_drugs['score'])):
            ax.text(score + 0.1, bar.get_y() + bar.get_height()/2, f'{score:.1f}',
                   va='center', fontsize=10, fontweight='bold')
        
        # Legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#2ECC71', edgecolor='black', label='FDA Approved'),
            Patch(facecolor='#F39C12', edgecolor='black', label='Clinical Trial'),
            Patch(facecolor='#E74C3C', edgecolor='black', label='Experimental'),
        ]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
        
        ax.grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.dirs['figures'], f"{self.cancer_type}_top_drugs.png"),
            dpi=300, bbox_inches='tight'
        )
        plt.close()
        
        print(f"   ✓ Saved: Top drugs visualization")
    
    def run_complete_analysis(self):
        """Execute complete analysis pipeline for this cancer type."""
        print(f"\n\n{'='*70}")
        print(f"🔬 INDIVIDUAL CANCER ANALYSIS: {self.cancer_type.upper().replace('_', ' ')}")
        print(f"{'='*70}")
        
        print("\n1️⃣ PREPROCESSING & DATA GENERATION")
        print("-" * 60)
        self.generate_data(n_patients=50, n_cells=5000)
        
        # Run IGF pathway analysis
        igf_summary = self.run_igf_pathway_analysis()
        
        # Identify dormancy signatures
        dormancy_sig = self.identify_dormancy_signatures()
        
        # Window of risk
        risk_scores = self.run_window_of_risk_analysis()
        print(f"   ✓ Saved risk scores")
        
        # FDA drug recommendations
        drugs = self.suggest_fda_drugs()
        
        # Generate all figures
        print("\n6️⃣ GENERATING PUBLICATION-QUALITY FIGURES")
        print("-" * 60)
        
        self.generate_igf_expression_heatmap()
        self.generate_pathway_activation_distribution()
        self.generate_dormancy_stacked_bar()
        self.generate_window_of_risk_plots()
        self.generate_top_drugs_visualization()
        
        print(f"\n✅ COMPLETE ANALYSIS DONE FOR {self.cancer_type.upper()}")
        print(f"   📁 Results saved to: {self.output_dir}/")
        
        return {
            'igf_summary': igf_summary,
            'dormancy_signature': dormancy_sig,
            'risk_scores': risk_scores,
            'drug_recommendations': drugs
        }


if __name__ == "__main__":
    cancer_types = ["breast_cancer", "lung_cancer", "prostate_cancer"]
    
    for cancer in cancer_types:
        analyzer = IndividualCancerAnalyzer(cancer_type=cancer)
        results = analyzer.run_complete_analysis()
