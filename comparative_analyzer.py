"""
COMPARATIVE ANALYSIS ACROSS CANCER TYPES
Compares IGF expression, dormancy, risk, and FDA drug recommendations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

class ComparativeAnalyzer:
    """
    Analytical comparison across breast, lung, and prostate cancer types
    """
    
    def __init__(self, output_base: str = "output"):
        """
        Initialize comparative analyzer.
        
        Args:
            output_base: Base output directory
        """
        self.output_base = output_base
        self.comparative_dir = os.path.join(output_base, "04_COMPARATIVE_ANALYSIS")
        os.makedirs(self.comparative_dir, exist_ok=True)
        
        self.cancer_types = ["breast_cancer", "lung_cancer", "prostate_cancer"]
        self.results = {}
        
    def load_cancer_results(self):
        """Load analysis results for all cancer types."""
        print("\n📂 Loading analysis results for all cancer types...")
        print("-" * 60)
        
        for cancer_type in self.cancer_types:
            cancer_dir = os.path.join(self.output_base, cancer_type)
            
            # Load IGF summary
            igf_file = os.path.join(cancer_dir, '01_igf_pathway_analysis', f"{cancer_type}_igf_summary.csv")
            if os.path.exists(igf_file):
                igf_summary = pd.read_csv(igf_file)
                self.results[cancer_type] = {'igf_summary': igf_summary}
            
            # Load dormancy signature
            dormancy_file = os.path.join(cancer_dir, '01_dormancy_signatures', f"{cancer_type}_dormancy_signature.csv")
            if os.path.exists(dormancy_file):
                dormancy = pd.read_csv(dormancy_file)
                if cancer_type in self.results:
                    self.results[cancer_type]['dormancy'] = dormancy
            
            # Load risk scores
            risk_file = os.path.join(cancer_dir, '02_window_of_risk', f"{cancer_type}_risk_scores.csv")
            if os.path.exists(risk_file):
                risk_scores = pd.read_csv(risk_file)
                if cancer_type in self.results:
                    self.results[cancer_type]['risk_scores'] = risk_scores
            
            # Load drug recommendations
            drug_file = os.path.join(cancer_dir, '03_drug_repurposing', f"{cancer_type}_drug_recommendations.csv")
            if os.path.exists(drug_file):
                drugs = pd.read_csv(drug_file)
                if cancer_type in self.results:
                    self.results[cancer_type]['drugs'] = drugs
            
            print(f"   ✓ Loaded {cancer_type}")
        
        return self.results
    
    def generate_igf_comparison_table(self):
        """Generate IGF pathway comparison table."""
        print("\n1️⃣  IGF PATHWAY COMPARISON")
        print("-" * 60)
        
        comparison_data = []
        
        for cancer_type in self.cancer_types:
            if cancer_type not in self.results or 'igf_summary' not in self.results[cancer_type]:
                continue
            
            igf = self.results[cancer_type]['igf_summary'].iloc[0]
            
            comparison_data.append({
                'Cancer_Type': cancer_type.replace('_', ' ').title(),
                'IGF_Ligand': igf['IGF_Ligand_Mean'],
                'IGF_Receptor': igf['IGF_Receptor_Mean'],
                'PI3K_AKT': igf['PI3K_AKT_Mean'],
                'mTOR': igf['mTOR_Mean'],
                'MAPK': igf['MAPK_Mean'],
                'MYC': igf['MYC_Mean'],
                'Overall_IGF_Activation': igf['IGF_Pathway_Activation'],
            })
        
        comparison_df = pd.DataFrame(comparison_data)
        comparison_df.to_csv(
            os.path.join(self.comparative_dir, 'IGF_Pathway_Comparison.csv'),
            index=False
        )
        
        print("\n   IGF Pathway Activation Levels:")
        print("   " + "-" * 56)
        for _, row in comparison_df.iterrows():
            print(f"   {row['Cancer_Type']:<20} | Overall IGF: {row['Overall_IGF_Activation']:.3f}")
        
        return comparison_df
    
    def generate_dormancy_comparison_table(self):
        """Generate dormancy signature comparison."""
        print("\n2️⃣ DORMANCY SIGNATURE COMPARISON")
        print("-" * 60)
        
        comparison_data = []
        
        for cancer_type in self.cancer_types:
            if cancer_type not in self.results or 'dormancy' not in self.results[cancer_type]:
                continue
            
            dorm = self.results[cancer_type]['dormancy'].iloc[0]
            
            comparison_data.append({
                'Cancer_Type': cancer_type.replace('_', ' ').title(),
                'Dormant_Pct': dorm['dormant_pct'],
                'Transitional_Pct': dorm['transitional_pct'],
                'Proliferative_Pct': dorm['proliferative_pct'],
                'Mean_IGF_Activation': dorm['mean_igf_activation'],
                'Mean_Dormancy_Score': dorm['mean_dormancy_score'],
                'Mean_Proliferation': dorm['mean_proliferation'],
            })
        
        comparison_df = pd.DataFrame(comparison_data)
        comparison_df.to_csv(
            os.path.join(self.comparative_dir, 'Dormancy_Signature_Comparison.csv'),
            index=False
        )
        
        print("\n   Cell Population Composition:")
        print("   " + "-" * 56)
        for _, row in comparison_df.iterrows():
            print(f"   {row['Cancer_Type']:<20}")
            print(f"      Dormant:       {row['Dormant_Pct']:>6.1f}%")
            print(f"      Transitional:  {row['Transitional_Pct']:>6.1f}%")
            print(f"      Proliferative: {row['Proliferative_Pct']:>6.1f}%")
        
        return comparison_df
    
    def generate_risk_comparison_table(self):
        """Generate window of risk comparison."""
        print("\n3️⃣ WINDOW OF RISK COMPARISON")
        print("-" * 60)
        
        comparison_data = []
        
        for cancer_type in self.cancer_types:
            if cancer_type not in self.results or 'risk_scores' not in self.results[cancer_type]:
                continue
            
            risk = self.results[cancer_type]['risk_scores']
            
            comparison_data.append({
                'Cancer_Type': cancer_type.replace('_', ' ').title(),
                'Risk_6mo_Mean': risk['risk_6mo'].mean(),
                'Risk_6mo_Median': risk['risk_6mo'].median(),
                'Risk_12mo_Mean': risk['risk_12mo'].mean(),
                'Risk_12mo_Median': risk['risk_12mo'].median(),
                'Risk_24mo_Mean': risk['risk_24mo'].mean(),
                'Risk_24mo_Median': risk['risk_24mo'].median(),
            })
        
        comparison_df = pd.DataFrame(comparison_data)
        comparison_df.to_csv(
            os.path.join(self.comparative_dir, 'Risk_Window_Comparison.csv'),
            index=False
        )
        
        print("\n   Relapse Risk Windows:")
        print("   " + "-" * 56)
        for _, row in comparison_df.iterrows():
            print(f"   {row['Cancer_Type']:<20}")
            print(f"      6-month:   {row['Risk_6mo_Mean']:.1%} (median: {row['Risk_6mo_Median']:.1%})")
            print(f"      12-month:  {row['Risk_12mo_Mean']:.1%} (median: {row['Risk_12mo_Median']:.1%})")
            print(f"      24-month:  {row['Risk_24mo_Mean']:.1%} (median: {row['Risk_24mo_Median']:.1%})")
        
        return comparison_df
    
    def generate_igf_pathway_comparison_plot(self, igf_comparison: pd.DataFrame):
        """Generate IGF pathway comparison heatmap."""
        if igf_comparison.empty:
            return
        
        # Create heatmap data
        heatmap_data = igf_comparison.set_index('Cancer_Type')[
            ['IGF_Ligand', 'IGF_Receptor', 'PI3K_AKT', 'mTOR', 'MAPK', 'MYC']
        ]
        
        # Normalize for visualization
        heatmap_normalized = (heatmap_data - heatmap_data.min()) / (heatmap_data.max() - heatmap_data.min())
        
        fig, ax = plt.subplots(figsize=(12, 5))
        
        sns.heatmap(
            heatmap_normalized.T,
            annot=heatmap_data.T,
            fmt='.2f',
            cmap='RdYlGn',
            ax=ax,
            cbar_kws={'label': 'Normalized Expression Level'},
            linewidths=2,
            linecolor='black'
        )
        
        ax.set_title('IGF Pathway Component Expression Across Cancer Types',
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('')
        ax.set_ylabel('Pathway Components', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.comparative_dir, 'COMPARATIVE_IGF_Pathway_Heatmap.png'),
            dpi=300, bbox_inches='tight'
        )
        plt.close()
        
        print("   ✓ Saved: IGF pathway heatmap")
    
    def generate_dormancy_comparison_plot(self, dormancy_comparison: pd.DataFrame):
        """Generate dormancy distribution comparison."""
        if dormancy_comparison.empty:
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Stacked bar chart
        ax1 = axes[0]
        x = np.arange(len(dormancy_comparison))
        width = 0.25
        
        dormancy_comparison_sorted = dormancy_comparison.sort_values('Dormant_Pct', ascending=False)
        
        ax1.bar(x, dormancy_comparison_sorted['Dormant_Pct'], width, 
               label='Dormant', color='#FF6B6B', edgecolor='black')
        ax1.bar(x, dormancy_comparison_sorted['Transitional_Pct'], width, 
               bottom=dormancy_comparison_sorted['Dormant_Pct'],
               label='Transitional', color='#FFA500', edgecolor='black')
        ax1.bar(x, dormancy_comparison_sorted['Proliferative_Pct'], width,
               bottom=dormancy_comparison_sorted['Dormant_Pct'] + dormancy_comparison_sorted['Transitional_Pct'],
               label='Proliferative', color='#4ECDC4', edgecolor='black')
        
        ax1.set_ylabel('Cell Percentage (%)', fontsize=12, fontweight='bold')
        ax1.set_title('Cell Population Composition\nAcross Cancer Types', 
                     fontsize=13, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels(dormancy_comparison_sorted['Cancer_Type'], rotation=0)
        ax1.legend(loc='upper right')
        ax1.set_ylim([0, 100])
        ax1.grid(axis='y', alpha=0.3)
        
        # IGF activation comparison
        ax2 = axes[1]
        colors = ['#3498DB', '#E74C3C', '#2ECC71']
        bars = ax2.bar(dormancy_comparison_sorted['Cancer_Type'], 
                       dormancy_comparison_sorted['Mean_IGF_Activation'],
                       color=colors, edgecolor='black', linewidth=2)
        
        ax2.set_ylabel('Mean IGF Pathway Activation', fontsize=12, fontweight='bold')
        ax2.set_title('IGF Pathway Activation\nby Cancer Type',
                     fontsize=13, fontweight='bold')
        ax2.set_ylim([0, max(dormancy_comparison_sorted['Mean_IGF_Activation']) * 1.2])
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.3f}',
                    ha='center', va='bottom', fontsize=11, fontweight='bold')
        
        ax2.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.comparative_dir, 'COMPARATIVE_Dormancy_Distribution.png'),
            dpi=300, bbox_inches='tight'
        )
        plt.close()
        
        print("   ✓ Saved: Dormancy distribution comparison")
    
    def generate_risk_comparison_plot(self, risk_comparison: pd.DataFrame):
        """Generate window of risk comparison."""
        if risk_comparison.empty:
            return
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        windows = [6, 12, 24]
        x = np.arange(len(risk_comparison))
        width = 0.25
        
        for idx, window in enumerate(windows):
            key = f'Risk_{window}mo_Mean'
            ax.bar(x + idx*width, risk_comparison[key], width, 
                  label=f'{window}-month', edgecolor='black', linewidth=1.5)
        
        ax.set_ylabel('Relapse Risk Probability', fontsize=12, fontweight='bold')
        ax.set_title('Window of Relapse Risk\nAcross Cancer Types',
                    fontsize=13, fontweight='bold')
        ax.set_xticks(x + width)
        ax.set_xticklabels(risk_comparison['Cancer_Type'], rotation=0)
        ax.legend(fontsize=11)
        ax.set_ylim([0, 0.4])
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.comparative_dir, 'COMPARATIVE_Window_of_Risk.png'),
            dpi=300, bbox_inches='tight'
        )
        plt.close()
        
        print("   ✓ Saved: Window of risk comparison")
    
    def generate_igf_vs_dormancy_scatter(self):
        """Generate IGF activation vs dormancy scatter plot."""
        scatter_data = []
        
        for cancer_type in self.cancer_types:
            if cancer_type not in self.results:
                continue
            
            if 'igf_summary' in self.results[cancer_type] and 'dormancy' in self.results[cancer_type]:
                igf = self.results[cancer_type]['igf_summary'].iloc[0]
                dorm = self.results[cancer_type]['dormancy'].iloc[0]
                
                scatter_data.append({
                    'Cancer_Type': cancer_type,
                    'IGF_Activation': igf['IGF_Pathway_Activation'],
                    'Dormancy_Pct': dorm['dormant_pct'],
                    'Cancer_Label': cancer_type.replace('_', ' ').title()
                })
        
        if not scatter_data:
            return
        
        scatter_df = pd.DataFrame(scatter_data)
        
        fig, ax = plt.subplots(figsize=(10, 7))
        
        colors = {'Breast Cancer': '#FF6B6B', 'Lung Cancer': '#FFA500', 'Prostate Cancer': '#4ECDC4'}
        
        for _, row in scatter_df.iterrows():
            ax.scatter(row['IGF_Activation'], row['Dormancy_Pct'], s=1000,
                      color=colors[row['Cancer_Label']], alpha=0.7,
                      edgecolor='black', linewidth=2, label=row['Cancer_Label'])
        
        ax.set_xlabel('IGF Pathway Activation', fontsize=12, fontweight='bold')
        ax.set_ylabel('Dormant Cell Percentage (%)', fontsize=12, fontweight='bold')
        ax.set_title('IGF Pathway Activation vs. Dormancy Status\nAcross Cancer Types',
                    fontsize=13, fontweight='bold')
        ax.legend(fontsize=11, loc='best')
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.comparative_dir, 'COMPARATIVE_IGF_vs_Dormancy.png'),
            dpi=300, bbox_inches='tight'
        )
        plt.close()
        
        print("   ✓ Saved: IGF vs dormancy scatter")
    
    def generate_summary_report(self):
        """Generate comprehensive summary report."""
        print("\n" + "="*70)
        print("📊 COMPARATIVE ANALYSIS SUMMARY")
        print("="*70)
        
        # Load results
        self.load_cancer_results()
        
        # Generate comparisons
        igf_comp = self.generate_igf_comparison_table()
        dorm_comp = self.generate_dormancy_comparison_table()
        risk_comp = self.generate_risk_comparison_table()
        
        # Generate visualizations
        print("\n4️⃣ GENERATING COMPARATIVE FIGURES")
        print("-" * 60)
        
        self.generate_igf_pathway_comparison_plot(igf_comp)
        self.generate_dormancy_comparison_plot(dorm_comp)
        self.generate_risk_comparison_plot(risk_comp)
        self.generate_igf_vs_dormancy_scatter()
        
        # Key findings
        print("\n5️⃣ KEY FINDINGS")
        print("-" * 60)
        
        # Find cancer with highest IGF
        if not igf_comp.empty:
            max_igf_idx = igf_comp['Overall_IGF_Activation'].idxmax()
            max_igf_cancer = igf_comp.loc[max_igf_idx, 'Cancer_Type']
            max_igf_value = igf_comp.loc[max_igf_idx, 'Overall_IGF_Activation']
            print(f"\n   🔴 HIGHEST IGF EXPRESSION:")
            print(f"      {max_igf_cancer:<25} IGF Activation: {max_igf_value:.3f}")
        
        # Find cancer with highest dormancy
        if not dorm_comp.empty:
            max_dorm_idx = dorm_comp['Dormant_Pct'].idxmax()
            max_dorm_cancer = dorm_comp.loc[max_dorm_idx, 'Cancer_Type']
            max_dorm_pct = dorm_comp.loc[max_dorm_idx, 'Dormant_Pct']
            print(f"\n   💤 HIGHEST DORMANCY:")
            print(f"      {max_dorm_cancer:<25} Dormant Cells: {max_dorm_pct:.1f}%")
        
        # Find cancer with highest relapse risk (24mo)
        if not risk_comp.empty:
            max_risk_idx = risk_comp['Risk_24mo_Mean'].idxmax()
            max_risk_cancer = risk_comp.loc[max_risk_idx, 'Cancer_Type']
            max_risk_value = risk_comp.loc[max_risk_idx, 'Risk_24mo_Mean']
            print(f"\n   ⚠️  HIGHEST 24-MONTH RELAPSE RISK:")
            print(f"      {max_risk_cancer:<25} Risk: {max_risk_value:.1%}")
        
        print(f"\n✅ COMPARATIVE ANALYSIS COMPLETE")
        print(f"   📁 Results saved to: {self.comparative_dir}/")


if __name__ == "__main__":
    analyzer = ComparativeAnalyzer()
    analyzer.generate_summary_report()
