"""
Professional Visualization Suite for Comparative Cancer Analysis
Generates publication-quality figures in Nature-style format (300 DPI)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class ProfessionalComparativeVisualizer:
    """Generates professional publication-quality graphs for cancer comparative analysis"""
    
    def __init__(self, output_base="output"):
        self.output_base = Path(output_base)
        self.comparative_dir = self.output_base / "04_COMPARATIVE_ANALYSIS"
        self.comparative_dir.mkdir(parents=True, exist_ok=True)
        
        # Publication-grade styling
        plt.style.use('seaborn-v0_8-whitegrid')
        sns.set_palette("husl")
        self.dpi = 300
        
    def load_all_cancer_data(self):
        """Load IGF and dormancy data from individual cancer analyses"""
        cancer_types = ['breast_cancer', 'lung_cancer', 'prostate_cancer']
        igf_data = {}
        dormancy_data = {}
        
        for cancer in cancer_types:
            cancer_dir = self.output_base / cancer
            
            # Load IGF summary
            igf_file = cancer_dir / "01_igf_pathway_analysis" / f"{cancer}_igf_summary.csv"
            if igf_file.exists():
                igf_data[cancer] = pd.read_csv(igf_file)
            
            # Load dormancy signature
            dorm_file = cancer_dir / "01_dormancy_signatures" / f"{cancer}_dormancy_signature.csv"
            if dorm_file.exists():
                dormancy_data[cancer] = pd.read_csv(dorm_file)
        
        return igf_data, dormancy_data
    
    def create_igf_comparison_heatmap(self, igf_data):
        """Create publication-quality IGF pathway comparison heatmap"""
        cancer_names = ['Breast', 'Lung', 'Prostate']
        
        # Extract pathway components for each cancer
        pathway_components = [
            'IGF_Ligand_Mean', 'IGF_Receptor_Mean', 'PI3K_AKT_Mean',
            'mTOR_Mean', 'MAPK_Mean', 'MYC_Mean', 
            'Dormancy_Mean', 'Proliferation_Mean'
        ]
        
        comparison_matrix = []
        for cancer in sorted(igf_data.keys()):
            row = []
            data = igf_data[cancer]
            for component in pathway_components:
                if component in data.columns:
                    row.append(data[component].values[0])
                else:
                    row.append(0)
            comparison_matrix.append(row)
        
        comparison_matrix = np.array(comparison_matrix)
        
        # Create professional heatmap
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # Normalize for better color visualization
        sns.heatmap(comparison_matrix.T, 
                   xticklabels=cancer_names,
                   yticklabels=[c.replace('_', ' ') for c in pathway_components],
                   cmap='RdYlGn',
                   annot=True,
                   fmt='.3f',
                   cbar_kws={'label': 'Expression Level'},
                   linewidths=0.5,
                   linecolor='gray',
                   ax=ax,
                   vmin=0,
                   vmax=1)
        
        ax.set_title('IGF Pathway Component Comparison Across Cancer Types',
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Cancer Type', fontsize=12, fontweight='bold')
        ax.set_ylabel('Pathway Component', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        fig.savefig(self.comparative_dir / 'IGF_Pathway_Comparison_Heatmap.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print("[OK] IGF Pathway Comparison Heatmap")
    
    def create_igf_activation_comparison(self, igf_data):
        """Create side-by-side IGF activation comparison bar chart"""
        cancer_names = ['Breast', 'Lung', 'Prostate']
        
        # Calculate overall IGF activation
        igf_activations = []
        for cancer in sorted(igf_data.keys()):
            data = igf_data[cancer]
            main_components = ['IGF_Ligand_Mean', 'IGF_Receptor_Mean', 'PI3K_AKT_Mean', 'mTOR_Mean']
            values = [data[comp].values[0] for comp in main_components if comp in data.columns]
            activation = np.mean(values) if values else 0.0
            igf_activations.append(max(0.001, activation))
        
        fig, ax = plt.subplots(figsize=(10, 7))
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
        bars = ax.bar(cancer_names, igf_activations, color=colors, alpha=0.8, 
                     edgecolor='black', linewidth=2)
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}', ha='center', va='bottom', 
                   fontsize=11, fontweight='bold')
        
        ax.set_ylabel('IGF Pathway Activation Score', fontsize=12, fontweight='bold')
        ax.set_xlabel('Cancer Type', fontsize=12, fontweight='bold')
        ax.set_title('Overall IGF Pathway Activation',
                    fontsize=14, fontweight='bold', pad=20)
        max_val = max(igf_activations)
        ax.set_ylim(0, max_val * 1.15)
        ax.grid(axis='y', alpha=0.3)
        ax.set_axisbelow(True)
        
        plt.tight_layout()
        fig.savefig(self.comparative_dir / 'IGF_Activation_Comparison.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print("[OK] IGF Activation Comparison")
    
    def create_dormancy_comparison(self, dormancy_data):
        """Create professional dormancy signature comparison"""
        cancer_names = ['Breast', 'Lung', 'Prostate']
        dormancy_percentages = []
        transitional_percentages = []
        proliferative_percentages = []
        
        for cancer in sorted(dormancy_data.keys()):
            data = dormancy_data[cancer]
            dormancy_pct = (data['cell_type'] == 'Dormant').sum() / len(data) * 100
            transitional_pct = (data['cell_type'] == 'Transitional').sum() / len(data) * 100
            proliferative_pct = (data['cell_type'] == 'Proliferative').sum() / len(data) * 100
            
            dormancy_percentages.append(dormancy_pct)
            transitional_percentages.append(transitional_pct)
            proliferative_percentages.append(proliferative_pct)
        
        fig, ax = plt.subplots(figsize=(12, 7))
        
        x = np.arange(len(cancer_names))
        width = 0.25
        
        colors = ['#E74C3C', '#F39C12', '#27AE60']
        
        ax.bar(x - width, dormancy_percentages, width, label='Dormant', 
              color=colors[0], alpha=0.8, edgecolor='black', linewidth=1.5)
        ax.bar(x, transitional_percentages, width, label='Transitional', 
              color=colors[1], alpha=0.8, edgecolor='black', linewidth=1.5)
        ax.bar(x + width, proliferative_percentages, width, label='Proliferative', 
              color=colors[2], alpha=0.8, edgecolor='black', linewidth=1.5)
        
        ax.set_ylabel('Percentage of Cells', fontsize=12, fontweight='bold')
        ax.set_xlabel('Cancer Type', fontsize=12, fontweight='bold')
        ax.set_title('Cell Population Composition - Dormancy Signatures',
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xticks(x)
        ax.set_xticklabels(cancer_names)
        ax.legend(fontsize=11, loc='upper right')
        ax.grid(axis='y', alpha=0.3)
        ax.set_axisbelow(True)
        
        plt.tight_layout()
        fig.savefig(self.comparative_dir / 'Dormancy_Signature_Comparison.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print("[OK] Dormancy Signature Comparison")
    
    def create_igf_dormancy_correlation(self, igf_data, dormancy_data):
        """Create scatter plot showing IGF activation vs dormancy"""
        cancer_names = ['Breast', 'Lung', 'Prostate']
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
        
        igf_activations = []
        dormancy_percentages = []
        
        for cancer in sorted(igf_data.keys()):
            # IGF activation
            igf_df = igf_data[cancer]
            main_components = ['IGF_Ligand_Mean', 'IGF_Receptor_Mean', 'PI3K_AKT_Mean', 'mTOR_Mean']
            values = [igf_df[comp].values[0] for comp in main_components if comp in igf_df.columns]
            activation = np.mean(values) if values else 0.0
            igf_activations.append(activation)
            
            # Dormancy percentage
            dorm_df = dormancy_data[cancer]
            dormancy_pct = (dorm_df['cell_type'] == 'Dormant').sum() / len(dorm_df) * 100
            dormancy_percentages.append(dormancy_pct)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        for i, cancer_name in enumerate(cancer_names):
            ax.scatter(igf_activations[i], dormancy_percentages[i],
                      s=400, color=colors[i], alpha=0.7, edgecolor='black',
                      linewidth=2, label=cancer_name, zorder=3)
        
        # Add trend line
        if len(igf_activations) > 1:
            z = np.polyfit(igf_activations, dormancy_percentages, 1)
            p = np.poly1d(z)
            x_trend = np.linspace(min(igf_activations)-0.05, max(igf_activations)+0.05, 100)
            ax.plot(x_trend, p(x_trend), "k--", alpha=0.5, linewidth=2, label='Trend')
        
        ax.set_xlabel('IGF Pathway Activation Score', fontsize=12, fontweight='bold')
        ax.set_ylabel('Dormancy Percentage', fontsize=12, fontweight='bold')
        ax.set_title('IGF Activation vs Dormancy Relationship',
                    fontsize=14, fontweight='bold', pad=20)
        ax.legend(fontsize=11, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_axisbelow(True)
        
        plt.tight_layout()
        fig.savefig(self.comparative_dir / 'IGF_Dormancy_Correlation.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print("[OK] IGF-Dormancy Correlation")
    
    def generate_all_visualizations(self):
        """Generate all professional comparative analysis visualizations"""
        print("\n" + "="*60)
        print("GENERATING PROFESSIONAL VISUALIZATIONS")
        print("="*60)
        
        igf_data, dormancy_data = self.load_all_cancer_data()
        
        if not igf_data:
            print("[ERROR] No IGF data found.")
            return
        
        print("\nGenerating publication-quality figures...")
        
        try:
            self.create_igf_comparison_heatmap(igf_data)
            self.create_igf_activation_comparison(igf_data)
            self.create_dormancy_comparison(dormancy_data)
            self.create_igf_dormancy_correlation(igf_data, dormancy_data)
        except Exception as e:
            print(f"[ERROR] {str(e)}")
            return
        
        print("\n" + "="*60)
        print("[SUCCESS] ALL PROFESSIONAL VISUALIZATIONS COMPLETE")
        print(f"[OUTPUT] Location: {self.comparative_dir}")
        print("="*60 + "\n")


if __name__ == "__main__":
    visualizer = ProfessionalComparativeVisualizer()
    visualizer.generate_all_visualizations()
