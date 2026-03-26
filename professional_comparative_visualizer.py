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
        plt.style.use('seaborn-v0_8-darkgrid')
        sns.set_palette("husl")
        self.dpi = 300
        self.font_family = 'Arial'
        self.font_sizes = {
            'title': 16,
            'label': 12,
            'tick': 10,
            'legend': 10
        }
        
    def load_all_cancer_data(self):
        """Load IGF and dormancy data from individual cancer analyses"""
        cancer_types = ['breast_cancer', 'lung_cancer', 'prostate_cancer']
        igf_data = {}
        dormancy_data = {}
        risk_data = {}
        
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
            
            # Load risk scores
            risk_file = cancer_dir / "02_window_of_risk" / f"{cancer}_risk_scores.csv"
            if risk_file.exists():
                risk_data[cancer] = pd.read_csv(risk_file)
        
        return igf_data, dormancy_data, risk_data
    
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
                    fontsize=self.font_sizes['title'],
                    fontweight='bold',
                    pad=20)
        ax.set_xlabel('Cancer Type', fontsize=self.font_sizes['label'], fontweight='bold')
        ax.set_ylabel('Pathway Component', fontsize=self.font_sizes['label'], fontweight='bold')
        
        plt.tight_layout()
        fig.savefig(self.comparative_dir / 'IGF_Pathway_Comparison_Heatmap.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print("✓ IGF Pathway Comparison Heatmap generated")
    
    def create_igf_activation_comparison(self, igf_data):
        """Create side-by-side IGF activation comparison bar chart"""
        cancer_names = ['Breast', 'Lung', 'Prostate']
        
        # Calculate overall IGF activation (mean of main components)
        igf_activations = []
        for cancer in sorted(igf_data.keys()):
            data = igf_data[cancer]
            main_components = ['IGF_Ligand_Mean', 'IGF_Receptor_Mean', 'PI3K_AKT_Mean', 'mTOR_Mean']
            values = [data[comp].values[0] for comp in main_components if comp in data.columns]
            activation = np.mean(values) if values else 0.0
            igf_activations.append(max(0.001, activation))  # Ensure no zero values
        
        fig, ax = plt.subplots(figsize=(10, 7))
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
        bars = ax.bar(cancer_names, igf_activations, color=colors, alpha=0.8, edgecolor='black', linewidth=2)
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}',
                   ha='center', va='bottom', fontsize=11, fontweight='bold')
        
        ax.set_ylabel('IGF Pathway Activation Score', fontsize=self.font_sizes['label'], fontweight='bold')
        ax.set_xlabel('Cancer Type', fontsize=self.font_sizes['label'], fontweight='bold')
        ax.set_title('Overall IGF Pathway Activation by Cancer Type',
                    fontsize=self.font_sizes['title'],
                    fontweight='bold',
                    pad=20)
        max_val = max(igf_activations)
        ax.set_ylim(0, max_val * 1.15)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        
        plt.tight_layout()
        fig.savefig(self.comparative_dir / 'IGF_Activation_Comparison.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print("✓ IGF Activation Comparison generated")
    
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
        
        bars1 = ax.bar(x - width, dormancy_percentages, width, label='Dormant', color=colors[0], alpha=0.8, edgecolor='black', linewidth=1.5)
        bars2 = ax.bar(x, transitional_percentages, width, label='Transitional', color=colors[1], alpha=0.8, edgecolor='black', linewidth=1.5)
        bars3 = ax.bar(x + width, proliferative_percentages, width, label='Proliferative', color=colors[2], alpha=0.8, edgecolor='black', linewidth=1.5)
        
        ax.set_ylabel('Percentage of Cells (%)', fontsize=self.font_sizes['label'], fontweight='bold')
        ax.set_xlabel('Cancer Type', fontsize=self.font_sizes['label'], fontweight='bold')
        ax.set_title('Cell Population Composition by Cancer Type\n(Dormancy Signatures)',
                    fontsize=self.font_sizes['title'],
                    fontweight='bold',
                    pad=20)
        ax.set_xticks(x)
        ax.set_xticklabels(cancer_names)
        ax.legend(fontsize=self.font_sizes['legend'], loc='upper right', framealpha=0.9)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        
        plt.tight_layout()
        fig.savefig(self.comparative_dir / 'Dormancy_Signature_Comparison.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print("✓ Dormancy Signature Comparison generated")
    
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
            activation = np.mean([igf_df[comp].values[0] for comp in main_components if comp in igf_df.columns])
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
        z = np.polyfit(igf_activations, dormancy_percentages, 1)
        p = np.poly1d(z)
        x_trend = np.linspace(min(igf_activations)-0.05, max(igf_activations)+0.05, 100)
        ax.plot(x_trend, p(x_trend), "k--", alpha=0.5, linewidth=2, label='Trend')
        
        ax.set_xlabel('IGF Pathway Activation Score', fontsize=self.font_sizes['label'], fontweight='bold')
        ax.set_ylabel('Dormancy Percentage (%)', fontsize=self.font_sizes['label'], fontweight='bold')
        ax.set_title('Relationship Between IGF Activation and Dormancy',
                    fontsize=self.font_sizes['title'],
                    fontweight='bold',
                    pad=20)
        ax.legend(fontsize=self.font_sizes['legend'], loc='best', framealpha=0.9)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        
        plt.tight_layout()
        fig.savefig(self.comparative_dir / 'IGF_Dormancy_Correlation.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print("✓ IGF-Dormancy Correlation plot generated")
    
    def create_pathway_component_ranking(self, igf_data):
        """Create ranking of pathway components across cancer types"""
        cancer_names = ['Breast', 'Lung', 'Prostate']
        components = ['IGF_Ligand_Mean', 'IGF_Receptor_Mean', 'PI3K_AKT_Mean',
                     'mTOR_Mean', 'MAPK_Mean', 'MYC_Mean']
        
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        axes = axes.flatten()
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
        
        for idx, component in enumerate(components):
            ax = axes[idx]
            values = []
            
            for cancer in sorted(igf_data.keys()):
                data = igf_data[cancer]
                if component in data.columns:
                    values.append(data[component].values[0])
                else:
                    values.append(0)
            
            bars = ax.bar(cancer_names, values, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
            
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{height:.3f}',
                       ha='center', va='bottom', fontsize=9, fontweight='bold')
            
            ax.set_title(component.replace('_', ' '), fontsize=11, fontweight='bold')
            ax.set_ylabel('Expression Level', fontsize=10)
            ax.set_ylim(0, max(values) * 1.2 if values else 1)
            ax.grid(axis='y', alpha=0.3, linestyle='--')
            ax.set_axisbelow(True)
        
        fig.suptitle('Pathway Component Expression Ranking Across Cancer Types',
                    fontsize=self.font_sizes['title'],
                    fontweight='bold',
                    y=0.995)
        plt.tight_layout()
        fig.savefig(self.comparative_dir / 'Pathway_Components_Ranking.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print("✓ Pathway Components Ranking generated")
    
    def create_summary_statistics_table(self, igf_data, dormancy_data):
        """Generate professional statistics table as figure"""
        cancer_names = list(sorted(igf_data.keys()))
        stats_data = []
        
        for cancer in cancer_names:
            igf_df = igf_data[cancer]
            dorm_df = dormancy_data[cancer]
            
            dormancy_pct = (dorm_df['cell_type'] == 'Dormant').sum() / len(dorm_df) * 100
            
            main_components = ['IGF_Ligand_Mean', 'IGF_Receptor_Mean', 'PI3K_AKT_Mean', 'mTOR_Mean']
            igf_activation = np.mean([igf_df[comp].values[0] for comp in main_components if comp in igf_df.columns])
            
            stats_data.append({
                'Cancer Type': cancer.replace('_', ' ').title(),
                'IGF Activation': f'{igf_activation:.3f}',
                'Dormant Cells': f'{dormancy_pct:.1f}%',
                'mTOR Level': f'{igf_df["mTOR_Mean"].values[0]:.3f}',
                'MAPK Level': f'{igf_df["MAPK_Mean"].values[0]:.3f}',
                'MYC Level': f'{igf_df["MYC_Mean"].values[0]:.3f}'
            })
        
        stats_df = pd.DataFrame(stats_data)
        
        fig, ax = plt.subplots(figsize=(14, 4))
        ax.axis('tight')
        ax.axis('off')
        
        table = ax.table(cellText=stats_df.values,
                        colLabels=stats_df.columns,
                        cellLoc='center',
                        loc='center',
                        colWidths=[0.18, 0.13, 0.15, 0.13, 0.13, 0.13])
        
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1, 2.5)
        
        # Style header
        for i in range(len(stats_df.columns)):
            table[(0, i)].set_facecolor('#3498DB')
            table[(0, i)].set_text_props(weight='bold', color='black')
        
        # Alternate row colors
        for i in range(1, len(stats_df) + 1):
            for j in range(len(stats_df.columns)):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#ECF0F1')
                else:
                    table[(i, j)].set_facecolor('#FFFFFF')
        
        plt.title('Comparative Cancer Analysis Summary Statistics',
                 fontsize=self.font_sizes['title'],
                 fontweight='bold',
                 pad=20)
        
        plt.tight_layout()
        fig.savefig(self.comparative_dir / 'Summary_Statistics_Table.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print("✓ Summary Statistics Table generated")
    
    def generate_all_visualizations(self):
        """Generate all professional comparative analysis visualizations"""
        print("\n" + "="*60)
        print("GENERATING PROFESSIONAL COMPARATIVE ANALYSIS VISUALIZATIONS")
        print("="*60)
        
        igf_data, dormancy_data, risk_data = self.load_all_cancer_data()
        
        if not igf_data:
            print("❌ No IGF data found. Ensure individual cancer analyses are complete.")
            return
        
        print("\nGenerating publication-quality figures (300 DPI)...")
        
        self.create_igf_comparison_heatmap(igf_data)
        self.create_igf_activation_comparison(igf_data)
        self.create_dormancy_comparison(dormancy_data)
        self.create_igf_dormancy_correlation(igf_data, dormancy_data)
        self.create_pathway_component_ranking(igf_data)
        self.create_summary_statistics_table(igf_data, dormancy_data)
        
        print("\n" + "="*60)
        print("✓ ALL PROFESSIONAL VISUALIZATIONS COMPLETED")
        print(f"✓ Output location: {self.comparative_dir}")
        print("="*60 + "\n")


if __name__ == "__main__":
    visualizer = ProfessionalComparativeVisualizer()
    visualizer.generate_all_visualizations()
