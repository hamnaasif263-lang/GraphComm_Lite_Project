"""
Lung Cancer Correlation Analysis
Analysis of IGF1R-MYC-CCND-p27 relationships in real lung cancer data
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

class LungCancerCorrelationAnalysis:
    """Perform correlation analysis on lung cancer single-cell data"""
    
    def __init__(self, output_base="output"):
        self.output_base = Path(output_base)
        self.output_base.mkdir(parents=True, exist_ok=True)
        self.dpi = 150
    
    def load_lung_cancer_data(self):
        """Load lung cancer data"""
        lc_path = self.output_base / "lung_cancer/01_preprocessing/lung_cancer_single_cell_data.csv"
        print(f"[LOADING] {lc_path.name}")
        sc_data = pd.read_csv(lc_path)
        print(f"  Shape: {sc_data.shape}")
        return sc_data
    
    def compute_correlations(self, sc_data):
        """Compute all pairwise correlations"""
        
        corr_dict = {}
        
        # All cells
        print("\n[ALL CELLS - n={}]".format(len(sc_data)))
        igf_myc_r, igf_myc_p = pearsonr(sc_data['igf_score'], sc_data['myc_score'])
        igf_prolif_r, igf_prolif_p = pearsonr(sc_data['igf_score'], sc_data['proliferation_score'])
        igf_cdkn1b_r, igf_cdkn1b_p = pearsonr(sc_data['igf_score'], sc_data['cdkn1b_score'])
        myc_prolif_r, myc_prolif_p = pearsonr(sc_data['myc_score'], sc_data['proliferation_score'])
        myc_cdkn1b_r, myc_cdkn1b_p = pearsonr(sc_data['myc_score'], sc_data['cdkn1b_score'])
        prolif_cdkn1b_r, prolif_cdkn1b_p = pearsonr(sc_data['proliferation_score'], sc_data['cdkn1b_score'])
        
        print(f"IGF1R vs MYC:          r = {igf_myc_r:+.4f}  (p = {igf_myc_p:.2e})")
        print(f"IGF1R vs CCND (Prolif):r = {igf_prolif_r:+.4f}  (p = {igf_prolif_p:.2e})")
        print(f"IGF1R vs p27 (CDKN1B): r = {igf_cdkn1b_r:+.4f}  (p = {igf_cdkn1b_p:.2e})")
        print(f"MYC vs CCND:           r = {myc_prolif_r:+.4f}  (p = {myc_prolif_p:.2e})")
        print(f"MYC vs p27:            r = {myc_cdkn1b_r:+.4f}  (p = {myc_cdkn1b_p:.2e})")
        print(f"CCND vs p27:           r = {prolif_cdkn1b_r:+.4f}  (p = {prolif_cdkn1b_p:.2e})")
        
        corr_dict['all'] = {
            'igf_myc': (igf_myc_r, igf_myc_p),
            'igf_prolif': (igf_prolif_r, igf_prolif_p),
            'igf_cdkn1b': (igf_cdkn1b_r, igf_cdkn1b_p),
            'myc_prolif': (myc_prolif_r, myc_prolif_p),
            'myc_cdkn1b': (myc_cdkn1b_r, myc_cdkn1b_p),
            'prolif_cdkn1b': (prolif_cdkn1b_r, prolif_cdkn1b_p),
        }
        
        # By cell type
        print("\n[BY CELL TYPE]")
        for cell_type in sorted(sc_data['cell_type'].unique()):
            ct_data = sc_data[sc_data['cell_type'] == cell_type]
            print(f"\n{cell_type.upper()} (n={len(ct_data):,}):")
            
            try:
                igf_myc_r, _ = pearsonr(ct_data['igf_score'], ct_data['myc_score'])
                igf_prolif_r, _ = pearsonr(ct_data['igf_score'], ct_data['proliferation_score'])
                igf_cdkn1b_r, _ = pearsonr(ct_data['igf_score'], ct_data['cdkn1b_score'])
                myc_prolif_r, _ = pearsonr(ct_data['myc_score'], ct_data['proliferation_score'])
                myc_cdkn1b_r, _ = pearsonr(ct_data['myc_score'], ct_data['cdkn1b_score'])
                prolif_cdkn1b_r, _ = pearsonr(ct_data['proliferation_score'], ct_data['cdkn1b_score'])
                
                print(f"  IGF→MYC:   r = {igf_myc_r:+.4f}")
                print(f"  IGF→CCND:  r = {igf_prolif_r:+.4f}")
                print(f"  IGF→p27:   r = {igf_cdkn1b_r:+.4f}")
                print(f"  MYC↔CCND:  r = {myc_prolif_r:+.4f}")
                print(f"  MYC↔p27:   r = {myc_cdkn1b_r:+.4f}")
                print(f"  CCND↔p27:  r = {prolif_cdkn1b_r:+.4f}")
                
                corr_dict[cell_type] = {
                    'igf_myc': igf_myc_r,
                    'igf_prolif': igf_prolif_r,
                    'igf_cdkn1b': igf_cdkn1b_r,
                    'myc_prolif': myc_prolif_r,
                    'myc_cdkn1b': myc_cdkn1b_r,
                    'prolif_cdkn1b': prolif_cdkn1b_r,
                }
            except:
                print(f"  Error computing correlations for {cell_type}")
        
        return corr_dict
    
    def create_correlation_heatmap(self, sc_data):
        """Create 5x5 correlation heatmap"""
        
        variables = ['igf_score', 'myc_score', 'proliferation_score', 'cdkn1b_score', 'dormancy_score']
        corr_matrix = sc_data[variables].corr()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create heatmap
        sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                    cbar_kws={'label': 'Pearson r'}, vmin=-1, vmax=1, ax=ax,
                    square=True, linewidths=0.5)
        
        # Labels
        labels = ['IGF1R', 'MYC', 'CCND', 'p27', 'Dormancy']
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_yticklabels(labels, rotation=0)
        
        ax.set_title('Lung Cancer: Correlation Matrix (n=300,000 cells)', 
                     fontsize=14, fontweight='bold', pad=20)
        
        plt.tight_layout()
        fig.savefig(self.output_base / 'lc_correlation_heatmap.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white', format='png')
        print(f"  [SAVED] lc_correlation_heatmap.png ({fig.get_size_inches()[0]:.1f}x{fig.get_size_inches()[1]:.1f}\")")
        plt.close()
        
        return corr_matrix
    
    def create_scatter_plots(self, sc_data):
        """Create 2x3 scatter plots with correlations"""
        
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        fig.suptitle('Lung Cancer: Phenotype Correlation Analysis (n=300,000 cells)',
                     fontsize=16, fontweight='bold', y=0.995)
        
        pairs = [
            ('igf_score', 'myc_score', 'IGF1R', 'MYC', 0, 0),
            ('igf_score', 'proliferation_score', 'IGF1R', 'CCND', 0, 1),
            ('myc_score', 'proliferation_score', 'MYC', 'CCND', 0, 2),
            ('myc_score', 'cdkn1b_score', 'MYC', 'p27', 1, 0),
            ('proliferation_score', 'cdkn1b_score', 'CCND', 'p27', 1, 1),
            ('igf_score', 'cdkn1b_score', 'IGF1R', 'p27', 1, 2),
        ]
        
        for x_col, y_col, x_label, y_label, i, j in pairs:
            ax = axes[i, j]
            
            # Hexbin density
            try:
                hb = ax.hexbin(sc_data[x_col], sc_data[y_col], gridsize=25,
                              cmap='YlOrRd', mincnt=1, alpha=0.75, edgecolors='none')
                
                # Trend line
                z = np.polyfit(sc_data[x_col], sc_data[y_col], 1)
                p = np.poly1d(z)
                x_trend = np.linspace(sc_data[x_col].min(), sc_data[x_col].max(), 100)
                ax.plot(x_trend, p(x_trend), 'r-', linewidth=2.5, alpha=0.8, label='Trend')
                
                # Correlation
                r, p_val = pearsonr(sc_data[x_col], sc_data[y_col])
                
                # Labels
                ax.set_xlabel(x_label, fontsize=11, fontweight='bold')
                ax.set_ylabel(y_label, fontsize=11, fontweight='bold')
                ax.set_title(f'r = {r:+.4f} (p < {p_val:.2e})',
                           fontsize=10, fontweight='bold', color='darkred')
                
                # Colorbar
                cb = plt.colorbar(hb, ax=ax)
                cb.set_label('Count', fontsize=9)
                
                ax.grid(True, alpha=0.3)
                ax.legend(loc='upper left', fontsize=9)
                
            except Exception as e:
                print(f"  [ERROR {x_label} vs {y_label}] {e}")
        
        plt.tight_layout()
        try:
            fig.savefig(self.output_base / 'lc_scatter_correlations.png',
                       dpi=150, bbox_inches='tight', facecolor='white', format='png')
            print(f"  [SAVED] lc_scatter_correlations.png")
        except Exception as e:
            print(f"  [ERROR saving scatter] {e}")
        plt.close()
    
    def save_correlation_matrix(self, sc_data):
        """Save correlation matrix as CSV"""
        
        variables = ['igf_score', 'myc_score', 'proliferation_score', 'cdkn1b_score', 'dormancy_score']
        corr_matrix = sc_data[variables].corr()
        
        output_path = self.output_base / 'lc_correlation_matrix.csv'
        corr_matrix.to_csv(output_path)
        print(f"  [SAVED] lc_correlation_matrix.csv")
        
        return corr_matrix
    
    def run_analysis(self):
        """Run complete analysis pipeline"""
        
        print("\n" + "="*80)
        print("LUNG CANCER CORRELATION ANALYSIS - REAL DATA")
        print("="*80)
        
        # Load data
        sc_data = self.load_lung_cancer_data()
        
        # Compute correlations
        corr_results = self.compute_correlations(sc_data)
        
        # Generate visualizations
        print("\n[GENERATING VISUALIZATIONS]")
        self.create_correlation_heatmap(sc_data)
        self.create_scatter_plots(sc_data)
        self.save_correlation_matrix(sc_data)
        
        print("\nAnalysis complete")


if __name__ == "__main__":
    analyzer = LungCancerCorrelationAnalysis()
    analyzer.run_analysis()
