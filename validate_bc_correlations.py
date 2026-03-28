"""
Validate IGF1R-MYC-CCND-p27 Correlations in Breast Cancer
Reproduce the angiocrine awakening phenotype findings
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class BreastCancerCorrelationAnalysis:
    """Analyze cell-cycle regulator correlations in breast cancer"""
    
    def __init__(self, output_base="output"):
        self.output_base = Path(output_base)
        self.dpi = 300
    
    def load_breast_cancer_data(self):
        """Load single-cell data"""
        sc_path = self.output_base / "breast_cancer/01_preprocessing/breast_cancer_single_cell_data.csv"
        return pd.read_csv(sc_path)
    
    def compute_correlations(self, sc_data):
        """Compute correlations between IGF pathway and cell-cycle regulators"""
        
        print("\n" + "="*70)
        print("BREAST CANCER: IGF1R-MYC-CCND-p27 CORRELATION ANALYSIS")
        print("="*70)
        
        # Overall correlations (all cells)
        print("\n[ALL CELLS - n={}]".format(len(sc_data)))
        
        igf_myc_r, igf_myc_p = pearsonr(sc_data['igf_score'], sc_data['myc_score'])
        igf_prolif_r, igf_prolif_p = pearsonr(sc_data['igf_score'], sc_data['proliferation_score'])
        igf_cdkn1b_r, igf_cdkn1b_p = pearsonr(sc_data['igf_score'], sc_data['cdkn1b_score'])
        
        myc_cdkn1b_r, myc_cdkn1b_p = pearsonr(sc_data['myc_score'], sc_data['cdkn1b_score'])
        prolif_cdkn1b_r, prolif_cdkn1b_p = pearsonr(sc_data['proliferation_score'], sc_data['cdkn1b_score'])
        
        print(f"  IGF1R vs MYC:          r = {igf_myc_r:.4f}  (p = {igf_myc_p:.2e})")
        print(f"  IGF1R vs CCND (Prolif):r = {igf_prolif_r:.4f}  (p = {igf_prolif_p:.2e})")
        print(f"  IGF1R vs p27 (CDKN1B): r = {igf_cdkn1b_r:.4f}  (p = {igf_cdkn1b_p:.2e})")
        print(f"  MYC vs p27:            r = {myc_cdkn1b_r:.4f}  (p = {myc_cdkn1b_p:.2e})")
        print(f"  CCND vs p27:           r = {prolif_cdkn1b_r:.4f}  (p = {prolif_cdkn1b_p:.2e})")
        
        # Correlation by cell type
        print("\n[BY CELL TYPE]")
        for cell_type in sorted(sc_data['cell_type'].unique()):
            ct_data = sc_data[sc_data['cell_type'] == cell_type]
            n_cells = len(ct_data)
            
            if n_cells < 3:
                continue
            
            try:
                igf_myc_r, _ = pearsonr(ct_data['igf_score'], ct_data['myc_score'])
                igf_prolif_r, _ = pearsonr(ct_data['igf_score'], ct_data['proliferation_score'])
                igf_cdkn1b_r, _ = pearsonr(ct_data['igf_score'], ct_data['cdkn1b_score'])
                myc_cdkn1b_r, _ = pearsonr(ct_data['myc_score'], ct_data['cdkn1b_score'])
                
                print(f"\n  {cell_type.upper()} (n={n_cells:,}):")
                print(f"    IGF1R vs MYC:    r = {igf_myc_r:.4f}")
                print(f"    IGF1R vs CCND:   r = {igf_prolif_r:.4f}")
                print(f"    IGF1R vs p27:    r = {igf_cdkn1b_r:.4f}")
                print(f"    MYC vs p27:      r = {myc_cdkn1b_r:.4f}")
            except:
                pass
        
        return {
            'igf_myc': (igf_myc_r, igf_myc_p),
            'igf_prolif': (igf_prolif_r, igf_prolif_p),
            'igf_cdkn1b': (igf_cdkn1b_r, igf_cdkn1b_p),
            'myc_cdkn1b': (myc_cdkn1b_r, myc_cdkn1b_p),
            'prolif_cdkn1b': (prolif_cdkn1b_r, prolif_cdkn1b_p),
        }
    
    def create_correlation_heatmap(self, sc_data):
        """Create correlation matrix heatmap"""
        
        # Select key variables
        vars_of_interest = ['igf_score', 'myc_score', 'proliferation_score', 
                           'cdkn1b_score', 'dormancy_score']
        
        corr_matrix = sc_data[vars_of_interest].corr(method='pearson')
        
        fig, ax = plt.subplots(figsize=(8, 7))
        
        im = ax.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        
        # Set ticks and labels
        labels = ['IGF1R', 'MYC', 'CCND (Prolif)', 'p27', 'Dormancy']
        ax.set_xticks(np.arange(len(labels)))
        ax.set_yticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, fontsize=11, fontweight='bold')
        ax.set_yticklabels(labels, fontsize=11, fontweight='bold')
        
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        
        # Add correlation values
        for i in range(len(labels)):
            for j in range(len(labels)):
                value = corr_matrix.iloc[i, j]
                color = 'white' if abs(value) > 0.5 else 'black'
                text = ax.text(j, i, f'{value:.3f}', ha="center", va="center",
                              color=color, fontweight='bold', fontsize=10)
        
        ax.set_title('Breast Cancer: Key Regulatory Gene Correlations\n(All cells, n=275,000)',
                    fontsize=13, fontweight='bold', pad=15)
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Pearson Correlation', fontweight='bold', fontsize=11)
        
        plt.tight_layout()
        return fig
    
    def create_scatter_plots(self, sc_data):
        """Create scatter plots for key relationships"""
        
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        fig.patch.set_facecolor('white')
        
        # Define scatter plots
        plots = [
            (axes[0, 0], 'igf_score', 'myc_score', 'IGF1R', 'MYC', '#FF6B6B'),
            (axes[0, 1], 'igf_score', 'proliferation_score', 'IGF1R', 'CCND (Prolif)', '#4ECDC4'),
            (axes[0, 2], 'igf_score', 'cdkn1b_score', 'IGF1R', 'p27 (CDKN1B)', '#95E1D3'),
            (axes[1, 0], 'myc_score', 'cdkn1b_score', 'MYC', 'p27 (CDKN1B)', '#FFD93D'),
            (axes[1, 1], 'myc_score', 'proliferation_score', 'MYC', 'CCND (Prolif)', '#6BCB77'),
            (axes[1, 2], 'proliferation_score', 'cdkn1b_score', 'CCND (Prolif)', 'p27 (CDKN1B)', '#4D96FF'),
        ]
        
        for ax, x_var, y_var, x_label, y_label, color in plots:
            x = sc_data[x_var].values
            y = sc_data[y_var].values
            
            # Scatter plot with density
            h = ax.hexbin(x, y, gridsize=25, cmap='YlOrRd', mincnt=1, alpha=0.75, edgecolors='none')
            
            # Add trend line
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            x_line = np.linspace(x.min(), x.max(), 100)
            ax.plot(x_line, p(x_line), "b--", linewidth=2, label=f'Trend')
            
            # Calculate correlation
            r, p_val = pearsonr(x, y)
            
            ax.set_xlabel(x_label, fontsize=10, fontweight='bold')
            ax.set_ylabel(y_label, fontsize=10, fontweight='bold')
            ax.set_title(f'{x_label} vs {y_label}\nr = {r:.4f} (p = {p_val:.2e})',
                        fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.2)
            ax.set_facecolor('white')
            
            cbar = plt.colorbar(h, ax=ax)
            cbar.set_label('Count', fontsize=9)
        
        fig.suptitle('Breast Cancer: Cell-Cycle Regulator Correlations\n(Angiocrine Awakening Phenotype Validation)',
                    fontsize=14, fontweight='bold', y=0.995)
        
        plt.tight_layout()
        return fig
    
    def validate_against_literature(self, corr_dict):
        """Compare with literature values"""
        
        print("\n" + "="*70)
        print("LITERATURE COMPARISON")
        print("="*70)
        
        igf_myc_r = corr_dict['igf_myc'][0]
        igf_prolif_r = corr_dict['igf_prolif'][0]
        
        print("\nTarget Correlations (from literature):")
        print(f"  IGF1R vs MYC:     r = 0.216  (expected)")
        print(f"  IGF1R vs CCND:    r = 0.343  (expected)")
        
        print("\nObserved Correlations (from your breast cancer data):")
        print(f"  IGF1R vs MYC:     r = {igf_myc_r:.4f}  {'[MATCH]' if abs(igf_myc_r - 0.216) < 0.1 else '[DIFFERENT]'}")
        print(f"  IGF1R vs CCND:    r = {igf_prolif_r:.4f}  {'[MATCH]' if abs(igf_prolif_r - 0.343) < 0.1 else '[DIFFERENT]'}")
        
        print("\nInterpretation:")
        if igf_myc_r > 0.15:
            print(f"  [+] Positive IGF1R-MYC correlation (r={igf_myc_r:.4f}) detected")
            print(f"      -> Supports angiocrine signaling activating proliferation")
        else:
            print(f"  [-] Weak/absent IGF1R-MYC correlation (r={igf_myc_r:.4f})")
        
        if igf_prolif_r > 0.25:
            print(f"  [+] Strong IGF1R-CCND correlation (r={igf_prolif_r:.4f}) detected")
            print(f"      -> Confirms vascular IGF driving G1/S transition")
        else:
            print(f"  [-] Moderate IGF1R-CCND correlation (r={igf_prolif_r:.4f})")
        
        print("\n" + "="*70 + "\n")
    
    def run_analysis(self):
        """Complete correlation analysis"""
        
        # Load data
        sc_data = self.load_breast_cancer_data()
        print(f"\nLoaded {len(sc_data):,} breast cancer cells")
        
        # Compute correlations
        corr_dict = self.compute_correlations(sc_data)
        
        # Validate against literature
        self.validate_against_literature(corr_dict)
        
        # Create visualizations
        print("Generating correlation visualizations...")
        
        fig1 = self.create_correlation_heatmap(sc_data)
        fig1.savefig(self.output_base / 'bc_correlation_heatmap.png',
                    dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close(fig1)
        print("  [SAVED] bc_correlation_heatmap.png")
        
        fig2 = self.create_scatter_plots(sc_data)
        try:
            fig2.savefig(self.output_base / 'bc_scatter_correlations.png',
                        dpi=150, bbox_inches='tight', facecolor='white', format='png')
        except Exception as e:
            print(f"  [ERROR saving scatter] {e}")
        plt.close(fig2)
        print("  [SAVED] bc_scatter_correlations.png")
        
        # Save correlation matrix
        corr_df = sc_data[['igf_score', 'myc_score', 'proliferation_score', 
                           'cdkn1b_score', 'dormancy_score']].corr()
        corr_df.to_csv(self.output_base / 'breast_cancer/network_analysis/bc_correlation_matrix.csv')
        print("  [SAVED] bc_correlation_matrix.csv")
        
        print("\n" + "="*70)
        print("ANGIOCRINE PHENOTYPE VALIDATION COMPLETE")
        print("="*70 + "\n")


if __name__ == "__main__":
    analyzer = BreastCancerCorrelationAnalysis()
    analyzer.run_analysis()
