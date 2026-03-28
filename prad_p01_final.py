#!/usr/bin/env python
"""
PRAD_P01 Final Pipeline - One-shot read with fast processing
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

print("Starting PRAD_P01 Pipeline...", flush=True)

DATA_DIR = Path('data/prostate/raw')
RESULTS_DIR = Path('graphcomm/results/prostate')
PLOTS_DIR = Path('graphcomm/plots/prostate')

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

count_file = DATA_DIR / 'GSM4203181_data.raw.matrix.txt'
if not count_file.exists():
    count_file = DATA_DIR / 'GSM4203181_data.matrix.txt'

print(f"Loading: {count_file}", flush=True)
print(f"File size: {count_file.stat().st_size / 1e9:.2f} GB", flush=True)

try:
    print("Reading entire matrix into memory...", flush=True)
    counts = pd.read_csv(count_file, sep='\t', index_col=0, low_memory=False)
    print(f"Loaded: {counts.shape}", flush=True)
    print(f"Memory: {counts.memory_usage(deep=True).sum() / 1e9:.2f} GB", flush=True)
    
    # Quick convert to float32
    print("Converting to float32...", flush=True)
    counts = counts.astype(np.float32)
    
    # Filter
    print("Filtering...", flush=True)
    genes_per_cell = (counts > 0).sum(axis=0)
    counts = counts.loc[:, genes_per_cell >= 200]
    
    cells_per_gene = (counts > 0).sum(axis=1)
    counts = counts.loc[cells_per_gene >= 3, :]
    print(f"After filtering: {counts.shape}", flush=True)
    
    # Normalize
    print("Normalizing...", flush=True)
    lib_size = counts.sum(axis=0)
    counts_norm = counts.mul(1e4 / lib_size, axis=1)
    counts_log = np.log1p(counts_norm)
    
    # Save
    print("Saving normalized matrix...", flush=True)
    counts_log.to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_norm.csv')
    print(f"✓ Saved to {RESULTS_DIR / 'scRNA_PRAD_P01_norm.csv'}", flush=True)
    
    # Samples
    print("Creating visualization sample...", flush=True)
    n_vis = min(3000, counts_log.shape[1])
    vis_idx = np.random.choice(counts_log.shape[1], n_vis, replace=False)
    counts_log.iloc[:, vis_idx].to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_vis_sample.csv')
    print(f"✓ Saved sample {counts_log.iloc[:, vis_idx].shape}", flush=True)
    
    # Summary
    with open(RESULTS_DIR / 'scRNA_PRAD_P01_summary.txt', 'w') as f:
        f.write(f"""PRAD_P01 PREPROCESSING SUMMARY
Cells: {counts.shape[1]:,}
Genes: {counts.shape[0]:,}
Library Mean: {lib_size.mean():.0f}
Library Median: {lib_size.median():.0f}
""")
    print("✓ Saved summary", flush=True)
    
    # IGF
    print("\nIGF Pathway Analysis...", flush=True)
    IGF_GENES = ["IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1", "MTOR", "FOXO1"]
    igf_found = [g for g in IGF_GENES if g in counts_log.index]
    print(f"IGF genes found: {len(igf_found)}", flush=True)
    
    if igf_found:
        igf_scores = counts_log.loc[igf_found].mean(axis=0)
        
        pd.DataFrame({'cell_barcode': counts_log.columns, 'igf_activity_score': igf_scores.values}).to_csv(
            RESULTS_DIR / 'PRAD_P01_IGF_cell_scores.csv', index=False)
        
        pd.DataFrame({
            'metric': ['n_genes', 'n_cells', 'mean', 'median', 'std', 'min', 'max'],
            'value': [len(igf_found), len(igf_scores), igf_scores.mean(), 
                      igf_scores.median(), igf_scores.std(), igf_scores.min(), igf_scores.max()]
        }).to_csv(RESULTS_DIR / 'PRAD_P01_IGF_summary.csv', index=False)
        
        print(f"✓ IGF Scores: mean={igf_scores.mean():.4f}", flush=True)
        
        # Plots
        print("\nGenerating plots...", flush=True)
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            
            # Histogram
            fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
            ax.hist(igf_scores.values, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
            ax.axvline(igf_scores.mean(), color='red', linestyle='--', linewidth=2, 
                       label=f'Mean: {igf_scores.mean():.3f}')
            ax.axvline(igf_scores.median(), color='green', linestyle='--', linewidth=2, 
                       label=f'Median: {igf_scores.median():.3f}')
            ax.set_xlabel('IGF Pathway Activity Score')
            ax.set_ylabel('Number of Cells')
            ax.set_title('Distribution of IGF Pathway Activity in PRAD_P01')
            ax.legend()
            ax.grid(alpha=0.3)
            plt.savefig(PLOTS_DIR / 'PRAD_P01_IGF_activity_histogram.png', dpi=300, bbox_inches='tight')
            plt.close()
            print(f"✓ Histogram saved", flush=True)
            
            # Top 30
            top_30_idx = np.argsort(igf_scores.values)[-30:][::-1]
            fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
            y_pos = np.arange(30)
            ax.barh(y_pos, igf_scores.values[top_30_idx], color='steelblue', edgecolor='black')
            ax.set_yticks(y_pos)
            ax.set_yticklabels([str(c)[:15] for c in counts_log.columns[top_30_idx]], fontsize=8)
            ax.set_xlabel('IGF Activity Score')
            ax.set_title('Top 30 Cells by IGF Activity')
            ax.invert_yaxis()
            ax.grid(alpha=0.3, axis='x')
            plt.savefig(PLOTS_DIR / 'PRAD_P01_IGF_top30_cells.png', dpi=300, bbox_inches='tight')
            plt.close()
            print(f"✓ Top 30 plot saved", flush=True)
            
        except Exception as e:
            print(f"Warning: Could not generate plots: {e}", flush=True)
    
    print("\n" + "="*60)
    print("✅ PIPELINE COMPLETED SUCCESSFULLY")
    print("="*60)
    print(f"\nResults: {RESULTS_DIR}")
    print(f"Plots: {PLOTS_DIR}")
    
except Exception as e:
    print(f"\n❌ ERROR: {e}", flush=True)
    import traceback
    traceback.print_exc()
    sys.exit(1)
