"""
PRAD_P01 Smart Sampling Pipeline
Works with a strategic sampling of cells to complete analysis quickly
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging

# Setup
log_file = Path('prad_p01_sampling.log')
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

DATA_DIR = Path('data/prostate/raw')
RESULTS_DIR = Path('graphcomm/results/prostate')
PLOTS_DIR = Path('graphcomm/plots/prostate')
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

logger.info("PRAD_P01 SMART SAMPLING PIPELINE")
logger.info("="*60)

count_file = DATA_DIR / 'GSM4203181_data.raw.matrix.txt'
if not count_file.exists():
    count_file = DATA_DIR / 'GSM4203181_data.matrix.txt'

# Quick header parse
logger.info("Reading header...")
with open(count_file) as f:
    header = f.readline().strip().split('\t')

n_cells_total = len(header) - 1
logger.info(f"Total cells in file: {n_cells_total}")

# Sample cells evenly across the dataset
logger.info("\nSampling 15,000 cells...")
sample_size = min(15000, n_cells_total)
cell_indices = np.linspace(0, n_cells_total - 1, sample_size, dtype=int)
cell_indices = np.unique(cell_indices)

logger.info(f"Sampling {len(cell_indices)} cells")

# Read only sampled cells
logger.info("Reading sampled data...")
chunks = []
for i, chunk in enumerate(pd.read_csv(count_file, sep='\t', index_col=0, usecols=[0] + list(cell_indices + 1), 
                                       chunksize=2000, low_memory=False)):
    chunks.append(chunk.astype(np.float32))
    if (i + 1) % 5 == 0:
        logger.info(f"  Loaded {i+1} chunks...")

counts = pd.concat(chunks)
logger.info(f"Sampled matrix shape: {counts.shape}")

# Filter
logger.info("\nFiltering (<200 genes/cell, <3 cells/gene)...")
genes_per_cell = (counts > 0).sum(axis=0)
cells_to_keep = genes_per_cell >= 200
counts = counts.loc[:, cells_to_keep]

cells_per_gene = (counts > 0).sum(axis=1)
genes_to_keep = cells_per_gene >= 3
counts = counts.loc[genes_to_keep, :]

logger.info(f"Filtered shape: {counts.shape}")

# Normalize & Log
logger.info("\nNormalizing...")
lib_size = counts.sum(axis=0)
counts_norm = counts.mul(1e4 / lib_size, axis=1)
counts_log = np.log1p(counts_norm)

logger.info(f"Log range: {counts_log.min().min():.4f} to {counts_log.max().max():.4f}")

# Save
logger.info("\nSaving outputs...")
counts_log.to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_norm.csv')
logger.info(f"  Normalized matrix: {counts_log.shape}")

# Summary
with open(RESULTS_DIR / 'scRNA_PRAD_P01_summary.txt', 'w') as f:
    f.write(f"""PRAD_P01 SAMPLING SUMMARY
{'='*50}

Dataset (sampled):
  Cells: {counts.shape[1]:,}
  Genes: {counts.shape[0]:,}
  Sampling rate: {len(cell_indices)/n_cells_total*100:.1f}%

Library size:
  Mean: {lib_size.mean():.0f}
  Median: {lib_size.median():.0f}
""")

logger.info("Summary saved")

# Visualization sample
logger.info("\nCreating visualization sample...")
n_vis = min(3000, counts_log.shape[1])
vis_idx = np.random.choice(counts_log.shape[1], n_vis, replace=False)
counts_log.iloc[:, vis_idx].to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_vis_sample.csv')
logger.info(f"  Created sample: {counts_log.iloc[:, vis_idx].shape}")

# IGF Analysis
logger.info("\nIGF Pathway Analysis...")
IGF_GENES = ["IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1", "MTOR", "FOXO1"]
igf_found = [g for g in IGF_GENES if g in counts_log.index]
logger.info(f"IGF genes found: {len(igf_found)}/{len(IGF_GENES)}")

if igf_found:
    igf_scores = counts_log.loc[igf_found].mean(axis=0)
    
    igf_df = pd.DataFrame({
        'cell_barcode': counts_log.columns,
        'igf_activity_score': igf_scores.values
    })
    igf_df.to_csv(RESULTS_DIR / 'PRAD_P01_IGF_cell_scores.csv', index=False)
    
    igf_summary = pd.DataFrame({
        'metric': ['n_genes', 'n_cells', 'mean', 'median', 'std', 'min', 'max'],
        'value': [len(igf_found), len(igf_scores), igf_scores.mean(), 
                  igf_scores.median(), igf_scores.std(), igf_scores.min(), igf_scores.max()]
    })
    igf_summary.to_csv(RESULTS_DIR / 'PRAD_P01_IGF_summary.csv', index=False)
    
    logger.info(f"  Mean activity: {igf_scores.mean():.4f}")
    logger.info(f"  Saved: IGF scores and summary")
    
    # Plots
    logger.info("\nGenerating plots...")
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
        logger.info("  Histogram saved")
        
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
        logger.info("  Top 30 plot saved")
        
    except Exception as e:
        logger.warning(f"Plots: {e}")

logger.info("\n" + "="*60)
logger.info("COMPLETED SUCCESSFULLY")
logger.info("="*60)
logger.info(f"\nResults: {RESULTS_DIR}")
logger.info(f"Plots: {PLOTS_DIR}")

print("\n✅ Pipeline completed!")
