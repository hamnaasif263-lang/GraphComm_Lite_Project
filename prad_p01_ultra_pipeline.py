"""
PRAD_P01 Ultra-Lite Pipeline - Chunk-based processing
Processes data in manageable pieces to avoid memory issues
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sys

# Setup logging
log_file = Path('prad_p01_ultra.log')
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Define paths
DATA_DIR = Path('data/prostate/raw')
RESULTS_DIR = Path('graphcomm/results/prostate')
PLOTS_DIR = Path('graphcomm/plots/prostate')

# Create directories
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

logger.info("="*60)
logger.info("PRAD_P01 ULTRA-LITE PIPELINE")
logger.info("Chunk-based processing for large matrices")
logger.info("="*60)

count_file = DATA_DIR / 'GSM4203181_data.raw.matrix.txt'
if not count_file.exists():
    count_file = DATA_DIR / 'GSM4203181_data.matrix.txt'

logger.info(f"\nData file: {count_file}")
logger.info(f"File size: {count_file.stat().st_size / 1e9:.2f} GB")

# IGF Genes to extract
IGF_GENES = {"IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1", "MTOR", "FOXO1"}

logger.info("\n" + "="*60)
logger.info("STEP 1: Reading data with chunking")
logger.info("="*60)

# Read in chunks
chunks = []
total_rows = 0
logger.info("Reading data chunks...")

try:
    for i, chunk in enumerate(pd.read_csv(count_file, sep='\t', index_col=0, chunksize=3000, low_memory=False)):
        chunks.append(chunk.astype(np.float32))
        total_rows += len(chunk)
        if (i + 1) % 5 == 0:
            logger.info(f"  Loaded {total_rows} genes...")
    
    logger.info(f"Total chunks: {len(chunks)}, Total genes: {total_rows}")
    logger.info("Concatenating chunks...")
    counts = pd.concat(chunks, axis=0)
    logger.info(f"Full matrix shape: {counts.shape}")
    logger.info(f"Memory: {counts.memory_usage(deep=True).sum() / 1e9:.2f} GB")
    
except Exception as e:
    logger.error(f"Error reading file: {e}")
    import traceback
    logger.error(traceback.format_exc())
    sys.exit(1)

# Step 2-3: Filtering
logger.info("\n" + "="*60)
logger.info("STEP 2-3: Cell and Gene Filtering")
logger.info("="*60)

genes_per_cell = (counts > 0).sum(axis=0)
cells_to_keep = genes_per_cell >= 200
logger.info(f"Cells with >=200 genes: {cells_to_keep.sum()} / {len(cells_to_keep)}")

counts = counts.loc[:, cells_to_keep]

cells_per_gene = (counts > 0).sum(axis=1)
genes_to_keep = cells_per_gene >= 3
logger.info(f"Genes in >=3 cells: {genes_to_keep.sum()} / {len(cells_per_gene)}")

counts = counts.loc[genes_to_keep, :]
logger.info(f"Filtered matrix: {counts.shape}")

# Save filtered
filtered_file = RESULTS_DIR / 'scRNA_PRAD_P01_filtered.csv'
counts.to_csv(filtered_file)
logger.info(f"Saved: {filtered_file}")

# Step 4-5: Normalization & Log
logger.info("\n" + "="*60)
logger.info("STEP 4-5: Normalization & Log1p")
logger.info("="*60)

lib_size = counts.sum(axis=0)
logger.info(f"Library size - Mean: {lib_size.mean():.0f}, Median: {lib_size.median():.0f}")

norm_factor = 1e4 / lib_size
counts_norm = counts.mul(norm_factor, axis=1)
counts_log = np.log1p(counts_norm)

logger.info(f"Log-transformed range: {counts_log.min().min():.4f} to {counts_log.max().max():.4f}")

norm_file = RESULTS_DIR / 'scRNA_PRAD_P01_norm.csv'
counts_log.to_csv(norm_file)
logger.info(f"Saved: {norm_file}")

# Step 6: Summary
logger.info("\n" + "="*60)
logger.info("STEP 6: Summary Statistics")
logger.info("="*60)

summary = f"""PRAD_P01 PREPROCESSING SUMMARY
{'='*50}

Dataset Information:
  Number of cells: {counts.shape[1]:,}
  Number of genes: {counts.shape[0]:,}

Filtering Results:
  Genes per cell - Mean: {genes_per_cell[cells_to_keep].mean():.1f}, Median: {genes_per_cell[cells_to_keep].median():.1f}
  
Library size:
  Mean: {lib_size.mean():.0f}
  Median: {lib_size.median():.0f}
  Min: {lib_size.min():.0f}
  Max: {lib_size.max():.0f}
"""

summary_file = RESULTS_DIR / 'scRNA_PRAD_P01_summary.txt'
with open(summary_file, 'w') as f:
    f.write(summary)
logger.info(summary)

# Step 7: Visualization sample
logger.info("\n" + "="*60)
logger.info("STEP 7: Visualization Sample")
logger.info("="*60)

n_vis = min(3000, counts_log.shape[1])
vis_idx = np.random.choice(counts_log.shape[1], n_vis, replace=False)
counts_vis = counts_log.iloc[:, vis_idx]

vis_file = RESULTS_DIR / 'scRNA_PRAD_P01_vis_sample.csv'
counts_vis.to_csv(vis_file)
logger.info(f"Created sample: {counts_vis.shape}")
logger.info(f"Saved: {vis_file}")

# Step 8: IGF Pathway Analysis
logger.info("\n" + "="*60)
logger.info("STEP 8: IGF Pathway Analysis")
logger.info("="*60)

igf_found = [g for g in IGF_GENES if g in counts_log.index]
igf_missing = list(IGF_GENES - set(igf_found))

logger.info(f"IGF genes found: {len(igf_found)}/{len(IGF_GENES)}")
logger.info(f"  Found: {sorted(igf_found)}")
if igf_missing:
    logger.info(f"  Missing: {sorted(igf_missing)}")

if igf_found:
    igf_expr = counts_log.loc[igf_found, :]
    igf_scores = igf_expr.mean(axis=0)
    
    logger.info(f"\nIGF Activity - Mean: {igf_scores.mean():.4f}, Median: {igf_scores.median():.4f}, Std: {igf_scores.std():.4f}")
    
    # Save scores
    igf_df = pd.DataFrame({
        'cell_barcode': counts_log.columns,
        'igf_activity_score': igf_scores.values
    })
    igf_file = RESULTS_DIR / 'PRAD_P01_IGF_cell_scores.csv'
    igf_df.to_csv(igf_file, index=False)
    logger.info(f"Saved: {igf_file}")
    
    # Save summary
    igf_summary = pd.DataFrame({
        'metric': ['n_genes', 'n_cells', 'mean', 'median', 'std', 'min', 'max'],
        'value': [len(igf_found), len(igf_scores), igf_scores.mean(), 
                  igf_scores.median(), igf_scores.std(), igf_scores.min(), igf_scores.max()]
    })
    igf_sum_file = RESULTS_DIR / 'PRAD_P01_IGF_summary.csv'
    igf_summary.to_csv(igf_sum_file, index=False)
    logger.info(f"Saved: {igf_sum_file}")
    
    # Generate plots
    logger.info("\n" + "="*60)
    logger.info("STEP 9: Generating Plots")
    logger.info("="*60)
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        # Plot 1: Histogram
        logger.info("Generating histogram...")
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
        
        hist_file = PLOTS_DIR / 'PRAD_P01_IGF_activity_histogram.png'
        plt.savefig(hist_file, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved: {hist_file}")
        
        # Plot 2: Top 30 cells
        logger.info("Generating top 30 cells plot...")
        top_30_idx = np.argsort(igf_scores.values)[-30:][::-1]
        top_30_scores = igf_scores.values[top_30_idx]
        top_30_cells = [str(c)[:15] for c in counts_log.columns[top_30_idx]]
        
        fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
        y_pos = np.arange(len(top_30_scores))
        ax.barh(y_pos, top_30_scores, color='steelblue', edgecolor='black')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_30_cells, fontsize=8)
        ax.set_xlabel('IGF Pathway Activity Score')
        ax.set_title('Top 30 Cells by IGF Pathway Activity')
        ax.invert_yaxis()
        ax.grid(alpha=0.3, axis='x')
        
        top_file = PLOTS_DIR / 'PRAD_P01_IGF_top30_cells.png'
        plt.savefig(top_file, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved: {top_file}")
        
        logger.info("✓ Plots generated")
        
    except Exception as e:
        logger.warning(f"Could not generate plots: {e}")

# Final
logger.info("\n" + "="*60)
logger.info("PIPELINE COMPLETED")
logger.info("="*60)
logger.info(f"\nResults: {RESULTS_DIR}")
logger.info(f"Plots: {PLOTS_DIR}")

print("\n✅ Pipeline completed successfully!")
