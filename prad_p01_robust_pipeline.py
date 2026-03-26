"""
Robust PRAD_P01 Pipeline - Memory-efficient processing
Handles large scRNA-seq count matrices
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
import logging
from datetime import datetime

warnings.filterwarnings('ignore')

# Setup logging
log_file = Path('prad_p01_robust.log')
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
logger.info("PRAD_P01 ROBUST PIPELINE STARTED")
logger.info("="*60)

# Step 1: Load the count matrix (use the raw matrix)
logger.info("\nStep 1: Loading count matrix...")
count_file = DATA_DIR / 'GSM4203181_data.raw.matrix.txt'

if not count_file.exists():
    count_file = DATA_DIR / 'GSM4203181_data.matrix.txt'

logger.info(f"Loading from: {count_file}")
logger.info(f"File size: {count_file.stat().st_size / 1e9:.2f} GB")

# Use chunks to read efficiently
logger.info("Reading count matrix in chunks...")
chunk_size = 5000
chunks = []
total_rows = 0

for chunk in pd.read_csv(count_file, sep='\t', index_col=0, chunksize=chunk_size, 
                          dtype='float32', low_memory=False):
    chunks.append(chunk)
    total_rows += len(chunk)
    if total_rows % 10000 == 0:
        logger.info(f"  Loaded {total_rows} genes...")

logger.info(f"Concatenating {len(chunks)} chunks...")
counts = pd.concat(chunks, axis=0)

logger.info(f"Count matrix shape: {counts.shape} (genes x cells)")
logger.info(f"Data type: {counts.dtypes[0]}")
logger.info(f"Memory usage: {counts.memory_usage(deep=True).sum() / 1e9:.2f} GB")

# Step 2: Filter cells with <200 detected genes
logger.info("\nStep 2: Filtering cells with <200 detected genes...")
initial_cells = counts.shape[1]

# Calculate genes per cell (non-zero counts)
genes_per_cell = (counts > 0).sum(axis=0)
min_genes = 200
cells_to_keep = genes_per_cell >= min_genes

counts_filt = counts.loc[:, cells_to_keep]
logger.info(f"  Initial cells: {initial_cells}")
logger.info(f"  Cells after filtering: {counts_filt.shape[1]}")
logger.info(f"  Cells removed: {initial_cells - counts_filt.shape[1]} ({100*(initial_cells - counts_filt.shape[1])/initial_cells:.1f}%)")

# Step 3: Filter genes expressed in <3 cells
logger.info("\nStep 3: Filtering genes expressed in <3 cells...")
initial_genes = counts_filt.shape[0]

cells_per_gene = (counts_filt > 0).sum(axis=1)
min_cells = 3
genes_to_keep = cells_per_gene >= min_cells

counts_filt = counts_filt.loc[genes_to_keep, :]
logger.info(f"  Initial genes: {initial_genes}")
logger.info(f"  Genes after filtering: {counts_filt.shape[0]}")
logger.info(f"  Genes removed: {initial_genes - counts_filt.shape[0]} ({100*(initial_genes - counts_filt.shape[0])/initial_genes:.1f}%)")

# Save filtered matrix
logger.info("\nSaving filtered count matrix...")
filtered_file = RESULTS_DIR / 'scRNA_PRAD_P01_filtered.csv'
counts_filt.to_csv(filtered_file)
logger.info(f"  Saved: {filtered_file}")

# Step 4: Normalize to 10,000 counts per cell
logger.info("\nStep 4: Normalizing to 10,000 counts per cell...")
# Calculate library size for each cell
lib_size = counts_filt.sum(axis=0)
logger.info(f"  Library size range: {lib_size.min():.0f} - {lib_size.max():.0f}")
logger.info(f"  Library size mean: {lib_size.mean():.0f}")

# Normalize
norm_factor = 1e4 / lib_size
counts_norm = counts_filt.mul(norm_factor, axis=1)
logger.info(f"  Normalized data range: {counts_norm.min().min():.6f} - {counts_norm.max().max():.6f}")

# Step 5: Apply log1p transformation
logger.info("\nStep 5: Applying log1p transformation...")
counts_log = np.log1p(counts_norm)
logger.info(f"  Log-transformed data range: {counts_log.min().min():.6f} - {counts_log.max().max():.6f}")

# Save normalized matrix
logger.info("\nSaving normalized count matrix...")
norm_file = RESULTS_DIR / 'scRNA_PRAD_P01_norm.csv'
counts_log.to_csv(norm_file)
logger.info(f"  Saved: {norm_file}")

# Step 6: Create summary statistics
logger.info("\nStep 6: Generating summary statistics...")
summary_stats = {
    'n_cells': counts_filt.shape[1],
    'n_genes': counts_filt.shape[0],
    'genes_per_cell_mean': genes_per_cell[cells_to_keep].mean(),
    'genes_per_cell_median': genes_per_cell[cells_to_keep].median(),
    'genes_per_cell_min': genes_per_cell[cells_to_keep].min(),
    'genes_per_cell_max': genes_per_cell[cells_to_keep].max(),
    'lib_size_mean': lib_size.mean(),
    'lib_size_median': lib_size.median(),
    'lib_size_min': lib_size.min(),
    'lib_size_max': lib_size.max(),
    'cells_removed_genes': initial_cells - counts_filt.shape[1],
    'genes_removed_cells': initial_genes - counts_filt.shape[0]
}

summary_file = RESULTS_DIR / 'scRNA_PRAD_P01_summary.txt'
with open(summary_file, 'w') as f:
    f.write("PRAD_P01 PREPROCESSING SUMMARY\n")
    f.write("="*50 + "\n\n")
    for key, val in summary_stats.items():
        if isinstance(val, float):
            f.write(f"{key}: {val:.2f}\n")
        else:
            f.write(f"{key}: {val}\n")

logger.info(f"  Saved: {summary_file}")
logger.info(f"\n{summary_stats}")

# Step 7: Create visualization sample (3,000 random cells)
logger.info("\nStep 7: Creating visualization sample (3,000 random cells)...")
n_vis = min(3000, counts_log.shape[1])
vis_cells = np.random.choice(counts_log.shape[1], n_vis, replace=False)
counts_vis = counts_log.iloc[:, vis_cells]

vis_file = RESULTS_DIR / 'scRNA_PRAD_P01_vis_sample.csv'
counts_vis.to_csv(vis_file)
logger.info(f"  Saved: {vis_file}")
logger.info(f"  Visualization sample shape: {counts_vis.shape}")

# Step 8: IGF pathway analysis
logger.info("\n" + "="*60)
logger.info("STEP 8: IGF PATHWAY ANALYSIS")
logger.info("="*60)

IGF_GENES = [
    "IGF1", "IGF2", "IGF1R",
    "IRS1", "IRS2",
    "AKT1", "MTOR", "FOXO1"
]

logger.info(f"\nIGF Pathway Genes: {IGF_GENES}")

# Find IGF genes in the matrix
igf_genes_found = [g for g in IGF_GENES if g in counts_log.index]
igf_genes_missing = [g for g in IGF_GENES if g not in counts_log.index]

logger.info(f"\nGenes found: {len(igf_genes_found)}/{len(IGF_GENES)}")
logger.info(f"  {igf_genes_found}")

if igf_genes_missing:
    logger.info(f"Genes not found: {len(igf_genes_missing)}")
    logger.info(f"  {igf_genes_missing}")

if len(igf_genes_found) > 0:
    # Extract IGF genes
    igf_expr = counts_log.loc[igf_genes_found, :]
    
    # Calculate mean expression per cell (IGF pathway activity score)
    igf_scores = igf_expr.mean(axis=0)
    
    logger.info(f"\nIGF Pathway Activity Scores:")
    logger.info(f"  Mean: {igf_scores.mean():.4f}")
    logger.info(f"  Median: {igf_scores.median():.4f}")
    logger.info(f"  Std: {igf_scores.std():.4f}")
    logger.info(f"  Min: {igf_scores.min():.4f}")
    logger.info(f"  Max: {igf_scores.max():.4f}")
    
    # Save cell scores
    igf_cell_scores = pd.DataFrame({
        'cell_barcode': counts_log.columns,
        'igf_activity_score': igf_scores.values
    })
    
    igf_scores_file = RESULTS_DIR / 'PRAD_P01_IGF_cell_scores.csv'
    igf_cell_scores.to_csv(igf_scores_file, index=False)
    logger.info(f"\nSaved: {igf_scores_file}")
    
    # Save summary
    igf_summary = pd.DataFrame({
        'metric': ['n_genes_used', 'n_cells', 'mean_activity', 'median_activity', 'std_activity', 'min_activity', 'max_activity'],
        'value': [len(igf_genes_found), len(igf_scores), igf_scores.mean(), igf_scores.median(), 
                  igf_scores.std(), igf_scores.min(), igf_scores.max()]
    })
    
    igf_summary_file = RESULTS_DIR / 'PRAD_P01_IGF_summary.csv'
    igf_summary.to_csv(igf_summary_file, index=False)
    logger.info(f"Saved: {igf_summary_file}")
    
    # Step 9: Generate IGF pathway plots
    logger.info("\n" + "="*60)
    logger.info("STEP 9: GENERATING IGF PATHWAY PLOTS")
    logger.info("="*60)
    
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        # Plot 1: Histogram of IGF pathway activity
        logger.info("\nGenerating IGF pathway activity histogram...")
        fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
        
        ax.hist(igf_scores, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
        ax.axvline(igf_scores.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {igf_scores.mean():.3f}')
        ax.axvline(igf_scores.median(), color='green', linestyle='--', linewidth=2, label=f'Median: {igf_scores.median():.3f}')
        ax.set_xlabel('IGF Pathway Activity Score')
        ax.set_ylabel('Number of Cells')
        ax.set_title('Distribution of IGF Pathway Activity in PRAD_P01')
        ax.legend()
        ax.grid(alpha=0.3)
        
        hist_file = PLOTS_DIR / 'PRAD_P01_IGF_activity_histogram.png'
        plt.savefig(hist_file, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved: {hist_file}")
        
        # Plot 2: Top 30 cells by IGF activity
        logger.info("Generating top 30 cells by IGF activity...")
        top_30_idx = np.argsort(igf_scores.values)[-30:][::-1]
        top_30_scores = igf_scores.values[top_30_idx]
        top_30_cells = counts_log.columns[top_30_idx].str[:15]  # Truncate barcodes
        
        fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
        
        y_pos = np.arange(len(top_30_scores))
        ax.barh(y_pos, top_30_scores, color='steelblue', edgecolor='black')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_30_cells, fontsize=8)
        ax.set_xlabel('IGF Pathway Activity Score')
        ax.set_title('Top 30 Cells by IGF Pathway Activity')
        ax.invert_yaxis()
        ax.grid(alpha=0.3, axis='x')
        
        top30_file = PLOTS_DIR / 'PRAD_P01_IGF_top30_cells.png'
        plt.savefig(top30_file, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved: {top30_file}")
        
        logger.info("✓ IGF pathway plots generated successfully")
        
    except ImportError as e:
        logger.warning(f"Could not generate plots (matplotlib/seaborn not installed): {e}")
else:
    logger.warning("No IGF genes found in the matrix - skipping IGF analysis")

# Step 10: GraphComm-Lite Analysis Placeholder
logger.info("\n" + "="*60)
logger.info("STEP 10: GRAPHCOMM-LITE COMMUNICATION ANALYSIS")
logger.info("="*60)

logger.info("\nPreparing data for GraphComm-Lite analysis...")
logger.info(f"  Input matrix shape: {counts_log.shape}")
logger.info(f"  Ready for cell-cell communication inference")
logger.info(f"  Output directory: {PLOTS_DIR}")

# Save a sample for GraphComm
graphcomm_input = RESULTS_DIR / 'scRNA_PRAD_P01_graphcomm_input.csv'
counts_log.to_csv(graphcomm_input)
logger.info(f"  Saved GraphComm input: {graphcomm_input}")

logger.info("\n" + "="*60)
logger.info("PIPELINE COMPLETED SUCCESSFULLY")
logger.info("="*60)

logger.info(f"\nOutput Summary:")
logger.info(f"  Results directory: {RESULTS_DIR}")
logger.info(f"  Plots directory: {PLOTS_DIR}")
logger.info(f"  Log file: {log_file}")

logger.info("\nGenerated Files:")
logger.info(f"  - scRNA_PRAD_P01_filtered.csv")
logger.info(f"  - scRNA_PRAD_P01_norm.csv")
logger.info(f"  - scRNA_PRAD_P01_vis_sample.csv")
logger.info(f"  - scRNA_PRAD_P01_summary.txt")
logger.info(f"  - PRAD_P01_IGF_cell_scores.csv")
logger.info(f"  - PRAD_P01_IGF_summary.csv")
logger.info(f"  - PRAD_P01_IGF_activity_histogram.png (300 dpi)")
logger.info(f"  - PRAD_P01_IGF_top30_cells.png (300 dpi)")

logger.info("\nNext Steps:")
logger.info(f"  Run GraphComm-Lite on the normalized data")
logger.info(f"  Classify cells as Autocrine/Sender/Receiver/Inactive")
logger.info(f"  Generate communication network visualizations")
