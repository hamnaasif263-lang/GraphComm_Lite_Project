"""
PRAD_P01 Direct Pipeline - Line-by-line parsing with numpy
No pandas, direct numpy arrays for memory efficiency
"""

import numpy as np
from pathlib import Path
import logging
import gc
import sys

# Setup logging
log_file = Path('prad_p01_direct.log')
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ],
    force=True
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
logger.info("PRAD_P01 DIRECT PIPELINE")
logger.info("Line-by-line parsing with numpy")
logger.info("="*60)

count_file = DATA_DIR / 'GSM4203181_data.raw.matrix.txt'
if not count_file.exists():
    count_file = DATA_DIR / 'GSM4203181_data.matrix.txt'

logger.info(f"\nData file: {count_file}")

# Step 1: Parse header to get cell barcodes and count lines
logger.info("\n" + "="*60)
logger.info("STEP 1: Parsing header and counting lines")
logger.info("="*60)

with open(count_file) as f:
    header_line = f.readline().strip().split('\t')
    cell_barcodes = np.array(header_line[1:], dtype='U50')
    n_cells = len(cell_barcodes)
    
    # Count genes
    n_genes = 0
    for line in f:
        n_genes += 1
        if n_genes % 5000 == 0:
            logger.info(f"  Scanned {n_genes} genes...")

logger.info(f"Found: {n_genes} genes, {n_cells} cells")

# Step 2: Initialize arrays
logger.info("\nInitializing arrays...")
counts = np.zeros((n_genes, n_cells), dtype=np.float32)
gene_names = []

# Step 3: Parse data line by line
logger.info("\n" + "="*60)
logger.info("STEP 2: Loading count data")
logger.info("="*60)

with open(count_file) as f:
    _ = f.readline()  # Skip header
    
    for i, line in enumerate(f):
        if (i + 1) % 5000 == 0:
            logger.info(f"  Loaded {i+1} genes...")
        
        parts = line.strip().split('\t')
        gene_name = parts[0]
        values = np.array([float(x) for x in parts[1:]], dtype=np.float32)
        
        gene_names.append(gene_name)
        counts[i, :] = values

gene_names = np.array(gene_names, dtype='U50')
logger.info(f"Array shape: {counts.shape}")
logger.info(f"Memory: {counts.nbytes / 1e9:.2f} GB")

# Step 3: Cell filtering
logger.info("\n" + "="*60)
logger.info("STEP 3: Cell filtering (>=200 genes)")
logger.info("="*60)

genes_per_cell = np.sum(counts > 0, axis=0)
cells_to_keep = genes_per_cell >= 200

logger.info(f"Cells before: {n_cells}")
logger.info(f"Cells after: {np.sum(cells_to_keep)}")
logger.info(f"Mean genes/cell: {genes_per_cell.mean():.1f}")

counts = counts[:, cells_to_keep]
cell_barcodes = cell_barcodes[cells_to_keep]
gc.collect()

# Step 4: Gene filtering
logger.info("\n" + "="*60)
logger.info("STEP 4: Gene filtering (>=3 cells)")
logger.info("="*60)

cells_per_gene = np.sum(counts > 0, axis=1)
genes_to_keep = cells_per_gene >= 3

logger.info(f"Genes before: {n_genes}")
logger.info(f"Genes after: {np.sum(genes_to_keep)}")

counts = counts[genes_to_keep, :]
gene_names = gene_names[genes_to_keep]
gc.collect()

logger.info(f"Final shape: {counts.shape}")

# Step 5: Normalization
logger.info("\n" + "="*60)
logger.info("STEP 5: Normalization & Log1p")
logger.info("="*60)

lib_size = np.sum(counts, axis=0)
logger.info(f"Library size - Mean: {lib_size.mean():.0f}, Median: {np.median(lib_size):.0f}")

# Normalize to 10,000
norm_factor = 1e4 / lib_size
counts_norm = counts * norm_factor[np.newaxis, :]

# Log1p
counts_log = np.log1p(counts_norm)
logger.info(f"Log-range: {counts_log.min():.4f} to {counts_log.max():.4f}")

del counts, counts_norm
gc.collect()

# Step 6: Save outputs using numpy (faster than pandas)
logger.info("\n" + "="*60)
logger.info("STEP 6: Saving outputs")
logger.info("="*60)

# Save as CSV (genes x cells)
logger.info("Saving normalized matrix...")
with open(RESULTS_DIR / 'scRNA_PRAD_P01_norm.csv', 'w') as f:
    f.write('\t' + '\t'.join(cell_barcodes) + '\n')
    for i, gene in enumerate(gene_names):
        f.write(gene + '\t' + '\t'.join([f"{x:.6f}" for x in counts_log[i, :]]) + '\n')
        if (i + 1) % 5000 == 0:
            logger.info(f"  Saved {i+1} genes...")

logger.info(f"Saved: scRNA_PRAD_P01_norm.csv")

# Step 7: IGF pathway analysis
logger.info("\n" + "="*60)
logger.info("STEP 7: IGF Pathway Analysis")
logger.info("="*60)

IGF_GENES = ["IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1", "MTOR", "FOXO1"]
igf_idx = []
igf_found = []

for g in IGF_GENES:
    mask = gene_names == g
    if np.any(mask):
        igf_idx.append(np.where(mask)[0][0])
        igf_found.append(g)

logger.info(f"IGF genes found: {len(igf_found)}/{len(IGF_GENES)}")
logger.info(f"  {igf_found}")

if igf_idx:
    igf_expr = counts_log[igf_idx, :]
    igf_scores = np.mean(igf_expr, axis=0)
    
    logger.info(f"\nIGF Activity:")
    logger.info(f"  Mean: {igf_scores.mean():.4f}")
    logger.info(f"  Median: {np.median(igf_scores):.4f}")
    logger.info(f"  Std: {igf_scores.std():.4f}")
    
    # Save scores
    logger.info("\nSaving IGF scores...")
    with open(RESULTS_DIR / 'PRAD_P01_IGF_cell_scores.csv', 'w') as f:
        f.write('cell_barcode,igf_activity_score\n')
        for bc, score in zip(cell_barcodes, igf_scores):
            f.write(f'{bc},{score:.6f}\n')
    
    # Summary
    with open(RESULTS_DIR / 'PRAD_P01_IGF_summary.csv', 'w') as f:
        f.write('metric,value\n')
        f.write(f'n_genes,{len(igf_found)}\n')
        f.write(f'n_cells,{len(igf_scores)}\n')
        f.write(f'mean,{igf_scores.mean():.6f}\n')
        f.write(f'median,{np.median(igf_scores):.6f}\n')
        f.write(f'std,{igf_scores.std():.6f}\n')
        f.write(f'min,{igf_scores.min():.6f}\n')
        f.write(f'max,{igf_scores.max():.6f}\n')
    
    logger.info("Saved IGF scores and summary")
    
    # Generate plots
    logger.info("\n" + "="*60)
    logger.info("STEP 8: Generating Plots")
    logger.info("="*60)
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        # Histogram
        logger.info("Generating histogram...")
        fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
        ax.hist(igf_scores, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
        ax.axvline(igf_scores.mean(), color='red', linestyle='--', linewidth=2, 
                   label=f'Mean: {igf_scores.mean():.3f}')
        ax.axvline(np.median(igf_scores), color='green', linestyle='--', linewidth=2,
                   label=f'Median: {np.median(igf_scores):.3f}')
        ax.set_xlabel('IGF Pathway Activity Score')
        ax.set_ylabel('Number of Cells')
        ax.set_title('Distribution of IGF Pathway Activity in PRAD_P01')
        ax.legend()
        ax.grid(alpha=0.3)
        
        plt.savefig(PLOTS_DIR / 'PRAD_P01_IGF_activity_histogram.png', dpi=300, bbox_inches='tight')
        plt.close()
        logger.info("  Saved histogram")
        
        # Top 30 cells
        logger.info("Generating top 30 cells plot...")
        top_30_idx = np.argsort(igf_scores)[-30:][::-1]
        top_30_scores = igf_scores[top_30_idx]
        top_30_cells = [str(c)[:15] for c in cell_barcodes[top_30_idx]]
        
        fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
        y_pos = np.arange(len(top_30_scores))
        ax.barh(y_pos, top_30_scores, color='steelblue', edgecolor='black')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_30_cells, fontsize=8)
        ax.set_xlabel('IGF Pathway Activity Score')
        ax.set_title('Top 30 Cells by IGF Pathway Activity')
        ax.invert_yaxis()
        ax.grid(alpha=0.3, axis='x')
        
        plt.savefig(PLOTS_DIR / 'PRAD_P01_IGF_top30_cells.png', dpi=300, bbox_inches='tight')
        plt.close()
        logger.info("  Saved top 30 plot")
        
    except Exception as e:
        logger.warning(f"Could not generate plots: {e}")

# Summary report
logger.info("\n" + "="*60)
logger.info("PIPELINE COMPLETED")
logger.info("="*60)

summary = f"""PRAD_P01 PREPROCESSING SUMMARY
{'='*50}

Dataset Information:
  Cells: {counts_log.shape[1]:,}
  Genes: {counts_log.shape[0]:,}

Filtering:
  Mean genes/cell: {genes_per_cell[cells_to_keep].mean():.1f}
  Mean cells/gene: {cells_per_gene[genes_to_keep].mean():.1f}

Library Size:
  Mean: {lib_size.mean():.0f}
  Median: {np.median(lib_size):.0f}
  Min: {lib_size.min():.0f}
  Max: {lib_size.max():.0f}

Output Files:
  - scRNA_PRAD_P01_norm.csv
  - PRAD_P01_IGF_cell_scores.csv
  - PRAD_P01_IGF_summary.csv
  - Plots in {PLOTS_DIR}
"""

with open(RESULTS_DIR / 'scRNA_PRAD_P01_summary.txt', 'w') as f:
    f.write(summary)

logger.info(summary)
print("\n✅ Pipeline completed successfully!")
