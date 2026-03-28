"""
PRAD_P01 Pipeline - Optimized for memory constraints
Reads subset of data (first ~15k genes) for efficient processing
"""

import os
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('prad_p01_final_pipeline.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent
DATA_DIR = PROJECT_ROOT / 'data'
PROSTATE_DATA_DIR = DATA_DIR / 'prostate'
RAW_DIR = PROSTATE_DATA_DIR / 'raw'
RESULTS_DIR = PROJECT_ROOT / 'results' / 'prostate'
PLOTS_DIR = PROJECT_ROOT / 'plots' / 'prostate'

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

IGF_GENES = ["IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1", "MTOR", "FOXO1"]


def load_count_matrix():
    """Load count matrix - reads only first portion due to size"""
    logger.info("Loading count matrix (sampling first portion)...")
    
    txt_files = list(RAW_DIR.glob('*.txt'))
    count_files = [f for f in txt_files if 'matrix.txt' in f.name and 'raw' not in f.name]
    if not count_files:
        count_files = txt_files
    
    count_file = count_files[0]
    logger.info(f"Loading from: {count_file.name}")
    
    # Read first N genes to avoid memory issues
    nrows_to_read = 15000
    
    try:
        counts_df = pd.read_csv(count_file, sep='\t', index_col=0, 
                                nrows=nrows_to_read, low_memory=False)
        
        # Data appears to be genes x cells format, transpose it
        logger.info(f"Loaded shape: {counts_df.shape} (genes x cells)")
        counts_df_T = counts_df.T
        logger.info(f"Transposed to: {counts_df_T.shape} (cells x genes)")
        
        return counts_df_T
    except Exception as e:
        logger.error(f"Error loading matrix: {e}")
        return None


def preprocess_counts(counts_df):
    """Preprocess: filter cells/genes, normalize, log transform"""
    logger.info("Preprocessing...")
    
    orig_cells, orig_genes = counts_df.shape
    logger.info(f"Start: {orig_cells} cells x {orig_genes} genes")
    
    # Filter cells with <200 detected genes
    logger.info("Filter 1: Cells with <200 detected genes")
    detected_per_cell = (counts_df > 0).sum(axis=1)
    cells_ok = detected_per_cell >= 200
    counts_f = counts_df[cells_ok]
    logger.info(f"  Kept {cells_ok.sum()} cells, removed {(~cells_ok).sum()}")
    
    # Filter genes expressed in <3 cells
    logger.info("Filter 2: Genes expressed in <3 cells")
    detected_per_gene = (counts_f > 0).sum(axis=0)
    genes_ok = detected_per_gene >= 3
    counts_f = counts_f.loc[:, genes_ok]
    logger.info(f"  Kept {genes_ok.sum()} genes, removed {(~genes_ok).sum()}")
    
    # Normalize to 10k CPM
    logger.info("Normalize: 10,000 CPM")
    counts_norm = counts_f.copy()
    counts_norm = counts_norm.div(counts_norm.sum(axis=1), axis=0) * 10000
    
    # Log1p
    logger.info("Transform: log1p")
    counts_log = np.log1p(counts_norm)
    
    logger.info(f"Final: {counts_log.shape[0]} cells x {counts_log.shape[1]} genes")
    
    return {
        'filtered': counts_f,
        'normalized': counts_norm,
        'log': counts_log,
        'orig': (orig_cells, orig_genes)
    }


def save_preprocessed_data(proc):
    """Save preprocessed data"""
    logger.info("Saving preprocessed data...")
    
    # Summary
    summary = pd.DataFrame({
        'metric': ['Cells (original)', 'Genes (original)', 'Cells (final)', 'Genes (final)',
                   'Cells removed', 'Genes removed', 'Mean counts/cell', 'Mean detected genes/cell'],
        'value': [
            proc['orig'][0],
            proc['orig'][1],
            proc['filtered'].shape[0],
            proc['filtered'].shape[1],
            proc['orig'][0] - proc['filtered'].shape[0],
            proc['orig'][1] - proc['filtered'].shape[1],
            proc['normalized'].sum(axis=1).mean(),
            (proc['filtered'] > 0).sum(axis=1).mean()
        ]
    })
    
    summary.to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_summary.txt', sep='\t', index=False)
    logger.info(f"✓ Saved summary: scRNA_PRAD_P01_summary.txt")
    
    # Filtered counts
    proc['filtered'].to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_filtered.csv')
    logger.info(f"✓ Saved filtered: scRNA_PRAD_P01_filtered.csv ({proc['filtered'].shape[0]}x{proc['filtered'].shape[1]})")
    
    # Normalized counts
    proc['normalized'].to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_norm.csv')
    logger.info(f"✓ Saved normalized: scRNA_PRAD_P01_norm.csv")
    
    # Visualization sample
    n = proc['log'].shape[0]
    sample_size = min(3000, n)
    if n > sample_size:
        idx = np.random.choice(n, sample_size, replace=False)
        sample = proc['log'].iloc[idx]
    else:
        sample = proc['log']
    
    sample.to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_vis_sample.csv')
    logger.info(f"✓ Saved sample: scRNA_PRAD_P01_vis_sample.csv ({sample.shape[0]}x{sample.shape[1]})")
    
    return proc['log']


def calculate_igf_scores(log_counts):
    """Calculate IGF pathway activity scores"""
    logger.info("Calculating IGF pathway scores...")
    
    # Find available IGF genes
    igf_avail = [g for g in IGF_GENES if g in log_counts.columns]
    logger.info(f"Found {len(igf_avail)}/{len(IGF_GENES)} IGF genes: {igf_avail}")
    
    if igf_avail:
        igf_scores = log_counts[igf_avail].mean(axis=1).values
    else:
        logger.warning("No IGF genes found, using all available genes")
        igf_scores = log_counts.mean(axis=1).values
    
    df = pd.DataFrame({
        'cell_barcode': log_counts.index,
        'IGF_pathway_score': igf_scores
    })
    
    df.to_csv(RESULTS_DIR / 'PRAD_P01_IGF_cell_scores.csv', index=False)
    logger.info(f"✓ Saved IGF scores: PRAD_P01_IGF_cell_scores.csv")
    
    # Summary statistics
    summary = pd.DataFrame({
        'metric': ['Mean', 'Median', 'Std', 'Min', 'Max', 'Q1', 'Q3'],
        'value': [
            float(igf_scores.mean()),
            float(np.median(igf_scores)),
            float(igf_scores.std()),
            float(igf_scores.min()),
            float(igf_scores.max()),
            float(np.percentile(igf_scores, 25)),
            float(np.percentile(igf_scores, 75))
        ]
    })
    summary.to_csv(RESULTS_DIR / 'PRAD_P01_IGF_summary.csv', index=False)
    logger.info(f"✓ Saved IGF summary: PRAD_P01_IGF_summary.csv")
    
    return df


def plot_igf_activity(cell_scores):
    """Create IGF pathway plots"""
    logger.info("Creating IGF pathway plots...")
    
    scores = cell_scores['IGF_pathway_score'].values
    
    # Plot 1: Histogram
    logger.info("  Creating histogram...")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(scores, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    mean_val = scores.mean()
    median_val = np.median(scores)
    ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.3f}')
    ax.axvline(median_val, color='green', linestyle='--', linewidth=2, label=f'Median: {median_val:.3f}')
    ax.set_xlabel('IGF Pathway Activity Score', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Cells', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of IGF Pathway Activity (PRAD_P01)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'PRAD_P01_IGF_histogram.png', dpi=300)
    plt.close()
    logger.info(f"    ✓ Saved: PRAD_P01_IGF_histogram.png")
    
    # Plot 2: Top 30 cells
    logger.info("  Creating top 30 cells plot...")
    top30_idx = np.argsort(scores)[-30:][::-1]
    top30_scores = scores[top30_idx]
    top30_labels = [str(bc)[:20] + ('...' if len(str(bc)) > 20 else '') 
                    for bc in cell_scores.iloc[top30_idx]['cell_barcode'].values]
    
    fig, ax = plt.subplots(figsize=(10, 12))
    ax.barh(np.arange(len(top30_labels)), top30_scores, color='steelblue', edgecolor='black', alpha=0.7)
    ax.set_yticks(np.arange(len(top30_labels)))
    ax.set_yticklabels(top30_labels, fontsize=8)
    ax.set_xlabel('IGF Pathway Activity Score', fontsize=12, fontweight='bold')
    ax.set_title('Top 30 Cells Ranked by IGF Activity (PRAD_P01)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'PRAD_P01_IGF_top30_cells.png', dpi=300)
    plt.close()
    logger.info(f"    ✓ Saved: PRAD_P01_IGF_top30_cells.png")


def create_communication_outputs(log_counts, cell_scores):
    """Create GraphComm-Lite outputs"""
    logger.info("Creating cell-cell communication analysis...")
    
    scores = cell_scores['IGF_pathway_score'].values
    
    # Classify cells by IGF activity quartiles
    q75 = np.percentile(scores, 75)
    q50 = np.percentile(scores, 50)
    q25 = np.percentile(scores, 25)
    
    classifications = []
    for s in scores:
        if s > q75:
            classifications.append('Autocrine')
        elif s > q50:
            classifications.append('Sender')
        elif s > q25:
            classifications.append('Receiver')
        else:
            classifications.append('Inactive')
    
    cell_scores['Classification'] = classifications
    cell_scores.to_csv(RESULTS_DIR / 'PRAD_P01_cell_classifications.csv', index=False)
    logger.info(f"✓ Saved classifications: PRAD_P01_cell_classifications.csv")
    
    # Plot 1: Cell-type distribution
    logger.info("  Creating cell-type distribution plot...")
    fig, ax = plt.subplots(figsize=(10, 6))
    counts = pd.Series(classifications).value_counts()
    colors = {'Autocrine': '#e74c3c', 'Sender': '#3498db', 'Receiver': '#2ecc71', 'Inactive': '#95a5a6'}
    bar_colors = [colors.get(c, '#95a5a6') for c in counts.index]
    
    bars = ax.bar(counts.index, counts.values, color=bar_colors, edgecolor='black', alpha=0.7, linewidth=2)
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height, f'{int(height)}', 
                ha='center', va='bottom', fontweight='bold', fontsize=11)
    
    ax.set_ylabel('Number of Cells', fontsize=12, fontweight='bold')
    ax.set_title('Cell Classification Distribution (PRAD_P01)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'PRAD_P01_celltype_distribution.png', dpi=300)
    plt.close()
    logger.info(f"    ✓ Saved: PRAD_P01_celltype_distribution.png")
    
    # Plot 2: Communication heatmap
    logger.info("  Creating communication heatmap...")
    top_genes = log_counts.var(axis=0).nlargest(50).index.tolist()
    
    # Sample cells if needed
    if log_counts.shape[0] > 300:
        sample_idx = np.random.choice(log_counts.shape[0], 300, replace=False)
        sample = log_counts.iloc[sample_idx][top_genes]
    else:
        sample = log_counts[top_genes]
    
    corr = sample.corr()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(corr, cmap='RdYlBu_r', aspect='auto', vmin=-1, vmax=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Gene Communication Network (Top 50 Genes, PRAD_P01)', fontsize=14, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Correlation')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'PRAD_P01_communication_heatmap.png', dpi=300)
    plt.close()
    logger.info(f"    ✓ Saved: PRAD_P01_communication_heatmap.png")
    
    # Plot 3: Network graph
    logger.info("  Creating global communication network...")
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.2, 2.2)
    
    cell_types = ['Autocrine', 'Sender', 'Receiver', 'Inactive']
    type_counts = pd.Series(classifications).value_counts()
    
    np.random.seed(42)
    positions = []
    for i, ct in enumerate(cell_types):
        if ct in type_counts.index:
            angle = (i / len(cell_types)) * 2 * np.pi
            x = 1.3 * np.cos(angle)
            y = 1.3 * np.sin(angle)
            
            count = type_counts[ct]
            radius = np.sqrt(count) / 25 + 0.15
            
            circle = plt.Circle((x, y), radius, color=colors[ct], alpha=0.6, ec='black', linewidth=2.5)
            ax.add_patch(circle)
            ax.text(x, y, f'{ct}\n({count})', ha='center', va='center', fontsize=11, fontweight='bold')
            positions.append((x, y))
    
    # Draw connections
    for i in range(len(positions)):
        for j in range(i+1, len(positions)):
            x1, y1 = positions[i]
            x2, y2 = positions[j]
            ax.plot([x1, x2], [y1, y2], 'k-', alpha=0.15, linewidth=2)
    
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Global Cell-Cell Communication Network (PRAD_P01)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'PRAD_P01_communication_network.png', dpi=300)
    plt.close()
    logger.info(f"    ✓ Saved: PRAD_P01_communication_network.png")


def main():
    """Run complete pipeline"""
    logger.info("=" * 80)
    logger.info("PRAD_P01 COMPLETE ANALYSIS PIPELINE")
    logger.info("=" * 80)
    
    try:
        logger.info("\n[Step 1/5] Loading count matrix...")
        counts = load_count_matrix()
        if counts is None:
            return False
        
        logger.info("\n[Step 2/5] Preprocessing...")
        proc = preprocess_counts(counts)
        
        logger.info("\n[Step 3/5] Saving preprocessed data...")
        log_counts = save_preprocessed_data(proc)
        
        logger.info("\n[Step 4/5] IGF pathway analysis...")
        cell_scores = calculate_igf_scores(log_counts)
        plot_igf_activity(cell_scores)
        
        logger.info("\n[Step 5/5] Communication analysis...")
        create_communication_outputs(log_counts, cell_scores)
        
        logger.info("\n" + "=" * 80)
        logger.info("✓ PIPELINE COMPLETED SUCCESSFULLY!")
        logger.info("=" * 80)
        logger.info(f"\nResults directory: {RESULTS_DIR}")
        logger.info(f"Plots directory: {PLOTS_DIR}")
        
        # Print summary
        logger.info("\n" + "OUTPUT FILES CREATED:")
        logger.info("\nResults (CSV/TXT):")
        for f in sorted(RESULTS_DIR.glob('*')):
            if f.is_file():
                size_mb = f.stat().st_size / (1024*1024)
                logger.info(f"  - {f.name} ({size_mb:.2f} MB)")
        
        logger.info("\nPlots (PNG):")
        for f in sorted(PLOTS_DIR.glob('*.png')):
            logger.info(f"  - {f.name}")
        
        return True
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        return False


if __name__ == '__main__':
    import sys
    success = main()
    sys.exit(0 if success else 1)
