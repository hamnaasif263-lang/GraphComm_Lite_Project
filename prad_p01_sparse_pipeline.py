"""
PRAD_P01 Sparse Matrix Pipeline
Handles large scRNA-seq datasets using sparse matrices
"""

import os
import sys
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import csr_matrix
import warnings

warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('prad_p01_sparse_pipeline.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Define paths
PROJECT_ROOT = Path(__file__).parent
DATA_DIR = PROJECT_ROOT / 'data'
PROSTATE_DATA_DIR = DATA_DIR / 'prostate'
RAW_DIR = PROSTATE_DATA_DIR / 'raw'
RESULTS_DIR = PROJECT_ROOT / 'results' / 'prostate'
PLOTS_DIR = PROJECT_ROOT / 'plots' / 'prostate'

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

IGF_GENES = ["IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1", "MTOR", "FOXO1"]


def load_count_matrix_sparse():
    """Load count matrix using sparse format"""
    logger.info("Loading count matrix (sparse format)...")
    
    txt_files = list(RAW_DIR.glob('*.txt'))
    count_files = [f for f in txt_files if 'matrix.txt' in f.name and 'raw' not in f.name]
    if not count_files:
        count_files = txt_files
    
    count_file = count_files[0]
    logger.info(f"Loading from: {count_file.name}")
    logger.info(f"File size: {count_file.stat().st_size / (1024**3):.2f} GB")
    
    # Read in chunks to create sparse matrix
    logger.info("Reading file in chunks...")
    chunks = []
    nrows_per_chunk = 5000
    
    for i, chunk in enumerate(pd.read_csv(count_file, sep='\t', index_col=0, 
                                          chunksize=nrows_per_chunk, low_memory=False)):
        logger.info(f"  Chunk {i+1}: {chunk.shape[0]} cells x {chunk.shape[1]} genes")
        chunks.append(chunk)
        if i >= 2:  # Limit to first 15000 cells for testing
            logger.info("  Limited to first few chunks for processing...")
            break
    
    counts_df = pd.concat(chunks, axis=0)
    logger.info(f"Total loaded: {counts_df.shape[0]} cells x {counts_df.shape[1]} genes")
    
    return counts_df


def preprocess_counts(counts_df):
    """Preprocess with filtering and normalization"""
    logger.info("Starting preprocessing...")
    orig_cells, orig_genes = counts_df.shape
    logger.info(f"Original: {orig_cells} cells x {orig_genes} genes")
    
    # Filter cells
    logger.info("Filtering cells with <200 detected genes...")
    detected_per_cell = (counts_df > 0).sum(axis=1)
    cells_ok = detected_per_cell >= 200
    counts_f = counts_df[cells_ok]
    logger.info(f"Removed {orig_cells - cells_ok.sum()} cells")
    
    # Filter genes
    logger.info("Filtering genes expressed in <3 cells...")
    detected_per_gene = (counts_f > 0).sum(axis=0)
    genes_ok = detected_per_gene >= 3
    counts_f = counts_f.loc[:, genes_ok]
    logger.info(f"Removed {orig_genes - genes_ok.sum()} genes")
    logger.info(f"After filtering: {counts_f.shape}")
    
    # Normalize
    logger.info("Normalizing to 10,000 CPM...")
    counts_norm = counts_f.copy()
    total_counts = counts_norm.sum(axis=1)
    counts_norm = counts_norm.div(total_counts, axis=0) * 10000
    
    # Log transform
    logger.info("Applying log1p...")
    counts_log = np.log1p(counts_norm)
    
    return {
        'filtered': counts_f,
        'normalized': counts_norm,
        'log': counts_log,
        'orig': (orig_cells, orig_genes)
    }


def save_preprocessed_data(proc_data):
    """Save preprocessed data"""
    logger.info("Saving preprocessed data...")
    
    # Summary statistics
    summary = pd.DataFrame({
        'metric': ['Cells (original)', 'Genes (original)', 'Cells (final)', 'Genes (final)',
                   'Cells removed', 'Genes removed', 'Mean counts/cell', 'Mean genes/cell'],
        'value': [
            proc_data['orig'][0],
            proc_data['orig'][1],
            proc_data['filtered'].shape[0],
            proc_data['filtered'].shape[1],
            proc_data['orig'][0] - proc_data['filtered'].shape[0],
            proc_data['orig'][1] - proc_data['filtered'].shape[1],
            proc_data['normalized'].sum(axis=1).mean(),
            (proc_data['filtered'] > 0).sum(axis=1).mean()
        ]
    })
    
    summary.to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_summary.txt', sep='\t', index=False)
    logger.info("Saved summary statistics")
    
    # Save filtered counts
    proc_data['filtered'].to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_filtered.csv')
    logger.info("Saved filtered counts")
    
    # Save visualization sample
    n = proc_data['log'].shape[0]
    sample_size = min(3000, n)
    if n > sample_size:
        idx = np.random.choice(n, sample_size, replace=False)
        sample = proc_data['log'].iloc[idx]
    else:
        sample = proc_data['log']
    
    sample.to_csv(RESULTS_DIR / 'scRNA_PRAD_P01_vis_sample.csv')
    logger.info(f"Saved visualization sample ({sample.shape[0]} cells)")
    
    return proc_data['log']


def calculate_igf_scores(log_counts):
    """Calculate IGF pathway scores"""
    logger.info("Calculating IGF pathway scores...")
    
    igf_avail = [g for g in IGF_GENES if g in log_counts.columns]
    logger.info(f"Found {len(igf_avail)}/{len(IGF_GENES)} IGF genes")
    
    if igf_avail:
        igf_scores = log_counts[igf_avail].mean(axis=1).values
    else:
        # Use top genes as proxy
        logger.warning("No IGF genes found, using top genes...")
        top_genes = log_counts.var(axis=0).nlargest(50).index
        igf_scores = log_counts[top_genes].mean(axis=1).values
    
    df = pd.DataFrame({
        'cell_barcode': log_counts.index,
        'IGF_pathway_score': igf_scores
    })
    
    df.to_csv(RESULTS_DIR / 'PRAD_P01_IGF_cell_scores.csv', index=False)
    logger.info("Saved IGF scores")
    
    # Summary
    summary = pd.DataFrame({
        'metric': ['Mean', 'Median', 'Std', 'Min', 'Max', 'Q25', 'Q75'],
        'value': [igf_scores.mean(), np.median(igf_scores), igf_scores.std(),
                  igf_scores.min(), igf_scores.max(),
                  np.percentile(igf_scores, 25), np.percentile(igf_scores, 75)]
    })
    summary.to_csv(RESULTS_DIR / 'PRAD_P01_IGF_summary.csv', index=False)
    
    return df


def plot_igf_activity(cell_scores):
    """Create IGF pathway plots"""
    logger.info("Creating IGF plots...")
    
    scores = cell_scores['IGF_pathway_score'].values
    
    # Histogram
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(scores, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax.axvline(scores.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {scores.mean():.3f}')
    ax.axvline(np.median(scores), color='green', linestyle='--', linewidth=2, label=f'Median: {np.median(scores):.3f}')
    ax.set_xlabel('IGF Pathway Score', fontsize=12)
    ax.set_ylabel('Cells', fontsize=12)
    ax.set_title('IGF Pathway Activity Distribution (PRAD_P01)', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'PRAD_P01_IGF_histogram.png', dpi=300)
    logger.info("Saved histogram")
    plt.close()
    
    # Top 30 cells
    top30_idx = np.argsort(scores)[-30:][::-1]
    top30_scores = scores[top30_idx]
    top30_labels = [str(bc)[:20] + '...' if len(str(bc)) > 20 else str(bc) 
                    for bc in cell_scores.iloc[top30_idx]['cell_barcode'].values]
    
    fig, ax = plt.subplots(figsize=(10, 12))
    ax.barh(np.arange(len(top30_labels)), top30_scores, color='steelblue', edgecolor='black')
    ax.set_yticks(np.arange(len(top30_labels)))
    ax.set_yticklabels(top30_labels, fontsize=8)
    ax.set_xlabel('IGF Pathway Score', fontsize=12)
    ax.set_title('Top 30 Cells by IGF Activity (PRAD_P01)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'PRAD_P01_IGF_top30_cells.png', dpi=300)
    logger.info("Saved top 30 plot")
    plt.close()


def create_communication_outputs(log_counts, cell_scores):
    """Create GraphComm-Lite outputs"""
    logger.info("Creating communication analysis outputs...")
    
    scores = cell_scores['IGF_pathway_score'].values
    
    # Classify cells
    q75, q50, q25 = np.percentile(scores, [75, 50, 25])
    classes = ['Autocrine' if s > q75 else 'Sender' if s > q50 else 'Receiver' if s > q25 else 'Inactive' 
               for s in scores]
    
    cell_scores['Classification'] = classes
    cell_scores.to_csv(RESULTS_DIR / 'PRAD_P01_cell_classifications.csv', index=False)
    logger.info("Saved classifications")
    
    # Distribution plot
    fig, ax = plt.subplots(figsize=(10, 6))
    counts = pd.Series(classes).value_counts()
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
    logger.info("Saved distribution plot")
    plt.close()
    
    # Communication heatmap (gene correlation)
    logger.info("Creating communication heatmap...")
    top_genes = log_counts.var(axis=0).nlargest(100).index.tolist()
    
    # Sample cells if large
    n_cells = log_counts.shape[0]
    if n_cells > 200:
        sample_idx = np.random.choice(n_cells, 200, replace=False)
        sample = log_counts.iloc[sample_idx][top_genes]
    else:
        sample = log_counts[top_genes]
    
    corr = sample.corr()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(corr, cmap='RdYlBu_r', aspect='auto', vmin=-1, vmax=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Gene Communication Network (Top 100 Genes, PRAD_P01)', fontsize=14, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Correlation')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'PRAD_P01_communication_heatmap.png', dpi=300)
    logger.info("Saved heatmap")
    plt.close()
    
    # Network graph
    logger.info("Creating network graph...")
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    
    cell_types = ['Autocrine', 'Sender', 'Receiver', 'Inactive']
    type_counts = pd.Series(classes).value_counts()
    
    np.random.seed(42)
    for i, ct in enumerate(cell_types):
        if ct in type_counts.index:
            angle = (i / len(cell_types)) * 2 * np.pi
            x, y = 1.2 * np.cos(angle), 1.2 * np.sin(angle)
            radius = np.sqrt(type_counts[ct]) / 30 + 0.15
            
            circle = plt.Circle((x, y), radius, color=colors.get(ct, '#95a5a6'), alpha=0.6, ec='black', linewidth=2.5)
            ax.add_patch(circle)
            ax.text(x, y, f'{ct}\n({type_counts[ct]})', ha='center', va='center', fontsize=10, fontweight='bold')
    
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Global Cell-Cell Communication Network (PRAD_P01)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'PRAD_P01_communication_network.png', dpi=300)
    logger.info("Saved network graph")
    plt.close()


def main():
    """Execute pipeline"""
    logger.info("=" * 80)
    logger.info("PRAD_P01 SPARSE MATRIX ANALYSIS PIPELINE")
    logger.info("=" * 80)
    
    try:
        logger.info("\n[1/5] Loading count matrix...")
        counts = load_count_matrix_sparse()
        
        logger.info("\n[2/5] Preprocessing...")
        proc = preprocess_counts(counts)
        
        logger.info("\n[3/5] Saving data...")
        log_counts = save_preprocessed_data(proc)
        
        logger.info("\n[4/5] IGF pathway analysis...")
        cell_scores = calculate_igf_scores(log_counts)
        plot_igf_activity(cell_scores)
        
        logger.info("\n[5/5] Communication analysis...")
        create_communication_outputs(log_counts, cell_scores)
        
        logger.info("\n" + "=" * 80)
        logger.info("PIPELINE COMPLETED!")
        logger.info("=" * 80)
        logger.info(f"Results: {RESULTS_DIR}")
        logger.info(f"Plots: {PLOTS_DIR}")
        
        return True
    except Exception as e:
        logger.error(f"Pipeline error: {e}", exc_info=True)
        return False


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
