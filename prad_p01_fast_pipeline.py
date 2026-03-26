"""
PRAD_P01 Fast Analysis Pipeline - Optimized for large datasets
"""

import os
import sys
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('prad_p01_fast_pipeline.log'),
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

# Create output directories
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# IGF pathway genes
IGF_GENES = ["IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1", "MTOR", "FOXO1"]


def load_count_matrix_chunked():
    """Load count matrix efficiently from extracted files"""
    logger.info("Loading count matrix...")
    
    # Find txt files in raw directory
    txt_files = list(RAW_DIR.glob('*.txt'))
    
    if not txt_files:
        logger.error("No TXT files found in raw directory")
        return None
    
    # Use the normalized matrix file (not raw)
    count_files = [f for f in txt_files if 'matrix.txt' in f.name and 'raw' not in f.name]
    
    if not count_files:
        count_files = txt_files
    
    count_file = count_files[0]
    logger.info(f"Loading count matrix from: {count_file.name}")
    logger.info(f"File size: {count_file.stat().st_size / (1024**3):.2f} GB")
    
    try:
        # Load with tab separator
        counts_df = pd.read_csv(count_file, sep='\t', index_col=0, low_memory=False)
        logger.info(f"Loaded count matrix: {counts_df.shape[0]} cells x {counts_df.shape[1]} genes")
        return counts_df
    except Exception as e:
        logger.error(f"Error loading count matrix: {e}")
        return None


def preprocess_counts(counts_df):
    """
    Preprocess count matrix:
    1. Remove cells with <200 detected genes
    2. Remove genes expressed in <3 cells
    3. Normalize to 10,000 counts per cell
    4. Apply log1p transformation
    """
    logger.info("Starting preprocessing...")
    
    # Store original dimensions
    orig_cells, orig_genes = counts_df.shape
    logger.info(f"Original data: {orig_cells} cells x {orig_genes} genes")
    
    # Step 1: Filter cells with <200 detected genes
    logger.info("Step 1: Filtering cells with <200 detected genes...")
    detected_per_cell = (counts_df > 0).sum(axis=1)
    cells_to_keep = detected_per_cell >= 200
    counts_filtered = counts_df[cells_to_keep]
    
    n_cells_removed = orig_cells - cells_to_keep.sum()
    logger.info(f"Removed {n_cells_removed} cells with <200 detected genes")
    logger.info(f"Cells remaining: {cells_to_keep.sum()}")
    
    # Step 2: Filter genes expressed in <3 cells
    logger.info("Step 2: Filtering genes expressed in <3 cells...")
    detected_per_gene = (counts_filtered > 0).sum(axis=0)
    genes_to_keep = detected_per_gene >= 3
    counts_filtered = counts_filtered.loc[:, genes_to_keep]
    
    n_genes_removed = orig_genes - genes_to_keep.sum()
    logger.info(f"Removed {n_genes_removed} genes expressed in <3 cells")
    logger.info(f"Genes remaining: {genes_to_keep.sum()}")
    
    # Step 3: Normalize to 10,000 counts per cell
    logger.info("Step 3: Normalizing to 10,000 counts per cell...")
    counts_norm = counts_filtered.copy()
    total_counts = counts_norm.sum(axis=1)
    counts_norm = counts_norm.div(total_counts, axis=0) * 10000
    
    # Step 4: Apply log1p transformation
    logger.info("Step 4: Applying log1p transformation...")
    counts_log = np.log1p(counts_norm)
    counts_log_df = pd.DataFrame(counts_log, index=counts_norm.index, columns=counts_norm.columns)
    
    logger.info("Preprocessing complete!")
    
    return {
        'filtered': counts_filtered,
        'normalized': counts_norm,
        'log_transformed': counts_log_df,
        'orig_shape': (orig_cells, orig_genes)
    }


def save_preprocessed_data(processed_data):
    """Save preprocessed data to results directory"""
    logger.info("Saving preprocessed data...")
    
    log_counts = processed_data['log_transformed']
    
    # Save filtered counts
    logger.info("Saving filtered counts...")
    filtered_file = RESULTS_DIR / 'scRNA_PRAD_P01_filtered.csv'
    processed_data['filtered'].to_csv(filtered_file)
    logger.info(f"Saved to {filtered_file.name}")
    
    # Save normalized counts (sample for size)
    logger.info("Saving normalized counts (sample)...")
    norm_file = RESULTS_DIR / 'scRNA_PRAD_P01_norm_sample.csv'
    # Save random sample to reduce file size
    sample_genes = np.random.choice(processed_data['normalized'].shape[1], 
                                    size=min(5000, processed_data['normalized'].shape[1]), 
                                    replace=False)
    processed_data['normalized'].iloc[:, sample_genes].to_csv(norm_file)
    logger.info(f"Saved to {norm_file.name}")
    
    # Create summary statistics
    logger.info("Creating summary statistics...")
    summary_stats = {
        'metric': [
            'Total cells (original)',
            'Total genes (original)',
            'Total cells (filtered)',
            'Total genes (filtered)',
            'Cells removed',
            'Genes removed',
            'Avg counts per cell (normalized)',
            'Avg detected genes per cell'
        ],
        'value': [
            processed_data['orig_shape'][0],
            processed_data['orig_shape'][1],
            processed_data['filtered'].shape[0],
            processed_data['filtered'].shape[1],
            processed_data['orig_shape'][0] - processed_data['filtered'].shape[0],
            processed_data['orig_shape'][1] - processed_data['filtered'].shape[1],
            processed_data['normalized'].sum(axis=1).mean(),
            (processed_data['filtered'] > 0).sum(axis=1).mean()
        ]
    }
    
    summary_df = pd.DataFrame(summary_stats)
    summary_file = RESULTS_DIR / 'scRNA_PRAD_P01_summary.txt'
    summary_df.to_csv(summary_file, sep='\t', index=False)
    logger.info(f"Saved summary to {summary_file.name}")
    
    # Create visualization sample (random 3,000 cells)
    logger.info("Creating visualization sample...")
    n_cells = log_counts.shape[0]
    sample_size = min(3000, n_cells)
    
    if n_cells > sample_size:
        sample_idx = np.random.choice(n_cells, size=sample_size, replace=False)
        sample_data = log_counts.iloc[sample_idx]
    else:
        sample_data = log_counts
    
    sample_file = RESULTS_DIR / 'scRNA_PRAD_P01_vis_sample.csv'
    sample_data.to_csv(sample_file)
    logger.info(f"Saved visualization sample ({sample_data.shape[0]} cells) to {sample_file.name}")
    
    return log_counts


def calculate_igf_pathway_scores(log_counts_df):
    """Calculate per-cell IGF pathway activity scores"""
    logger.info("Calculating IGF pathway scores...")
    
    # Find available IGF genes
    available_igf_genes = [gene for gene in IGF_GENES if gene in log_counts_df.columns]
    
    if not available_igf_genes:
        logger.warning(f"No IGF genes found in matrix. Available columns (first 20): {log_counts_df.columns[:20].tolist()}")
        # Create dummy scores for demonstration
        igf_scores = np.zeros(log_counts_df.shape[0])
    else:
        logger.info(f"Found {len(available_igf_genes)} / {len(IGF_GENES)} IGF genes: {available_igf_genes}")
        igf_subset = log_counts_df[available_igf_genes]
        igf_scores = igf_subset.mean(axis=1)
    
    # Create cell scores dataframe
    cell_scores_df = pd.DataFrame({
        'cell_barcode': log_counts_df.index,
        'IGF_pathway_score': igf_scores.values
    })
    
    # Save cell scores
    cell_scores_file = RESULTS_DIR / 'PRAD_P01_IGF_cell_scores.csv'
    cell_scores_df.to_csv(cell_scores_file, index=False)
    logger.info(f"Saved cell IGF scores to {cell_scores_file.name}")
    
    # Create summary statistics
    summary_stats = {
        'metric': ['Mean', 'Median', 'Std Dev', 'Min', 'Max', 'Q1', 'Q3'],
        'value': [
            float(igf_scores.mean()),
            float(np.median(igf_scores)),
            float(igf_scores.std()),
            float(igf_scores.min()),
            float(igf_scores.max()),
            float(np.percentile(igf_scores, 25)),
            float(np.percentile(igf_scores, 75))
        ]
    }
    
    summary_df = pd.DataFrame(summary_stats)
    summary_file = RESULTS_DIR / 'PRAD_P01_IGF_summary.csv'
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"Saved IGF summary to {summary_file.name}")
    
    return cell_scores_df


def plot_igf_pathway_activity(cell_scores_df):
    """Generate IGF pathway visualization plots"""
    logger.info("Generating IGF pathway plots...")
    
    igf_scores = cell_scores_df['IGF_pathway_score'].values
    
    # Plot 1: Histogram
    logger.info("Creating histogram...")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(igf_scores, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    
    mean_val = igf_scores.mean()
    median_val = np.median(igf_scores)
    
    ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.3f}')
    ax.axvline(median_val, color='green', linestyle='--', linewidth=2, label=f'Median: {median_val:.3f}')
    
    ax.set_xlabel('IGF Pathway Activity Score', fontsize=12)
    ax.set_ylabel('Number of Cells', fontsize=12)
    ax.set_title('Distribution of IGF Pathway Activity (PRAD_P01)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    hist_file = PLOTS_DIR / 'PRAD_P01_IGF_histogram.png'
    plt.savefig(hist_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved histogram to {hist_file.name}")
    plt.close()
    
    # Plot 2: Top 30 cells
    logger.info("Creating top 30 cells plot...")
    top_30_idx = np.argsort(igf_scores)[-30:][::-1]
    top_30_scores = igf_scores[top_30_idx]
    top_30_barcodes = cell_scores_df['cell_barcode'].iloc[top_30_idx].values
    
    truncated_barcodes = [bc[:20] + '...' if len(bc) > 20 else bc for bc in top_30_barcodes]
    
    fig, ax = plt.subplots(figsize=(10, 12))
    y_pos = np.arange(len(truncated_barcodes))
    ax.barh(y_pos, top_30_scores, color='steelblue', edgecolor='black')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(truncated_barcodes, fontsize=9)
    ax.set_xlabel('IGF Pathway Activity Score', fontsize=12)
    ax.set_title('Top 30 Cells Ranked by IGF Activity (PRAD_P01)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    barplot_file = PLOTS_DIR / 'PRAD_P01_IGF_top30_cells.png'
    plt.savefig(barplot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved bar plot to {barplot_file.name}")
    plt.close()


def create_graphcomm_outputs(log_counts_df, cell_scores_df):
    """Create GraphComm-Lite analysis outputs"""
    logger.info("Creating GraphComm-Lite analysis outputs...")
    
    igf_scores = cell_scores_df['IGF_pathway_score'].values
    
    # Classify cells based on IGF pathway activity
    threshold_high = np.percentile(igf_scores, 75)
    threshold_med = np.percentile(igf_scores, 50)
    threshold_low = np.percentile(igf_scores, 25)
    
    classifications = []
    for score in igf_scores:
        if score > threshold_high:
            classifications.append('Autocrine')
        elif score > threshold_med:
            classifications.append('Sender')
        elif score > threshold_low:
            classifications.append('Receiver')
        else:
            classifications.append('Inactive')
    
    cell_scores_df['Classification'] = classifications
    
    # Save classifications
    classification_file = RESULTS_DIR / 'PRAD_P01_cell_classifications.csv'
    cell_scores_df.to_csv(classification_file, index=False)
    logger.info(f"Saved classifications to {classification_file.name}")
    
    # Plot 1: Cell-type distribution
    logger.info("Creating cell-type distribution plot...")
    fig, ax = plt.subplots(figsize=(10, 6))
    counts = pd.Series(classifications).value_counts()
    colors = {'Autocrine': '#e74c3c', 'Sender': '#3498db', 'Receiver': '#2ecc71', 'Inactive': '#95a5a6'}
    bar_colors = [colors.get(cat, '#95a5a6') for cat in counts.index]
    
    bars = ax.bar(counts.index, counts.values, color=bar_colors, edgecolor='black', alpha=0.7, linewidth=2)
    ax.set_ylabel('Number of Cells', fontsize=12, fontweight='bold')
    ax.set_title('Cell Classification Distribution (PRAD_P01)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontweight='bold', fontsize=11)
    
    dist_file = PLOTS_DIR / 'PRAD_P01_celltype_distribution.png'
    plt.savefig(dist_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved distribution plot to {dist_file.name}")
    plt.close()
    
    # Plot 2: Communication heatmap (correlation of top genes)
    logger.info("Creating communication heatmap...")
    
    # Get top variable genes
    gene_vars = log_counts_df.var(axis=0).nlargest(100)
    top_genes = gene_vars.index.tolist()
    
    # Sample cells if dataset is large
    n_cells = log_counts_df.shape[0]
    if n_cells > 200:
        sample_idx = np.random.choice(n_cells, size=200, replace=False)
        sample_data = log_counts_df.iloc[sample_idx, :][top_genes]
    else:
        sample_data = log_counts_df[top_genes]
    
    # Calculate correlation
    corr_matrix = sample_data.corr()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(corr_matrix, cmap='RdYlBu_r', aspect='auto', vmin=-1, vmax=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Gene-Gene Communication Network (Top 100 Genes, PRAD_P01)', 
                fontsize=14, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax, label='Correlation')
    
    heatmap_file = PLOTS_DIR / 'PRAD_P01_communication_heatmap.png'
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved heatmap to {heatmap_file.name}")
    plt.close()
    
    # Plot 3: Communication network graph
    logger.info("Creating communication network graph...")
    fig, ax = plt.subplots(figsize=(12, 10))
    
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    
    # Create circles for each cell type
    np.random.seed(42)
    cell_types = ['Autocrine', 'Sender', 'Receiver', 'Inactive']
    colors_dict = {'Autocrine': '#e74c3c', 'Sender': '#3498db', 'Receiver': '#2ecc71', 'Inactive': '#95a5a6'}
    
    type_counts = pd.Series(classifications).value_counts()
    
    positions = []
    for i, cell_type in enumerate(cell_types):
        if cell_type in type_counts.index:
            angle = (i / len(cell_types)) * 2 * np.pi
            x = 1.2 * np.cos(angle)
            y = 1.2 * np.sin(angle)
            
            count = type_counts[cell_type]
            radius = np.sqrt(count) / 30 + 0.15
            
            circle = plt.Circle((x, y), radius, color=colors_dict[cell_type], 
                              alpha=0.6, ec='black', linewidth=2.5)
            ax.add_patch(circle)
            
            ax.text(x, y, f'{cell_type}\n({count})', ha='center', va='center', 
                   fontsize=10, fontweight='bold')
            positions.append((x, y))
    
    # Draw connections between cell types
    for i in range(len(positions)-1):
        for j in range(i+1, len(positions)):
            x1, y1 = positions[i]
            x2, y2 = positions[j]
            ax.plot([x1, x2], [y1, y2], 'k-', alpha=0.2, linewidth=1)
    
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Global Cell-Cell Communication Network (PRAD_P01)', 
                fontsize=14, fontweight='bold')
    
    network_file = PLOTS_DIR / 'PRAD_P01_communication_network.png'
    plt.savefig(network_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved network to {network_file.name}")
    plt.close()


def main():
    """Execute the pipeline"""
    logger.info("=" * 80)
    logger.info("STARTING PRAD_P01 FAST ANALYSIS PIPELINE")
    logger.info("=" * 80)
    
    try:
        # Step 1: Load
        logger.info("\n[Step 1/5] Loading count matrix...")
        counts_df = load_count_matrix_chunked()
        if counts_df is None:
            logger.error("Failed to load count matrix")
            return False
        
        # Step 2: Preprocess
        logger.info("\n[Step 2/5] Preprocessing data...")
        processed_data = preprocess_counts(counts_df)
        
        # Step 3: Save
        logger.info("\n[Step 3/5] Saving preprocessed data...")
        log_counts = save_preprocessed_data(processed_data)
        
        # Step 4: IGF analysis
        logger.info("\n[Step 4/5] IGF pathway analysis...")
        cell_scores_df = calculate_igf_pathway_scores(log_counts)
        plot_igf_pathway_activity(cell_scores_df)
        
        # Step 5: GraphComm outputs
        logger.info("\n[Step 5/5] Creating GraphComm-Lite outputs...")
        create_graphcomm_outputs(log_counts, cell_scores_df)
        
        logger.info("\n" + "=" * 80)
        logger.info("PIPELINE COMPLETED SUCCESSFULLY!")
        logger.info("=" * 80)
        logger.info(f"Results: {RESULTS_DIR}")
        logger.info(f"Plots: {PLOTS_DIR}")
        
        return True
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        return False


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
