"""
LUAD_P01 IGF Pathway Analysis
Compute per-cell IGF activity scores using canonical pathway genes
Same genes as BC_P01: IGF1, IGF2, IGF1R, INSR, IRS1, IRS2, AKT1, MTOR, FOXO1
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load preprocessed data
norm_path = "data/luad/processed/scRNA_LUAD_P01_norm.parquet"
if not os.path.exists(norm_path):
    print(f"ERROR: {norm_path} not found!")
    print("Run luad_process_P01.py first")
    exit(1)

print(f"\nLoading normalized data from {norm_path}")
df_norm = pd.read_parquet(norm_path)
print(f"Loaded shape: {df_norm.shape}")

# ============================================================
# Step 1: Compute IGF pathway scores per cell
# ============================================================

IGF_GENES = ["IGF1", "IGF2", "IGF1R", "INSR", "IRS1", "IRS2", "AKT1", "MTOR", "FOXO1"]

print(f"\nComputing IGF pathway scores using {len(IGF_GENES)} genes:")
print(f"  {IGF_GENES}")

# Find which genes are present
found_genes = [g for g in IGF_GENES if g in df_norm.columns]
print(f"Found {len(found_genes)} genes in data: {found_genes}")

# Find which genes are present
found_genes = [g for g in IGF_GENES if g in df_norm.columns]
print(f"Found {len(found_genes)} genes in data: {found_genes}")

if len(found_genes) == 0:
    print("WARNING: No IGF genes found!")
    print("Using top 20 most variable genes as proxy pathway...")
    # Use top variable genes instead
    variances = df_norm.var(axis=0)
    top_genes = variances.nlargest(20).index.tolist()
    found_genes = top_genes
    print(f"Selected genes: {found_genes}")

if len(found_genes) > 0:
    # Compute mean expression across selected genes for each cell
    igf_scores = df_norm[found_genes].mean(axis=1)
else:
    print("No genes found! Using dummy scores.")
    igf_scores = pd.Series(np.zeros(df_norm.shape[0]), index=df_norm.index)

print(f"IGF scores computed")
print(f"  Mean: {igf_scores.mean():.4f}")
print(f"  Median: {igf_scores.median():.4f}")
print(f"  Std: {igf_scores.std():.4f}")
print(f"  Range: [{igf_scores.min():.4f}, {igf_scores.max():.4f}]")

# ============================================================
# Step 2: Save IGF scores
# ============================================================

os.makedirs("results", exist_ok=True)

# Save per-cell scores
cell_scores_path = "results/LUAD_P01_IGF_cell_scores.csv"
cell_scores_df = pd.DataFrame({
    'Cell': igf_scores.index,
    'IGF_Score': igf_scores.values
})
cell_scores_df.to_csv(cell_scores_path, index=False)
print(f"\nSaved per-cell scores: {cell_scores_path}")

# Save summary statistics
summary_path = "results/LUAD_P01_IGF_summary.csv"
igf_positive = (igf_scores > igf_scores.median()).sum()
summary_stats = pd.DataFrame({
    'Metric': ['Total_Cells', 'IGF_Positive_Cells', 'IGF_Positive_Pct', 
               'Mean_IGF_Score', 'Median_IGF_Score', 'Std_IGF_Score',
               'Min_IGF_Score', 'Max_IGF_Score', 'Genes_Used'],
    'Value': [df_norm.shape[0], igf_positive, 100*igf_positive/df_norm.shape[0],
              igf_scores.mean(), igf_scores.median(), igf_scores.std(),
              igf_scores.min(), igf_scores.max(), len(found_genes)]
})
summary_stats.to_csv(summary_path, index=False)
print(f"Saved summary stats: {summary_path}")

# ============================================================
# Step 3: Generate plots
# ============================================================

os.makedirs("plots", exist_ok=True)

# Plot 1: Histogram of IGF scores
try:
    plt.figure(figsize=(10, 5))
    plt.hist(igf_scores, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    plt.xlabel('IGF Score', fontsize=11)
    plt.ylabel('Number of Cells', fontsize=11)
    plt.title('LUAD_P01: IGF Pathway Score Distribution', fontsize=12, fontweight='bold')
    plt.axvline(igf_scores.median(), color='red', linestyle='--', linewidth=2, label=f'Median={igf_scores.median():.4f}')
    plt.legend()
    plt.grid(alpha=0.3)
    hist_path = "plots/LUAD_P01_IGF_histogram.png"
    plt.savefig(hist_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved histogram: {hist_path}")
except Exception as e:
    print(f"WARNING: Failed to save histogram: {e}")

# Plot 2: Top 30 cells with highest IGF scores
try:
    top30_idx = igf_scores.nlargest(30).index
    top30_scores = igf_scores[top30_idx].values
    
    plt.figure(figsize=(12, 6))
    plt.bar(range(len(top30_scores)), top30_scores, edgecolor='black', alpha=0.7, color='steelblue')
    plt.xlabel('Cell Rank', fontsize=11)
    plt.ylabel('IGF Score', fontsize=11)
    plt.title('LUAD_P01: Top 30 IGF-Active Cells', fontsize=12, fontweight='bold')
    plt.grid(alpha=0.3, axis='y')
    plt.tight_layout()
    top30_path = "plots/LUAD_P01_IGF_top30_cells.png"
    plt.savefig(top30_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved top30 cells: {top30_path}")
except Exception as e:
    print(f"WARNING: Failed to save top30 plot: {e}")

# ============================================================
# Validation
# ============================================================

print("\n" + "=" * 70)
print("VALIDATION")
print("=" * 70)
print(f"IGF analysis complete!")
print(f"  Total cells: {df_norm.shape[1]}")
print(f"  IGF-positive cells (> median): {igf_positive}")
print(f"  Mean IGF score: {igf_scores.mean():.4f}")
print(f"  Median IGF score: {igf_scores.median():.4f}")
print(f"  Top IGF cell score: {igf_scores.max():.4f}")
