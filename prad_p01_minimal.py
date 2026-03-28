"""
PRAD_P01 Minimal Processing Pipeline
Same as LUAD_P01 but for prostate cancer
Uses scRNA data from graphcomm/data/prostate/raw/
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from pathlib import Path

print("=" * 70)
print("PRAD_P01 MINIMAL ANALYSIS PIPELINE")
print("=" * 70)

# ============================================================
# STEP 0: Load PRAD data
# ============================================================

print("\nSTEP 0: Loading PRAD data")
print("-" * 70)

data_dir = Path("graphcomm/data/prostate/raw")
data_files = list(data_dir.glob("*.csv"))

if not data_files:
    print(f"ERROR: No CSV files found in {data_dir}")
    exit(1)

data_path = data_files[0]
print(f"Loading: {data_path}")

df_raw = pd.read_csv(data_path, index_col=0)
print(f"Shape: {df_raw.shape}")
print(f"Cells (rows): {df_raw.shape[0]}, Genes (columns): {df_raw.shape[1]}")

# ============================================================
# STEP 1: Filter cells and genes
# ============================================================

print("\nSTEP 1: Filtering cells and genes")
print("-" * 70)

initial_cells = df_raw.shape[0]
initial_genes = df_raw.shape[1]

# Filter: cells with >= 200 detected genes
genes_per_cell = (df_raw > 0).sum(axis=1)
cells_keep = genes_per_cell >= 200
df_filt = df_raw.loc[cells_keep, :]
print(f"Cell filter (>=200 genes): {df_filt.shape[0]} cells (removed {initial_cells - df_filt.shape[0]})")

# Filter: genes expressed in >= 3 cells
cells_per_gene = (df_filt > 0).sum(axis=0)
genes_keep = cells_per_gene >= 3
df_filt = df_filt.loc[:, genes_keep]
print(f"Gene filter (>=3 cells): {df_filt.shape[1]} genes (removed {initial_genes - df_filt.shape[1]})")
print(f"Final shape: {df_filt.shape}")

# ============================================================
# STEP 2: Normalize to 10k library size
# ============================================================

print("\nSTEP 2: Normalizing (10k library)")
print("-" * 70)

lib_size = df_filt.sum(axis=1)
df_norm = df_filt.divide(lib_size, axis=0) * 10000
print(f"Mean library size: {df_norm.sum(axis=1).mean():.1f}")
print(f"Library range: {df_norm.sum(axis=1).min():.0f} - {df_norm.sum(axis=1).max():.0f}")

# ============================================================
# STEP 3: Log1p transformation
# ============================================================

print("\nSTEP 3: Log1p transformation")
print("-" * 70)

df_log = np.log1p(df_norm)
print(f"Mean log expression: {df_log.mean().mean():.4f}")
print(f"Log range: {df_log.min().min():.4f} - {df_log.max().max():.4f}")

# ============================================================
# STEP 4: Save preprocessing outputs
# ============================================================

print("\nSTEP 4: Saving outputs")
print("-" * 70)

results_dir = Path("graphcomm/results/prostate")
results_dir.mkdir(parents=True, exist_ok=True)

# Save filtered
filt_path = results_dir / "scRNA_PRAD_P01_filtered.csv"
df_filt.to_csv(filt_path)
print(f"✓ Saved: {filt_path}")

# Save normalized
norm_path = results_dir / "scRNA_PRAD_P01_norm.csv"
df_log.to_csv(norm_path)
print(f"✓ Saved: {norm_path}")

# Save summary
summary_path = results_dir / "scRNA_PRAD_P01_summary.txt"
summary_text = f"""PRAD_P01 PREPROCESSING SUMMARY
{'='*60}

Dataset: Prostate Cancer scRNA-seq (PRAD_P01)

Cells (post-filter): {df_filt.shape[0]:,}
Genes (post-filter): {df_filt.shape[1]:,}
Initial cells: {initial_cells:,}
Initial genes: {initial_genes:,}

Filtering:
  - Cell filter: >=200 detected genes
  - Gene filter: >=3 cells per gene

Library size statistics:
  - Mean: {lib_size.mean():.1f}
  - Median: {lib_size.median():.1f}
  - Min: {lib_size.min():.0f}
  - Max: {lib_size.max():.0f}

Normalization: 10,000 counts per cell
Transformation: log1p
"""
with open(summary_path, 'w') as f:
    f.write(summary_text)
print(f"✓ Saved: {summary_path}")

# ============================================================
# STEP 5: IGF Pathway Analysis
# ============================================================

print("\nSTEP 5: IGF Pathway Analysis")
print("-" * 70)

# IGF pathway genes
igf_genes = ["IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1", "MTOR", "FOXO1"]

# Find IGF genes in the data
igf_found = [g for g in igf_genes if g in df_log.columns]
print(f"IGF genes found: {len(igf_found)}/{len(igf_genes)}")
if igf_found:
    print(f"  {igf_found}")

if len(igf_found) > 0:
    # Calculate IGF pathway activity (mean expression)
    igf_expr = df_log[igf_found]
    igf_scores = igf_expr.mean(axis=1)
    
    print(f"\nIGF Activity Statistics:")
    print(f"  Mean: {igf_scores.mean():.4f}")
    print(f"  Median: {igf_scores.median():.4f}")
    print(f"  Std: {igf_scores.std():.4f}")
    print(f"  Range: {igf_scores.min():.4f} - {igf_scores.max():.4f}")
    
    # Save IGF cell scores
    igf_scores_path = results_dir / "PRAD_P01_IGF_cell_scores.csv"
    igf_df = pd.DataFrame({
        'cell': df_log.index,
        'igf_score': igf_scores.values
    })
    igf_df.to_csv(igf_scores_path, index=False)
    print(f"\n✓ Saved: {igf_scores_path}")
    
    # Save IGF summary
    igf_summary_path = results_dir / "PRAD_P01_IGF_summary.csv"
    igf_summary = pd.DataFrame({
        'metric': ['n_genes', 'n_cells', 'mean', 'median', 'std', 'min', 'max'],
        'value': [len(igf_found), len(igf_scores), igf_scores.mean(), 
                  igf_scores.median(), igf_scores.std(), igf_scores.min(), igf_scores.max()]
    })
    igf_summary.to_csv(igf_summary_path, index=False)
    print(f"✓ Saved: {igf_summary_path}")
    
    # ============================================================
    # STEP 6: Generate plots
    # ============================================================
    
    print("\nSTEP 6: Generating plots")
    print("-" * 70)
    
    plots_dir = Path("graphcomm/plots/prostate")
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Plot 1: Histogram of IGF activity
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    ax.hist(igf_scores, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(igf_scores.mean(), color='red', linestyle='--', linewidth=2, 
               label=f'Mean: {igf_scores.mean():.3f}')
    ax.axvline(igf_scores.median(), color='green', linestyle='--', linewidth=2,
               label=f'Median: {igf_scores.median():.3f}')
    ax.set_xlabel('IGF Pathway Activity Score')
    ax.set_ylabel('Number of Cells')
    ax.set_title('IGF Pathway Activity Distribution (PRAD_P01)')
    ax.legend()
    ax.grid(alpha=0.3)
    
    hist_path = plots_dir / "PRAD_P01_IGF_histogram.png"
    plt.savefig(hist_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {hist_path}")
    
    # Plot 2: Top 30 cells
    top_30_idx = np.argsort(igf_scores.values)[-30:][::-1]
    top_30_scores = igf_scores.values[top_30_idx]
    top_30_cells = [str(c)[:20] for c in igf_scores.index[top_30_idx]]
    
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    y_pos = np.arange(len(top_30_scores))
    ax.barh(y_pos, top_30_scores, color='steelblue', edgecolor='black')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_30_cells, fontsize=9)
    ax.set_xlabel('IGF Pathway Activity Score')
    ax.set_title('Top 30 Cells by IGF Pathway Activity (PRAD_P01)')
    ax.invert_yaxis()
    ax.grid(alpha=0.3, axis='x')
    
    top30_path = plots_dir / "PRAD_P01_IGF_top30_cells.png"
    plt.savefig(top30_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {top30_path}")

# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("PRAD_P01 PIPELINE COMPLETE")
print("=" * 70)
print(f"\nOutput directory: {results_dir}")
print(f"Plot directory: {plots_dir}")
print(f"\nGenerated files:")
print(f"  - scRNA_PRAD_P01_filtered.csv")
print(f"  - scRNA_PRAD_P01_norm.csv")
print(f"  - scRNA_PRAD_P01_summary.txt")
print(f"  - PRAD_P01_IGF_cell_scores.csv")
print(f"  - PRAD_P01_IGF_summary.csv")
print(f"  - PRAD_P01_IGF_histogram.png")
print(f"  - PRAD_P01_IGF_top30_cells.png")
print("\n✅ Analysis complete!")
