"""
LUAD_P01 Complete Processing Pipeline
Same as BC_P01 but for lung adenocarcinoma (LUAD)
Uses scRNA_LUNG_N01.csv as LUAD_P01
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ============================================================
# STEP 0: Load LUAD_N01 data
# ============================================================

print("=" * 70)
print("STEP 0: Loading LUAD data")
print("=" * 70)

data_path = "data/scRNA_LUNG_N01.csv"
print(f"Loading: {data_path}")
df_raw = pd.read_csv(data_path, index_col=0)
print(f"Shape: {df_raw.shape}")
print(f"Cells (rows): {df_raw.shape[0]}, Genes (columns): {df_raw.shape[1]}")
print(f"Sample genes: {df_raw.columns[:10].tolist()}")

# ============================================================
# STEP 1: Filter cells and genes
# ============================================================

print("\n" + "=" * 70)
print("STEP 1: Filtering cells and genes")
print("=" * 70)

initial_cells = df_raw.shape[0]
initial_genes = df_raw.shape[1]

# Filter: cells with >= 10 genes (genes with counts > 0)
genes_per_cell = (df_raw > 0).sum(axis=1)
cells_keep = genes_per_cell >= 10
df_filt = df_raw.loc[cells_keep, :]
print(f"After cell filter (>=10 genes): {df_filt.shape[0]} cells (was {initial_cells})")

# Filter: genes expressed in >= 2 cells
cells_per_gene = (df_filt > 0).sum(axis=0)
genes_keep = cells_per_gene >= 2
df_filt = df_filt.loc[:, genes_keep]
print(f"After gene filter (>=2 cells): {df_filt.shape[1]} genes (was {initial_genes})")

# ============================================================
# STEP 2: Normalize to 10k library size
# ============================================================

print("\n" + "=" * 70)
print("STEP 2: Normalizing (10k library)")
print("=" * 70)

lib_size = df_filt.sum(axis=1)  # Sum across genes (row sums)
df_norm = df_filt.divide(lib_size, axis=0) * 10000
print(f"Mean library size: {df_norm.sum(axis=1).mean():.1f}")

# ============================================================
# STEP 3: Log1p transformation
# ============================================================

print("\n" + "=" * 70)
print("STEP 3: Log1p transformation")
print("=" * 70)

df_log = np.log1p(df_norm)
print(f"Transformed to log1p scale")
print(f"Mean expression: {df_log.mean().mean():.4f}")

# ============================================================
# STEP 4: Save preprocessing outputs
# ============================================================

print("\n" + "=" * 70)
print("STEP 4: Saving outputs")
print("=" * 70)

os.makedirs("data/luad/processed", exist_ok=True)

# Save filtered
filt_path = "data/luad/processed/scRNA_LUAD_P01_filtered.parquet"
df_filt.to_parquet(filt_path)
print(f"Saved: {filt_path}")

# Save normalized
norm_path = "data/luad/processed/scRNA_LUAD_P01_norm.parquet"
df_log.to_parquet(norm_path)
print(f"Saved: {norm_path}")

# Save visualization sample (<=3000 cells)
n_vis = min(3000, df_log.shape[0])
vis_path = "data/luad/processed/scRNA_LUAD_P01_vis_sample.parquet"
df_log.iloc[:n_vis, :].to_parquet(vis_path)
print(f"Saved: {vis_path} ({n_vis} cells)")

# Save summary
summary_path = "data/luad/processed/scRNA_LUAD_P01_summary.txt"
summary_text = f"""Patient: LUAD_P01 (LUNG_N01)
Total cells (post-filter): {df_filt.shape[0]}
Total genes (post-filter): {df_filt.shape[1]}
Initial cells: {initial_cells}
Initial genes: {initial_genes}
Count range (raw): {int(df_filt.min().min())} - {int(df_filt.max().max())}
Mean counts/cell: {df_filt.sum(axis=1).mean():.2f}
Mean counts/gene: {df_filt.sum(axis=0).mean():.2f}
Normalization: 10k library
Transformation: log1p
Visualization sample: {n_vis} cells
"""
with open(summary_path, 'w') as f:
    f.write(summary_text)
print(f"Saved: {summary_path}")

print("\n" + "=" * 70)
print("PREPROCESSING COMPLETE")
print("=" * 70)
print(f"Ready for IGF pathway analysis")
