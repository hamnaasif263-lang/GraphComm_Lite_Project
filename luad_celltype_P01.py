"""
LUAD_P01 Cell-Type Annotation
Annotate cells as:
- Tumor epithelial cells
- Fibroblasts / CAFs
- Macrophages  
- T cells
- Endothelial cells
Using marker gene expression
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

print("=" * 70)
print("LUAD_P01: Cell-Type Annotation")
print("=" * 70)

# Load normalized data
norm_path = "data/luad/processed/scRNA_LUAD_P01_norm.parquet"
if not os.path.exists(norm_path):
    print(f"ERROR: {norm_path} not found!")
    print("Run luad_process_P01.py first")
    exit(1)

print(f"\nLoading normalized data from {norm_path}")
df_norm = pd.read_parquet(norm_path)
print(f"Loaded shape: {df_norm.shape}")

# ============================================================
# Define marker gene signatures for each cell type
# ============================================================

CELL_TYPE_MARKERS = {
    'Tumor_Epithelial': ['EPCAM', 'KRT7', 'KRT19', 'CDH1', 'TP63'],
    'Fibroblasts_CAF': ['COL1A1', 'COL1A2', 'FN1', 'VIM', 'THY1'],
    'Macrophages': ['CD14', 'CD68', 'CD163', 'MARCO', 'MSR1'],
    'T_Cells': ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A'],
    'Endothelial': ['PECAM1', 'CDH5', 'VWF', 'PLVAP', 'ENG']
}

print(f"\n[MARKERS] Cell type signatures:")
for ctype, genes in CELL_TYPE_MARKERS.items():
    print(f"  {ctype}: {genes}")

# ============================================================
# Compute cell-type scores
# ============================================================

cell_type_scores = {}

for cell_type, genes in CELL_TYPE_MARKERS.items():
    found_genes = [g for g in genes if g in df_norm.columns]
    
    if len(found_genes) > 0:
        # Mean expression of markers for this cell type
        score = df_norm[found_genes].mean(axis=0)
    else:
        # No markers found, use zero
        score = pd.Series(np.zeros(df_norm.shape[1]), index=df_norm.columns)
    
    cell_type_scores[cell_type] = score
    print(f"[{cell_type}] Found {len(found_genes)}/{len(genes)} markers")

# Convert to DataFrame
scores_df = pd.DataFrame(cell_type_scores)

# Assign each cell to the cell type with highest score
cell_assignments = scores_df.idxmax(axis=1)
cell_scores_max = scores_df.max(axis=1)

print(f"\nCell-type scores computed")
print(f"  Cells: {len(cell_assignments)}")

# ============================================================
# Save annotations
# ============================================================

os.makedirs("results", exist_ok=True)

# Save per-cell annotations
annotations_path = "results/LUAD_P01_celltype_annotations.csv"
annotations_df = pd.DataFrame({
    'Cell': df_norm.columns,
    'Cell_Type': cell_assignments.values,
    'Score': cell_scores_max.values
})
annotations_df.to_csv(annotations_path, index=False)
print(f"\n[OK] Saved annotations: {annotations_path}")

# ============================================================
# Generate summary and plot
# ============================================================

print(f"\n[SUMMARY] Cell-type distribution:")
counts = cell_assignments.value_counts()
for ctype in counts.index:
    pct = 100 * counts[ctype] / len(cell_assignments)
    print(f"  {ctype}: {counts[ctype]:6d} cells ({pct:5.1f}%)")

# Plot cell-type distribution
try:
    os.makedirs("plots", exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Bar plot
    counts.sort_values(ascending=False).plot(
        kind='bar',
        ax=ax,
        color='steelblue',
        edgecolor='black',
        alpha=0.7
    )
    
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Number of Cells')
    ax.set_title('LUAD_P01: Cell-Type Distribution')
    ax.grid(alpha=0.3, axis='y')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    dist_path = "plots/LUAD_P01_celltype_distribution.png"
    plt.savefig(dist_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n[PLOT] Saved distribution: {dist_path}")
except Exception as e:
    print(f"[WARNING] Failed to save plot: {e}")

# ============================================================
# Validation
# ============================================================

print("\n" + "=" * 70)
print("VALIDATION")
print("=" * 70)
print(f"[OK] Cell-type annotation complete!")
print(f"     Total cells: {len(cell_assignments)}")
print(f"     Cell types: {cell_assignments.nunique()}")
print(f"     Most abundant: {counts.idxmax()} ({counts.max()} cells)")
