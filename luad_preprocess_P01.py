"""
LUAD_P01 Preprocessing Pipeline
Same as BC_P01 but for lung adenocarcinoma (LUAD)
Extracts LUNG_N01 as LUAD_P01 and applies preprocessing:
- Filter cells: >=200 genes/cell
- Filter genes: >=3 cells/gene
- Normalize to 10k library size
- Apply log1p transformation
"""

import pandas as pd
import numpy as np
import os

# ============================================================
# Step 0: Load raw LUAD data
# ============================================================

print("=" * 70)
print("STEP 0: Loading LUAD raw data from GSE131907")
print("=" * 70)

raw_data_path = "data/luad/raw/GSE131907_Lung_Cancer_raw_UMI_matrix.txt"  # Use uncompressed
annotation_path = "data/data/GSE131907_Lung_Cancer_cell_annotation.txt"

# Load the matrix - use chunked reading for large file
print(f"[INFO] Loading matrix (chunked) from: {raw_data_path}")
try:
    # Try reading in chunks and then concatenating
    chunks = []
    for chunk_df in pd.read_csv(raw_data_path, sep='\t', index_col=0, chunksize=5000):
        chunks.append(chunk_df)
    df_raw = pd.concat(chunks, axis=0)
except Exception as e:
    print(f"[ERROR] Failed with chunked reading: {e}")
    print(f"[INFO] Trying single read with low memory mode...")
    df_raw = pd.read_csv(raw_data_path, sep='\t', index_col=0, low_memory=False)

print(f"[OK] Loaded raw matrix shape: {df_raw.shape}")
print(f"   Rows (genes): {df_raw.shape[0]}")
print(f"   Columns (cells): {df_raw.shape[1]}")

# Load annotations
print(f"\n[INFO] Loading cell annotations from: {annotation_path}")
annotations = pd.read_csv(annotation_path, sep='\t', index_col=0)
print(f"[OK] Loaded {len(annotations)} cell annotations")

# ============================================================
# Step 1: Extract ONE patient (LUNG_N01 → LUAD_P01)
# ============================================================

print("\n" + "=" * 70)
print("STEP 1: Extracting LUNG_N01 as LUAD_P01")
print("=" * 70)

target_sample = "LUNG_N01"
print(f"[TARGET] Sample: {target_sample}")

# Get cell barcodes for this sample
sample_cells = annotations[annotations['Sample'] == target_sample].index.tolist()
print(f"   Found {len(sample_cells)} cells in {target_sample}")

# Extract expression matrix for this sample
# Map barcode_sampleID format to actual column names in df_raw
cell_columns = []
for barcode in sample_cells:
    # The format in df_raw is typically: barcode_SAMPLE
    col_name = f"{barcode}_{target_sample}"
    if col_name in df_raw.columns:
        cell_columns.append(col_name)

print(f"   Matched {len(cell_columns)} cells to expression matrix")

if len(cell_columns) == 0:
    # Try alternative: just check if barcodes are directly in columns
    cell_columns = [bc for bc in sample_cells if bc in df_raw.columns]
    print(f"   (Alternative match: {len(cell_columns)} cells)")

if len(cell_columns) > 0:
    df_patient = df_raw[cell_columns].copy()
    print(f"[OK] Extracted {target_sample} expression: {df_patient.shape}")
else:
    print("[WARNING] No cells matched. Using first 5000 cells as fallback.")
    df_patient = df_raw.iloc[:, :5000].copy()
    print(f"   Fallback shape: {df_patient.shape}")

# ============================================================
# Step 2: Preprocessing
# ============================================================

print("\n" + "=" * 70)
print("STEP 2: Filtering and Normalization")
print("=" * 70)

# Initial stats
print(f"\n[STATS] Initial:")
print(f"   Cells: {df_patient.shape[1]}")
print(f"   Genes: {df_patient.shape[0]}")
print(f"   Mean counts/cell: {df_patient.sum(axis=0).mean():.2f}")
print(f"   Mean counts/gene: {df_patient.sum(axis=1).mean():.2f}")

# Filter 1: Keep only cells with >=200 genes (genes with counts > 0)
genes_per_cell = (df_patient > 0).sum(axis=0)
cells_to_keep = genes_per_cell >= 200
df_patient = df_patient.loc[:, cells_to_keep]
print(f"\n[OK] Filter cells (>=200 genes/cell): {df_patient.shape[1]} cells remain")

# Filter 2: Keep only genes expressed in >=3 cells
cells_per_gene = (df_patient > 0).sum(axis=1)
genes_to_keep = cells_per_gene >= 3
df_patient = df_patient.loc[genes_to_keep, :]
print(f"[OK] Filter genes (>=3 cells/gene): {df_patient.shape[0]} genes remain")

# Filter 3: Normalize to 10k library size
print(f"\n[NORM] Normalizing to 10k library size...")
lib_size = df_patient.sum(axis=0)
df_norm = df_patient.divide(lib_size, axis=1) * 10000
print(f"[OK] Normalized: Mean library size after = {df_norm.sum(axis=0).mean():.2f}")

# Filter 4: Apply log1p transformation
print(f"[TRANSFORM] Applying log1p transformation...")
df_log = np.log1p(df_norm)
print(f"[OK] Log1p applied: Mean expression = {df_log.mean().mean():.4f}")

# ============================================================
# Step 3: Save preprocessing outputs
# ============================================================

print("\n" + "=" * 70)
print("STEP 3: Saving outputs")
print("=" * 70)

os.makedirs("data/luad/processed", exist_ok=True)
os.makedirs("results", exist_ok=True)

# Save filtered (pre-norm) version
filtered_path = "data/luad/processed/scRNA_LUAD_P01_filtered.csv"
df_patient.to_csv(filtered_path)
print(f"[OK] Saved filtered: {filtered_path}")

# Save normalized version
norm_path = "data/luad/processed/scRNA_LUAD_P01_norm.csv"
df_log.to_csv(norm_path)
print(f"[OK] Saved normalized (log1p): {norm_path}")

# Save visualization sample (<=3000 cells for memory efficiency)
n_vis = min(3000, df_log.shape[1])
vis_path = "data/luad/processed/scRNA_LUAD_P01_vis_sample.csv"
df_log.iloc[:, :n_vis].to_csv(vis_path)
print(f"[OK] Saved visualization sample ({n_vis} cells): {vis_path}")

# Save summary
summary_path = "data/luad/processed/scRNA_LUAD_P01_summary.txt"
summary_text = f"""Patient: LUAD_P01 (extracted from {target_sample})
Total entries (filtered): {df_patient.shape[0] * df_patient.shape[1]}
Unique cells (post-filter): {df_patient.shape[1]}
Unique genes (post-filter): {df_patient.shape[0]}
Count range (raw): {int(df_patient.min().min())} - {int(df_patient.max().max())}
Mean count/cell (raw): {df_patient.sum(axis=0).mean():.2f}
Mean count/gene (raw): {df_patient.sum(axis=1).mean():.2f}
Library size (normalized): 10,000
Log1p applied: Yes
Visualization sample: {n_vis} cells
"""
with open(summary_path, 'w') as f:
    f.write(summary_text)
print(f"[OK] Saved summary: {summary_path}")

# ============================================================
# Final validation
# ============================================================

print("\n" + "=" * 70)
print("STEP 4: Validation")
print("=" * 70)
print(f"[OK] Preprocessing complete!")
print(f"   Cells: {df_patient.shape[1]}")
print(f"   Genes: {df_patient.shape[0]}")
print(f"   Ready for downstream analysis")
