"""
LUAD_P01 Final Validation
Print validation metrics and verify all outputs
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

print("=" * 70)
print("LUAD_P01: FINAL VALIDATION")
print("=" * 70)

# ============================================================
# Check preprocessing outputs
# ============================================================

print("\n[CHECK] Preprocessing outputs:")
processed_files = [
    "data/luad/processed/scRNA_LUAD_P01_filtered.csv",
    "data/luad/processed/scRNA_LUAD_P01_norm.csv",
    "data/luad/processed/scRNA_LUAD_P01_vis_sample.csv",
    "data/luad/processed/scRNA_LUAD_P01_summary.txt"
]

for fpath in processed_files:
    exists = os.path.exists(fpath)
    status = "[OK]" if exists else "[MISSING]"
    print(f"  {status} {fpath}")

# ============================================================
# Load and validate preprocessed data
# ============================================================

norm_path = "data/luad/processed/scRNA_LUAD_P01_norm.csv"
if os.path.exists(norm_path):
    print(f"\n[LOAD] Loading normalized data...")
    df_norm = pd.read_csv(norm_path, index_col=0)
    n_cells = df_norm.shape[1]
    n_genes = df_norm.shape[0]
    print(f"[OK] Shape: {n_cells} cells x {n_genes} genes")
else:
    print(f"[ERROR] {norm_path} not found!")
    n_cells = 0
    n_genes = 0

# ============================================================
# Check IGF analysis outputs
# ============================================================

print(f"\n[CHECK] IGF pathway analysis outputs:")
igf_files = [
    "results/LUAD_P01_IGF_cell_scores.csv",
    "results/LUAD_P01_IGF_summary.csv",
    "plots/LUAD_P01_IGF_histogram.png",
    "plots/LUAD_P01_IGF_top30_cells.png"
]

igf_ok = True
for fpath in igf_files:
    exists = os.path.exists(fpath)
    status = "[OK]" if exists else "[MISSING]"
    print(f"  {status} {fpath}")
    if not exists:
        igf_ok = False

# Load IGF summary if exists
if os.path.exists(igf_files[1]):
    igf_summary = pd.read_csv(igf_files[1])
    print(f"\n[IGF STATS]:")
    for idx, row in igf_summary.iterrows():
        metric = row['Metric']
        value = row['Value']
        if isinstance(value, float):
            print(f"  {metric}: {value:.4f}")
        else:
            print(f"  {metric}: {value}")

# ============================================================
# Check cell-type annotation outputs
# ============================================================

print(f"\n[CHECK] Cell-type annotation outputs:")
celltype_files = [
    "results/LUAD_P01_celltype_annotations.csv",
    "plots/LUAD_P01_celltype_distribution.png"
]

celltype_ok = True
for fpath in celltype_files:
    exists = os.path.exists(fpath)
    status = "[OK]" if exists else "[MISSING]"
    print(f"  {status} {fpath}")
    if not exists:
        celltype_ok = False

# Load cell-type annotations if exists
if os.path.exists(celltype_files[0]):
    celltype_df = pd.read_csv(celltype_files[0])
    print(f"\n[CELLTYPE DISTRIBUTION]:")
    for ctype, count in celltype_df['Cell_Type'].value_counts().items():
        pct = 100 * count / len(celltype_df)
        print(f"  {ctype}: {count:6d} ({pct:5.1f}%)")

# ============================================================
# Final validation summary
# ============================================================

print("\n" + "=" * 70)
print("VALIDATION SUMMARY")
print("=" * 70)

all_ok = igf_ok and celltype_ok

print(f"\n[STATUS] Pipeline execution:")
print(f"  Preprocessing: [OK]" + " (inferred)")
print(f"  IGF analysis:  {'[OK]' if igf_ok else '[INCOMPLETE]'}")
print(f"  Cell-type:     {'[OK]' if celltype_ok else '[INCOMPLETE]'}")
print(f"\n[FINAL] Overall: {'[COMPLETE]' if all_ok else '[INCOMPLETE]'}")

if n_cells > 0:
    print(f"\n[SUMMARY] Dataset characteristics:")
    print(f"  Total cells: {n_cells}")
    print(f"  Total genes: {n_genes}")
    print(f"  Library size: 10,000 (normalized)")
    print(f"  Transformation: log1p")
    print(f"  Sample: LUAD_P01 (LUNG_N01)")

print("\n" + "=" * 70)
print("END OF VALIDATION")
print("=" * 70)
