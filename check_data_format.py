import pandas as pd
import os

# Check BC_P01
print("="*60)
print("BC_P01 Data Check:")
path = "data/breast_cancer/raw/scRNA_BC_P01_norm.csv"
if os.path.exists(path):
    df = pd.read_csv(path, nrows=3, index_col=0)
    print(f"Shape (rows=cells, cols=genes): {df.shape}")
    print(f"Index (first 5): {df.index[:5]}")
    print(f"Columns (first 5): {list(df.columns[:5])}")
else:
    print(f"File not found: {path}")

# Check LUAD_P01
print("\n" + "="*60)
print("LUAD_P01 Data Check:")
path = "data/luad/processed/scRNA_LUAD_P01_norm.csv"
if os.path.exists(path):
    df = pd.read_csv(path, nrows=3, index_col=0)
    print(f"Shape (rows=cells, cols=genes): {df.shape}")
    print(f"Index (first 5): {df.index[:5]}")
    print(f"Columns (first 5): {list(df.columns[:5])}")
else:
    print(f"File not found: {path}")

# Check PRAD_P01
print("\n" + "="*60)
print("PRAD_P01 Data Check:")
path = "graphcomm/results/prostate/scRNA_PRAD_P01_norm.csv"
if os.path.exists(path):
    df = pd.read_csv(path, nrows=3, index_col=0)
    print(f"Shape (rows=cells, cols=genes): {df.shape}")
    print(f"Index (first 5): {df.index[:5]}")
    print(f"Columns (first 5): {list(df.columns[:5])}")
else:
    print(f"File not found: {path}")
