import pandas as pd
import os

print("Testing different file formats:")

# Try LUAD parquet
print("\n" + "="*60)
print("LUAD_P01 Parquet:")
try:
    df = pd.read_parquet("data/luad/processed/scRNA_LUAD_P01_norm.parquet")
    print(f"Shape: {df.shape}")
    print(f"Index (first 5): {df.index[:5].tolist()}")
    print(f"Columns (first 5): {list(df.columns[:5])}")
    # Check if it has our target genes
    target_genes = ["IGF1R", "INSR", "IRS1", "PIK3CA", "AKT1", "MTOR"]
    found = [g for g in target_genes if g in df.columns]
    print(f"Target genes found: {found}")
except Exception as e:
    print(f"Error: {e}")

# Try BC CSV as-is
print("\n" + "="*60)
print("BC_P01 CSV (Long format - raw):")
try:
    df = pd.read_csv("data/breast_cancer/raw/scRNA_BC_P01_norm.csv", nrows=1000)
    print(f"Shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    print("First 3 rows:")
    print(df.head(3))
except Exception as e:
    print(f"Error: {e}")

# Try BC filtered CSV
print("\n" + "="*60)
print("BC_P01 Filtered CSV (Long format - filtered):")
try:
    df = pd.read_csv("data/breast_cancer/raw/scRNA_BC_P01_filtered.csv", nrows=1000)
    print(f"Shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    print("First 3 rows:")
    print(df.head(3))
except Exception as e:
    print(f"Error: {e}")

# Try PRAD CSV
print("\n" + "="*60)
print("PRAD_P01 CSV (should be matrix):")
try:
    df = pd.read_csv("graphcomm/results/prostate/scRNA_PRAD_P01_norm.csv", nrows=10, index_col=0)
    print(f"Shape: {df.shape}")
    print(f"Index (first 3): {df.index[:3].tolist()}")
    print(f"Columns (first 5): {list(df.columns[:5])}")
    target_genes = ["IGF1R", "INSR", "IRS1", "PIK3CA", "AKT1", "MTOR"]
    found = [g for g in target_genes if g in df.columns]
    print(f"Target genes found: {found}")
except Exception as e:
    print(f"Error: {e}")
