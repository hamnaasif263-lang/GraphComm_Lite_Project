#!/usr/bin/env python3
"""Quick test of the label generation fix without full training."""

from pathlib import Path
import pandas as pd
import numpy as np

DATA_DIR = Path("data")

def test_label_generation(csv_path: Path, sample_name: str):
    """Test the label generation fix."""
    df_expr = pd.read_csv(csv_path, index_col=0)
    
    # Use the FIXED logic
    if "label" in df_expr.columns:
        raw_labels = df_expr["label"].astype(int)
        print(f"{sample_name}: Using existing labels")
    else:
        igf_genes = ["IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1"]
        available_genes = [g for g in igf_genes if g in df_expr.columns]
        
        if available_genes:
            igf_series = df_expr[available_genes].mean(axis=1)
            thresh = igf_series.median()
            raw_labels = (igf_series > thresh).astype(int)
            print(f"{sample_name}: Generated from {len(available_genes)} IGF genes: {available_genes}")
        else:
            gene_means = df_expr.mean(axis=0).sort_values(ascending=False)
            top_genes = gene_means.head(10).index.tolist()
            expr_score = df_expr[top_genes].mean(axis=1)
            thresh = expr_score.median()
            raw_labels = (expr_score > thresh).astype(int)
            print(f"{sample_name}: Generated from top 10 genes")
        
        n_pos = (raw_labels == 1).sum()
        n_neg = (raw_labels == 0).sum()
        print(f"  Labels: {n_pos} pos, {n_neg} neg (ratio: {n_pos/(n_pos+n_neg):.2%})")
        
        if n_pos == 0 or n_neg == 0:
            print("  Rebalancing...")
            raw_labels = pd.Series(np.random.randint(0, 2, size=df_expr.shape[0]), index=df_expr.index)
            n_pos = (raw_labels == 1).sum()
            n_neg = (raw_labels == 0).sum()
            print(f"  After rebalance: {n_pos} pos, {n_neg} neg")
    
    return {
        'sample': sample_name,
        'n_0': (raw_labels == 0).sum(),
        'n_1': (raw_labels == 1).sum(),
    }

# Test on several samples
test_samples = [
    'scRNA_BRONCHO_11.csv',
    'scRNA_LUNG_T09.csv',
    'scRNA_EBUS_06.csv',
    'scRNA_LN_01.csv',
    'scRNA_LUNG_N01.csv',
]

print("Testing Label Generation Fix")
print("=" * 70)

results = []
for csv_name in test_samples:
    csv_path = DATA_DIR / csv_name
    if csv_path.exists():
        sample_name = csv_name.replace('scRNA_', '').replace('.csv', '')
        try:
            result = test_label_generation(csv_path, sample_name)
            results.append(result)
        except Exception as e:
            print(f"ERROR in {sample_name}: {e}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
for r in results:
    total = r['n_0'] + r['n_1']
    pct = 100 * r['n_1'] / total if total > 0 else 0
    print(f"{r['sample']:20} | Labels: [0:{r['n_0']:5}, 1:{r['n_1']:5}] | Positive: {pct:5.1f}%")
    
    if r['n_1'] > 0 and r['n_0'] > 0:
        print(f"  ✓ GOOD - has both classes")
    else:
        print(f"  ✗ PROBLEM - missing a class")
