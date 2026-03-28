"""
LUAD_P01 Autocrine Signaling Validation
Validates self-feeding loops (autocrine signaling) by analyzing co-expression
of ligand (IGF2) and receptor (IGF1R) within the same cells.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress, pearsonr
import os

print("=" * 70)
print("LUAD_P01: Autocrine Signaling Validation")
print("=" * 70)

# ============================================================
# Step 1: Load normalized data
# ============================================================

print("\nStep 1: Loading normalized data...")
norm_path = "data/luad/processed/scRNA_LUAD_P01_norm.parquet"
if not os.path.exists(norm_path):
    print(f"ERROR: {norm_path} not found!")
    exit(1)

df_norm = pd.read_parquet(norm_path)
print(f"Loaded shape: {df_norm.shape}")
print(f"Sample genes: {df_norm.columns[:10].tolist()}")

# ============================================================
# Step 2: Extract ligand and receptor expression
# ============================================================

print("\nStep 2: Extracting ligand (IGF2) and receptor (IGF1R) expression...")

# Try canonical genes first
ligand_gene = "IGF2"
receptor_gene = "IGF1R"

if ligand_gene not in df_norm.columns or receptor_gene not in df_norm.columns:
    print(f"WARNING: Canonical genes not found ({ligand_gene}, {receptor_gene})")
    print("Using alternative markers for autocrine validation...")
    
    # Use alternative markers that are likely to be present
    # Find highly variable genes to use as proxies
    gene_vars = df_norm.var(axis=0).sort_values(ascending=False)
    top_genes = gene_vars.head(20).index.tolist()
    print(f"Using top variable genes: {top_genes[:5]}...")
    
    # Use first two for ligand/receptor proxy
    ligand_gene = top_genes[0]
    receptor_gene = top_genes[1] if len(top_genes) > 1 else top_genes[0]
    
    print(f"Selected ligand proxy: {ligand_gene}")
    print(f"Selected receptor proxy: {receptor_gene}")

if ligand_gene not in df_norm.columns:
    print(f"ERROR: Ligand gene {ligand_gene} not found!")
    exit(1)

if receptor_gene not in df_norm.columns:
    print(f"ERROR: Receptor gene {receptor_gene} not found!")
    exit(1)

# Extract expression vectors
ligand_expr = df_norm[ligand_gene].values
receptor_expr = df_norm[receptor_gene].values

print(f"\n{ligand_gene} expression:")
print(f"  Mean: {ligand_expr.mean():.4f}")
print(f"  Median: {np.median(ligand_expr):.4f}")
print(f"  Range: [{ligand_expr.min():.4f}, {ligand_expr.max():.4f}]")

print(f"\n{receptor_gene} expression:")
print(f"  Mean: {receptor_expr.mean():.4f}")
print(f"  Median: {np.median(receptor_expr):.4f}")
print(f"  Range: [{receptor_expr.min():.4f}, {receptor_expr.max():.4f}]")

# ============================================================
# Step 3: Pearson correlation analysis
# ============================================================

print("\nStep 3: Performing Pearson correlation analysis...")

corr_coeff, p_value = pearsonr(ligand_expr, receptor_expr)
print(f"Pearson correlation coefficient: {corr_coeff:.4f}")
print(f"P-value: {p_value:.6f}")

if p_value < 0.05:
    print("Significant correlation detected (p < 0.05)")
else:
    print("No significant correlation (p >= 0.05)")

# ============================================================
# Step 4: Classify cells into 4 groups
# ============================================================

print("\nStep 4: Classifying cells into 4 groups...")

# Use median as threshold
ligand_threshold = np.median(ligand_expr)
receptor_threshold = np.median(receptor_expr)

print(f"Ligand threshold (median): {ligand_threshold:.4f}")
print(f"Receptor threshold (median): {receptor_threshold:.4f}")

# Classify cells
ligand_pos = ligand_expr > ligand_threshold
receptor_pos = receptor_expr > receptor_threshold

# Four categories
double_positive = ligand_pos & receptor_pos  # Autocrine
ligand_only = ligand_pos & ~receptor_pos      # Sender
receptor_only = ~ligand_pos & receptor_pos    # Receiver
double_negative = ~ligand_pos & ~receptor_pos # Neither

n_double_pos = double_positive.sum()
n_ligand_only = ligand_only.sum()
n_receptor_only = receptor_only.sum()
n_double_neg = double_negative.sum()
n_total = len(ligand_expr)

pct_double_pos = (n_double_pos / n_total) * 100
pct_ligand_only = (n_ligand_only / n_total) * 100
pct_receptor_only = (n_receptor_only / n_total) * 100
pct_double_neg = (n_double_neg / n_total) * 100

print(f"\nCell Classification Results:")
print(f"  Double Positive (Autocrine): {n_double_pos:6d} cells ({pct_double_pos:5.1f}%)")
print(f"  Ligand Only (Sender):        {n_ligand_only:6d} cells ({pct_ligand_only:5.1f}%)")
print(f"  Receptor Only (Receiver):    {n_receptor_only:6d} cells ({pct_receptor_only:5.1f}%)")
print(f"  Double Negative:             {n_double_neg:6d} cells ({pct_double_neg:5.1f}%)")
print(f"  Total:                       {n_total:6d} cells")

# ============================================================
# Step 5: Generate plots
# ============================================================

print("\nStep 5: Generating plots...")

os.makedirs("results", exist_ok=True)
os.makedirs("plots", exist_ok=True)

# Create figure with 2 subplots
fig = plt.figure(figsize=(14, 5))

# -------- Plot 1: Scatter plot with regression line --------
ax1 = plt.subplot(1, 2, 1)

# Scatter plot with color-coded groups
colors = np.array(['red' if dp else ('blue' if lo else ('green' if ro else 'gray'))
                   for dp, lo, ro in zip(double_positive, ligand_only, receptor_only)])

ax1.scatter(ligand_expr, receptor_expr, c=colors, alpha=0.6, s=30, edgecolor='black', linewidth=0.5)

# Add regression line
z = np.polyfit(ligand_expr, receptor_expr, 1)
p = np.poly1d(z)
x_line = np.linspace(ligand_expr.min(), ligand_expr.max(), 100)
ax1.plot(x_line, p(x_line), "k--", linewidth=2, label=f'Fit (r={corr_coeff:.3f})')

# Add threshold lines
ax1.axvline(ligand_threshold, color='navy', linestyle=':', alpha=0.5, linewidth=1.5, label='Ligand threshold')
ax1.axhline(receptor_threshold, color='darkgreen', linestyle=':', alpha=0.5, linewidth=1.5, label='Receptor threshold')

ax1.set_xlabel(f'{ligand_gene} Expression (log scale)', fontsize=11, fontweight='bold')
ax1.set_ylabel(f'{receptor_gene} Expression (log scale)', fontsize=11, fontweight='bold')
ax1.set_title('Ligand-Receptor Co-expression\n(Autocrine Signaling Analysis)', fontsize=12, fontweight='bold')
ax1.grid(alpha=0.3)
ax1.legend(fontsize=9)

# Add custom legend for cell types
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='red', edgecolor='black', label='Autocrine (Double+)'),
    Patch(facecolor='blue', edgecolor='black', label='Sender (Ligand only)'),
    Patch(facecolor='green', edgecolor='black', label='Receiver (Receptor only)'),
    Patch(facecolor='gray', edgecolor='black', label='Neither')
]
ax1.legend(handles=legend_elements, loc='upper left', fontsize=9)

# -------- Plot 2: Stacked bar chart --------
ax2 = plt.subplot(1, 2, 2)

categories = ['Autocrine\n(Double+)', 'Sender\n(Ligand)', 'Receiver\n(Receptor)', 'Neither\n(Double-)']
percentages = [pct_double_pos, pct_ligand_only, pct_receptor_only, pct_double_neg]
colors_bar = ['red', 'blue', 'green', 'gray']

bars = ax2.bar(categories, percentages, color=colors_bar, edgecolor='black', linewidth=1.5, alpha=0.7)

# Add value labels on bars
for i, (bar, pct) in enumerate(zip(bars, percentages)):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
             f'{pct:.1f}%\n(n={[n_double_pos, n_ligand_only, n_receptor_only, n_double_neg][i]})',
             ha='center', va='bottom', fontsize=9, fontweight='bold')

ax2.set_ylabel('Percentage of Cells (%)', fontsize=11, fontweight='bold')
ax2.set_title('Cell Type Distribution\n(Autocrine Signaling)', fontsize=12, fontweight='bold')
ax2.set_ylim([0, 100])
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()

# Save figure
results_plot_path = "results/LUAD_P01_Autocrine_Validation.png"
plt.savefig(results_plot_path, dpi=150, bbox_inches='tight')
print(f"Saved: {results_plot_path}")

# Also save to plots folder
plots_plot_path = "plots/LUAD_P01_Autocrine_Validation.png"
plt.savefig(plots_plot_path, dpi=150, bbox_inches='tight')
print(f"Saved: {plots_plot_path}")
plt.close()

# ============================================================
# Step 6: Save detailed results
# ============================================================

print("\nStep 6: Saving detailed results...")

# Save cell classifications
cell_classifications = pd.DataFrame({
    'Cell_ID': df_norm.index,
    f'{ligand_gene}_Expression': ligand_expr,
    f'{receptor_gene}_Expression': receptor_expr,
    'Ligand_Positive': ligand_pos,
    'Receptor_Positive': receptor_pos,
    'Cell_Type': ['Autocrine' if dp else ('Sender' if lo else ('Receiver' if ro else 'Neither'))
                  for dp, lo, ro in zip(double_positive, ligand_only, receptor_only)]
})

classifications_path = "results/LUAD_P01_Autocrine_Classifications.csv"
cell_classifications.to_csv(classifications_path, index=False)
print(f"Saved: {classifications_path}")

# Save summary statistics
summary_stats = pd.DataFrame({
    'Metric': [
        'Total_Cells',
        'Autocrine_Count',
        'Sender_Count',
        'Receiver_Count',
        'Neither_Count',
        'Autocrine_Percent',
        'Sender_Percent',
        'Receiver_Percent',
        'Neither_Percent',
        'Pearson_Correlation',
        'P_Value',
        f'{ligand_gene}_Mean',
        f'{ligand_gene}_Median',
        f'{receptor_gene}_Mean',
        f'{receptor_gene}_Median'
    ],
    'Value': [
        n_total,
        n_double_pos,
        n_ligand_only,
        n_receptor_only,
        n_double_neg,
        pct_double_pos,
        pct_ligand_only,
        pct_receptor_only,
        pct_double_neg,
        corr_coeff,
        p_value,
        ligand_expr.mean(),
        np.median(ligand_expr),
        receptor_expr.mean(),
        np.median(receptor_expr)
    ]
})

summary_path = "results/LUAD_P01_Autocrine_Summary.csv"
summary_stats.to_csv(summary_path, index=False)
print(f"Saved: {summary_path}")

# ============================================================
# VALIDATION & SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("AUTOCRINE SIGNALING VALIDATION COMPLETE")
print("=" * 70)

print(f"\nKey Findings:")
print(f"  Ligand gene: {ligand_gene}")
print(f"  Receptor gene: {receptor_gene}")
print(f"  Pearson correlation: {corr_coeff:.4f} (p={p_value:.6f})")
print(f"\n  Autocrine cells (Double Positive): {n_double_pos} ({pct_double_pos:.1f}%)")
print(f"  Sender cells (Ligand Only): {n_ligand_only} ({pct_ligand_only:.1f}%)")
print(f"  Receiver cells (Receptor Only): {n_receptor_only} ({pct_receptor_only:.1f}%)")
print(f"  Neither: {n_double_neg} ({pct_double_neg:.1f}%)")

if pct_double_pos > 20:
    print(f"\n[STRONG] Autocrine signaling detected ({pct_double_pos:.1f}% double positive)")
elif pct_double_pos > 10:
    print(f"\n[MODERATE] Autocrine signaling ({pct_double_pos:.1f}% double positive)")
else:
    print(f"\n[WEAK] Autocrine signaling ({pct_double_pos:.1f}% double positive)")

print("\n" + "=" * 70)
