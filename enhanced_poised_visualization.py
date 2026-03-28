"""
Enhanced Multi-Patient Poised State Validation Visualization
Creates publication-quality plots for the Poised State analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from pathlib import Path

# Set professional style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# Load data
OUTPUT_DIR = Path('graphcomm/results/breast')
PLOTS_DIR = Path('graphcomm/plots/breast')

df_results = pd.read_csv(OUTPUT_DIR / 'Poised_State_Statistics.csv')
df_summary = pd.read_csv(OUTPUT_DIR / 'Poised_State_Summary.csv')

# Extract key metrics
mean_poised = float(df_summary[df_summary['Metric'] == 'Mean_Poised_Percentage']['Value'].values[0])
median_poised = float(df_summary[df_summary['Metric'] == 'Median_Poised_Percentage']['Value'].values[0])
std_poised = float(df_summary[df_summary['Metric'] == 'Std_Dev_Poised_Percentage']['Value'].values[0])

print("Creating enhanced visualization...")

# Create figure with professional layout
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

# ============================================================================
# PLOT 1: MAIN POISED PERCENTAGE WITH ANNOTATIONS
# ============================================================================

ax1 = fig.add_subplot(gs[0, :])

colors_poised = ['#E74C3C', '#3498DB', '#2ECC71']
x_pos = np.arange(len(df_results))
bars = ax1.bar(x_pos, df_results['Percentage_Poised'], 
               color=colors_poised, alpha=0.8, edgecolor='black', linewidth=2.5, width=0.6)

# Add mean and median lines
ax1.axhline(y=mean_poised, color='#E74C3C', linestyle='--', linewidth=3, 
            label=f'Mean: {mean_poised:.2f}% (σ={std_poised:.2f}%)', zorder=5)
ax1.axhline(y=median_poised, color='#8E44AD', linestyle='-.', linewidth=2.5, 
            label=f'Median: {median_poised:.2f}%', zorder=4)

# Customize appearance
ax1.set_xticks(x_pos)
ax1.set_xticklabels(df_results['Patient_ID'], fontsize=13, fontweight='bold')
ax1.set_ylabel('Poised Cells (%)', fontsize=14, fontweight='bold')
ax1.set_title('Multi-Patient Poised State Validation - Breast Cancer Cohort\n' + 
              'MYC+ & CDKN1B+ Double Positive Cells', 
              fontsize=15, fontweight='bold', pad=20)

# Add value labels on bars
for i, (bar, val, total) in enumerate(zip(bars, df_results['Percentage_Poised'], df_results['Total_Cells'])):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 0.3,
            f'{val:.2f}%\n(n={int(df_results.iloc[i]["Poised_Cells"])})',
            ha='center', va='bottom', fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.8, edgecolor='gray'))

# Add grid and legend
ax1.grid(axis='y', alpha=0.3, linestyle=':', linewidth=1.5)
ax1.legend(fontsize=12, loc='upper right', framealpha=0.95, edgecolor='black', fancybox=True)
ax1.set_ylim(0, max(df_results['Percentage_Poised']) + 5)

# ============================================================================
# PLOT 2: STACKED BAR - CELL POPULATIONS
# ============================================================================

ax2 = fig.add_subplot(gs[1, 0])

poised = df_results['Poised_Cells'].values
non_poised = (df_results['Total_Cells'] - df_results['Poised_Cells']).values
myc_only = (df_results['MYC_Positive'] - df_results['Poised_Cells']).values
cdkn1b_only = (df_results['CDKN1B_Positive'] - df_results['Poised_Cells']).values
neither = df_results['Total_Cells'] - poised - myc_only - cdkn1b_only

x_pos = np.arange(len(df_results))
width = 0.6

# Create stacked bars
p1 = ax2.bar(x_pos, poised, width, label='Poised (MYC+ & CDKN1B+)', 
             color='#E74C3C', alpha=0.85, edgecolor='black', linewidth=1.5)
p2 = ax2.bar(x_pos, myc_only, width, bottom=poised, label='MYC+ Only',
             color='#F39C12', alpha=0.85, edgecolor='black', linewidth=1.5)
p3 = ax2.bar(x_pos, cdkn1b_only, width, bottom=poised+myc_only, label='CDKN1B+ Only',
             color='#3498DB', alpha=0.85, edgecolor='black', linewidth=1.5)
p4 = ax2.bar(x_pos, neither, width, bottom=poised+myc_only+cdkn1b_only, label='Neither',
             color='#95A5A6', alpha=0.7, edgecolor='black', linewidth=1.5)

ax2.set_xticks(x_pos)
ax2.set_xticklabels(df_results['Patient_ID'], fontsize=12, fontweight='bold')
ax2.set_ylabel('Cell Count', fontsize=13, fontweight='bold')
ax2.set_title('Cell Population Breakdown', fontsize=13, fontweight='bold', pad=15)
ax2.legend(fontsize=10, loc='upper right', framealpha=0.95, edgecolor='black')
ax2.grid(axis='y', alpha=0.3, linestyle=':', linewidth=1)

# ============================================================================
# PLOT 3: MARKER PREVALENCE HEATMAP
# ============================================================================

ax3 = fig.add_subplot(gs[1, 1])

marker_data = np.array([
    [df_results.iloc[i]['MYC_Positive'] / df_results.iloc[i]['Total_Cells'] * 100 for i in range(len(df_results))],
    [df_results.iloc[i]['CDKN1B_Positive'] / df_results.iloc[i]['Total_Cells'] * 100 for i in range(len(df_results))],
    [df_results.iloc[i]['Percentage_Poised'] for i in range(len(df_results))]
])

im = ax3.imshow(marker_data, cmap='RdYlGn', aspect='auto', vmin=0, vmax=60)

ax3.set_xticks(np.arange(len(df_results)))
ax3.set_yticks(np.arange(3))
ax3.set_xticklabels(df_results['Patient_ID'], fontsize=12, fontweight='bold')
ax3.set_yticklabels(['MYC+', 'CDKN1B+', 'Poised\n(Both+)'], fontsize=11, fontweight='bold')
ax3.set_title('Marker Expression Heatmap (%)', fontsize=13, fontweight='bold', pad=15)

# Add text annotations
for i in range(3):
    for j in range(len(df_results)):
        value = marker_data[i, j]
        text = ax3.text(j, i, f'{value:.1f}%', ha="center", va="center",
                       color="black" if value > 30 else "black", 
                       fontsize=10, fontweight='bold')

# Add colorbar
cbar = plt.colorbar(im, ax=ax3, orientation='vertical', pad=0.02)
cbar.set_label('Expression (%)', fontsize=11, fontweight='bold')

# ============================================================================
# PLOT 4: DISTRIBUTION ANALYSIS
# ============================================================================

ax4 = fig.add_subplot(gs[2, 0])

# Violin plot
parts = ax4.violinplot([df_results['Percentage_Poised'].values], positions=[1], 
                        widths=0.5, showmeans=True, showmedians=True)

# Customize violin plot colors
for pc in parts['bodies']:
    pc.set_facecolor('#E74C3C')
    pc.set_alpha(0.7)
    pc.set_edgecolor('black')
    pc.set_linewidth(1.5)

for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians', 'cmeans'):
    if partname in parts:
        vp = parts[partname]
        vp.set_edgecolor('black')
        vp.set_linewidth(2)

# Add scatter points
y_pos = df_results['Percentage_Poised'].values
x_jitter = np.random.normal(1, 0.04, len(y_pos))
ax4.scatter(x_jitter, y_pos, alpha=0.6, s=200, color='#2C3E50', edgecolor='black', linewidth=1.5, zorder=3)

# Add annotations
ax4.axhline(y=mean_poised, color='#E74C3C', linestyle='--', linewidth=2, alpha=0.7, label='Mean')
ax4.axhline(y=median_poised, color='#8E44AD', linestyle='-.', linewidth=2, alpha=0.7, label='Median')

ax4.set_xlim(0.5, 1.5)
ax4.set_xticks([])
ax4.set_ylabel('Poised Percentage (%)', fontsize=13, fontweight='bold')
ax4.set_title('Distribution of Poised Cells', fontsize=13, fontweight='bold', pad=15)
ax4.legend(fontsize=11, loc='upper right')
ax4.grid(axis='y', alpha=0.3, linestyle=':', linewidth=1)

# Add statistics box
stats_text = f"n={len(df_results)}\nMean: {mean_poised:.2f}%\nMedian: {median_poised:.2f}%\nRange: {df_results['Percentage_Poised'].min():.2f}-{df_results['Percentage_Poised'].max():.2f}%"
ax4.text(0.98, 0.05, stats_text, transform=ax4.transAxes, fontsize=10,
        verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9, edgecolor='black', linewidth=1.5),
        fontfamily='monospace')

# ============================================================================
# PLOT 5: SUMMARY STATISTICS TABLE
# ============================================================================

ax5 = fig.add_subplot(gs[2, 1])
ax5.axis('off')

# Create summary table
summary_data = []
for i, row in df_results.iterrows():
    summary_data.append([
        row['Patient_ID'],
        f"{int(row['Total_Cells'])}",
        f"{int(row['Poised_Cells'])}",
        f"{row['Percentage_Poised']:.2f}%"
    ])

summary_data.append(['', '', '', ''])
summary_data.append([
    'COHORT',
    f"{int(df_results['Total_Cells'].sum())}",
    f"{int(df_results['Poised_Cells'].sum())}",
    f"{mean_poised:.2f}%"
])

table = ax5.table(cellText=summary_data,
                 colLabels=['Patient', 'Total Cells', 'Poised Cells', '% Poised'],
                 cellLoc='center',
                 loc='center',
                 colWidths=[0.25, 0.25, 0.25, 0.25])

table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1, 2.5)

# Style header
for i in range(4):
    table[(0, i)].set_facecolor('#34495E')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Style data rows
for i in range(1, len(summary_data) + 1):
    for j in range(4):
        if i == len(summary_data):  # Cohort row
            table[(i, j)].set_facecolor('#E8DAEF')
            table[(i, j)].set_text_props(weight='bold')
        else:
            table[(i, j)].set_facecolor('#ECF0F1' if i % 2 == 0 else 'white')

ax5.text(0.5, 0.95, 'Summary Statistics', transform=ax5.transAxes,
        fontsize=13, fontweight='bold', ha='center', va='top')

# ============================================================================
# ADD OVERALL TITLE AND CONCLUSION BOX
# ============================================================================

fig.suptitle('Multi-Patient Poised State Validation Analysis\nBreast Cancer Single-Cell RNA-seq',
            fontsize=16, fontweight='bold', y=0.98)

# Add conclusion box at bottom
conclusion_text = (f"✓ Poised Phenotype Detected in {(df_results['Percentage_Poised'] > 0).sum()}/{len(df_results)} patients | "
                  f"Mean Prevalence: {mean_poised:.2f}% ± {std_poised:.2f}% | "
                  f"Definition: MYC > 0.5 AND CDKN1B > 0.5")

fig.text(0.5, 0.01, conclusion_text, ha='center', fontsize=11, 
        bbox=dict(boxstyle='round,pad=0.8', facecolor='#D5F4E6', alpha=0.9, 
                 edgecolor='#27AE60', linewidth=2),
        fontweight='bold')

# Save figure
plot_path = PLOTS_DIR / 'Multi_Patient_Poised_Validation.png'
plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"✓ Enhanced visualization saved: {plot_path}")

plt.close()

# Create additional detailed comparison plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Comparison with reference ranges
ax = axes[0]
patient_names = df_results['Patient_ID'].values
poised_pct = df_results['Percentage_Poised'].values

# Create bars with gradient coloring
colors_gradient = plt.cm.RdYlGn(poised_pct / 30)  # Scale to 0-30% range
bars = ax.barh(patient_names, poised_pct, color=colors_gradient, edgecolor='black', linewidth=2)

# Add reference zones
ax.axvline(x=mean_poised, color='#E74C3C', linestyle='--', linewidth=2.5, label=f'Cohort Mean: {mean_poised:.1f}%', zorder=5)
ax.fill_betweenx([-0.5, len(patient_names)-0.5], mean_poised-std_poised, mean_poised+std_poised, 
                 alpha=0.2, color='#E74C3C', label=f'±1 SD: {std_poised:.1f}%')

ax.set_xlabel('Poised Cells (%)', fontsize=12, fontweight='bold')
ax.set_title('Patient-Level Poised State Comparison', fontsize=13, fontweight='bold')
ax.legend(fontsize=11, loc='lower right')
ax.grid(axis='x', alpha=0.3)

# Add percentage labels
for i, (bar, val) in enumerate(zip(bars, poised_pct)):
    ax.text(val + 0.3, bar.get_y() + bar.get_height()/2, f'{val:.2f}%',
           va='center', fontsize=11, fontweight='bold')

# Plot 2: Marker co-expression network
ax = axes[1]

marker_names = ['MYC\nOnly', 'CDKN1B\nOnly', 'Poised\n(Both+)', 'Neither']
counts = []

for i in range(len(df_results)):
    total = df_results.iloc[i]['Total_Cells']
    poised = df_results.iloc[i]['Poised_Cells']
    myc_only = df_results.iloc[i]['MYC_Positive'] - poised
    cdkn1b_only = df_results.iloc[i]['CDKN1B_Positive'] - poised
    neither = total - poised - myc_only - cdkn1b_only
    
    counts.append([myc_only, cdkn1b_only, poised, neither])

counts = np.array(counts).T
x_pos = np.arange(len(df_results))
width = 0.2

colors_markers = ['#F39C12', '#3498DB', '#E74C3C', '#95A5A6']
for i, (count, name, color) in enumerate(zip(counts, marker_names, colors_markers)):
    ax.bar(x_pos + i*width, count, width, label=name, color=color, alpha=0.8, edgecolor='black', linewidth=1)

ax.set_xlabel('Patient', fontsize=12, fontweight='bold')
ax.set_ylabel('Cell Count', fontsize=12, fontweight='bold')
ax.set_title('Marker Co-expression Pattern', fontsize=13, fontweight='bold')
ax.set_xticks(x_pos + width * 1.5)
ax.set_xticklabels(df_results['Patient_ID'], fontsize=11, fontweight='bold')
ax.legend(fontsize=10, loc='upper right')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
detail_path = PLOTS_DIR / 'Multi_Patient_Poised_Validation_Detailed.png'
plt.savefig(detail_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"✓ Detailed comparison saved: {detail_path}")

plt.close()

print("\n✅ Enhanced visualizations created successfully!")
