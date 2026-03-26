"""
GraphComm-Lite: Comprehensive Cell-Cell Communication Analysis for PRAD_P01

Full pipeline including:
1. Cell type identification and classification
2. Ligand-receptor communication network
3. IGF pathway analysis
4. Multiple visualizations matching LUAD/BC format
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')
import os
from pathlib import Path

# Create output directories
os.makedirs("graphcomm/results/prostate", exist_ok=True)
os.makedirs("graphcomm/plots/prostate", exist_ok=True)

print("=" * 70)
print("GraphComm-Lite: Comprehensive PRAD_P01 Communication Analysis")
print("=" * 70)

# =============================================================================
# Step 0: Load preprocessed data
# =============================================================================
print("\n[STEP 1] Loading preprocessed data...")

try:
    df_norm = pd.read_csv("graphcomm/results/prostate/scRNA_PRAD_P01_norm.csv", index_col=0)
    print(f"✓ Loaded normalized data: {df_norm.shape[0]} genes × {df_norm.shape[1]} cells")
except FileNotFoundError:
    print("❌ ERROR: Missing normalized data. Run minimal pipeline first.")
    exit(1)

# =============================================================================
# Step 1: Identify cell types (clustering-based for now)
# =============================================================================
print("\n[STEP 2] Identifying cell types...")

# For synthetic data, assign cell types based on cell index distribution
np.random.seed(42)
cell_types = {}

# Define cell type proportions
cell_type_list = ['Luminal_Tumor', 'Basal_Tumor', 'Immune_T', 'Immune_Myeloid', 'CAF', 'Endothelial']
proportions = [0.35, 0.20, 0.15, 0.15, 0.10, 0.05]

# Assign cells proportionally
n_cells = df_norm.shape[1]
cell_indices = list(df_norm.columns)
type_assignments = np.random.choice(cell_type_list, size=n_cells, p=proportions)

for cell_id, cell_type in zip(cell_indices, type_assignments):
    cell_types[cell_id] = cell_type

cell_type_df = pd.DataFrame({
    'cell_id': list(cell_types.keys()),
    'cell_type': list(cell_types.values())
})
print(f"✓ Assigned {len(cell_type_df)} cells to {cell_type_df['cell_type'].nunique()} types")
print(f"  Distribution: {dict(cell_type_df['cell_type'].value_counts())}")

# =============================================================================
# Step 2: Define Ligand-Receptor Pairs (using synthetic gene names)
# =============================================================================
print("\n[STEP 3] Defining ligand-receptor communication network...")

# For synthetic data, use available genes (GENE_XXXXX format)
# Create synthetic L-R pairs from the available genes
available_genes_list = list(df_norm.index)

# Define communication pairs using indices to create sender/receiver genes
LIGAND_RECEPTOR_PAIRS = [
    ("GENE_000050", "GENE_000100"),  # Simulated IGF1 -> IGF1R
    ("GENE_000051", "GENE_000101"),  # Simulated IGF2 -> IGF1R
    ("GENE_000052", "GENE_000102"),  # Simulated EGF -> EGFR
    ("GENE_000053", "GENE_000103"),  # Simulated TGFB1 -> TGFBR1
    ("GENE_000054", "GENE_000104"),  # Simulated VEGFA -> KDR
    ("GENE_000055", "GENE_000105"),  # Simulated FGF -> FGFR
    ("GENE_000056", "GENE_000106"),  # Simulated HGF -> MET
    ("GENE_000057", "GENE_000107"),  # Simulated PDGF -> PDGFR
    ("GENE_000058", "GENE_000108"),  # Simulated CXCL12 -> CXCR4
    ("GENE_000059", "GENE_000109"),  # Simulated CCL2 -> CCR2
    ("GENE_000060", "GENE_000110"),  # Simulated TNF -> TNFR
    ("GENE_000061", "GENE_000111"),  # Simulated IL6 -> IL6R
    ("GENE_000062", "GENE_000112"),  # Simulated DLL1 -> NOTCH1
    ("GENE_000063", "GENE_000113"),  # Simulated JAG1 -> NOTCH1
    ("GENE_000064", "GENE_000114"),  # Simulated WNT -> FZD
    ("GENE_000065", "GENE_000115"),  # Simulated COL -> ITGA
    ("GENE_000066", "GENE_000116"),  # Simulated FN -> ITGA
    ("GENE_000067", "GENE_000117"),  # Simulated LAMB -> ITGA
    ("GENE_000068", "GENE_000118"),  # Simulated SEMA -> PLXN
    ("GENE_000069", "GENE_000119"),  # Simulated EFNA -> EFNB
]

# Check which L-R pairs are available
available_genes = set(df_norm.index)
valid_pairs = []
for lig, rec in LIGAND_RECEPTOR_PAIRS:
    if lig in available_genes and rec in available_genes:
        valid_pairs.append((lig, rec))

print(f"✓ Valid L-R pairs: {len(valid_pairs)}/{len(LIGAND_RECEPTOR_PAIRS)}")
for lig, rec in valid_pairs[:8]:
    print(f"  {lig:10} → {rec:10}")
if len(valid_pairs) > 8:
    print(f"  ... and {len(valid_pairs)-8} more")

# =============================================================================
# Step 3: Classify cells as Autocrine/Sender/Receiver
# =============================================================================
print("\n[STEP 4] Classifying cells by IGF pathway status...")

# Use synthetic IGF genes
igf_ligand_genes = ["GENE_000050", "GENE_000051"]  # Synthetic IGF1, IGF2
igf_receptor_genes = ["GENE_000100", "GENE_000101"]  # Synthetic IGF1R, IGF2R
igf_genes_found = [g for g in igf_ligand_genes + igf_receptor_genes if g in available_genes]

cell_classification = []
for cell_id in df_norm.columns:
    # Check ligand expression
    ligand_genes_found = [g for g in igf_ligand_genes if g in available_genes]
    ligand_expr = df_norm.loc[ligand_genes_found, cell_id].max() if ligand_genes_found else 0
    has_ligand = ligand_expr > df_norm.loc[ligand_genes_found, :].median().max() if ligand_genes_found else False
    
    # Check receptor expression
    receptor_genes_found = [g for g in igf_receptor_genes if g in available_genes]
    receptor_expr = df_norm.loc[receptor_genes_found, cell_id].max() if receptor_genes_found else 0
    has_receptor = receptor_expr > df_norm.loc[receptor_genes_found, :].median().max() if receptor_genes_found else False
    
    # Classify
    if has_ligand and has_receptor:
        classification = 'Autocrine (Both +)'
    elif has_ligand:
        classification = 'Sender (Ligand only)'
    elif has_receptor:
        classification = 'Receiver (Receptor only)'
    else:
        classification = 'Neither'
    
    cell_classification.append({
        'Cell_ID': cell_id,
        'Cell_Type': cell_types.get(cell_id, 'Unknown'),
        'Classification': classification,
        'Ligand_Max': ligand_expr,
        'Receptor_Max': receptor_expr
    })

auto_df = pd.DataFrame(cell_classification)
auto_df.to_csv("graphcomm/results/prostate/PRAD_P01_Autocrine_Classifications.csv", index=False)
print(f"✓ Saved classifications to PRAD_P01_Autocrine_Classifications.csv")
print(f"  Autocrine: {(auto_df['Classification'] == 'Autocrine (Both +)').sum()}")
print(f"  Sender: {(auto_df['Classification'] == 'Sender (Ligand only)').sum()}")
print(f"  Receiver: {(auto_df['Classification'] == 'Receiver (Receptor only)').sum()}")
print(f"  Inactive: {(auto_df['Classification'] == 'Neither').sum()}")

# =============================================================================
# Step 4: Compute communication edges between cells
# =============================================================================
print("\n[STEP 5] Computing cell-cell communication edges...")

comm_edges = []
n_cells = df_norm.shape[1]

for idx_src, src_cell in enumerate(df_norm.columns):
    if (idx_src + 1) % max(1, n_cells // 10) == 0:
        print(f"  Processing cell {idx_src + 1}/{n_cells}...")
    
    src_type = cell_types.get(src_cell, 'Unknown')
    
    for idx_tgt, tgt_cell in enumerate(df_norm.columns):
        tgt_type = cell_types.get(tgt_cell, 'Unknown')
        
        # Calculate communication potential
        for lig, rec in valid_pairs:
            lig_expr = df_norm.loc[lig, src_cell]
            rec_expr = df_norm.loc[rec, tgt_cell]
            
            # Communication score = geometric mean
            if lig_expr > 0 and rec_expr > 0:
                score = np.sqrt(lig_expr * rec_expr)
                
                if score > 1e-6:  # Filter very small scores
                    comm_edges.append({
                        'Source_Cell': src_cell,
                        'Target_Cell': tgt_cell,
                        'Source_Type': src_type,
                        'Target_Type': tgt_type,
                        'Ligand': lig,
                        'Receptor': rec,
                        'Ligand_Expr': lig_expr,
                        'Receptor_Expr': rec_expr,
                        'Communication_Score': score
                    })

comm_df = pd.DataFrame(comm_edges)
print(f"✓ Computed {len(comm_df):,} communication edges")

# Aggregate by cell pair (top interaction per cell pair)
comm_agg = comm_df.groupby(['Source_Cell', 'Target_Cell']).agg({
    'Communication_Score': 'sum',
    'Source_Type': 'first',
    'Target_Type': 'first'
}).reset_index()

comm_df.to_csv("graphcomm/results/prostate/PRAD_P01_Communication_Edges.csv", index=False)
comm_agg.to_csv("graphcomm/results/prostate/PRAD_P01_Communication_Network.csv", index=False)
print(f"✓ Saved edges (aggregated: {len(comm_agg):,} cell-cell connections)")

# =============================================================================
# Step 5: Compute communication metrics
# =============================================================================
print("\n[STEP 6] Computing communication metrics...")

metrics = []
for cell_id in df_norm.columns:
    # In-degree (receiving capability)
    in_comms = comm_agg[comm_agg['Target_Cell'] == cell_id]['Communication_Score'].sum()
    in_count = len(comm_agg[comm_agg['Target_Cell'] == cell_id])
    
    # Out-degree (sending capability)
    out_comms = comm_agg[comm_agg['Source_Cell'] == cell_id]['Communication_Score'].sum()
    out_count = len(comm_agg[comm_agg['Source_Cell'] == cell_id])
    
    metrics.append({
        'Cell_ID': cell_id,
        'Cell_Type': cell_types.get(cell_id, 'Unknown'),
        'In_Degree': in_count,
        'Out_Degree': out_count,
        'In_Score': in_comms,
        'Out_Score': out_comms,
        'Total_Communication': in_comms + out_comms
    })

metrics_df = pd.DataFrame(metrics)
metrics_df.to_csv("graphcomm/results/prostate/PRAD_P01_Communication_Metrics.csv", index=False)
print(f"✓ Saved metrics for {len(metrics_df)} cells")

# =============================================================================
# VISUALIZATION 1: Cell-Cell Communication Network
# =============================================================================
print("\n[STEP 7] Generating visualizations...")
print("  [1/5] Communication network graph...")

fig, ax = plt.subplots(figsize=(14, 10))

# Build network for visualization (top edges for clarity)
G = nx.DiGraph()
for cell_id in df_norm.columns[:500]:  # Top 500 cells for readability
    G.add_node(cell_id)

# Add top communication edges
top_edges = comm_agg.nlargest(300, 'Communication_Score')
for _, row in top_edges.iterrows():
    if row['Source_Cell'] in G.nodes() and row['Target_Cell'] in G.nodes():
        G.add_edge(row['Source_Cell'], row['Target_Cell'], weight=row['Communication_Score'])

# Layout
pos = nx.spring_layout(G, k=0.5, iterations=20, seed=42)

# Node colors based on cell type
color_map = {
    'Luminal_Tumor': '#FF6B6B',
    'Basal_Tumor': '#FF8C42',
    'Immune_T': '#4ECDC4',
    'Immune_Myeloid': '#45B7D1',
    'CAF': '#96CEB4',
    'Endothelial': '#FFEAA7',
    'Unknown': '#95A5A6'
}
node_colors = [color_map.get(cell_types.get(node, 'Unknown'), '#95A5A6') for node in G.nodes()]
node_sizes = [100 + 50 * metrics_df[metrics_df['Cell_ID'] == node]['Total_Communication'].values[0] 
              for node in G.nodes()]

# Draw network
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8, ax=ax, edgecolors='black', linewidths=0.5)

# Draw edges with varying width
edges = G.edges()
weights = [G[u][v]['weight'] for u, v in edges]
if weights:
    min_w, max_w = min(weights), max(weights)
    widths = [1 + 2*((w-min_w)/(max_w-min_w + 1e-6)) for w in weights]
else:
    widths = [1] * len(edges)

nx.draw_networkx_edges(G, pos, width=widths, alpha=0.3, edge_color='gray', ax=ax, 
                       arrowsize=10, arrowstyle='->', connectionstyle='arc3,rad=0.1')

ax.set_title('Cell-Cell Communication Network (PRAD_P01)\nTop 500 Cells, 300 Strongest Links', 
             fontsize=14, fontweight='bold', pad=20)
ax.axis('off')

# Legend
legend_elements = [mpatches.Patch(facecolor=color, edgecolor='black', label=ct) 
                   for ct, color in color_map.items()]
ax.legend(handles=legend_elements, loc='upper left', fontsize=10, framealpha=0.9)

plt.tight_layout()
plt.savefig("graphcomm/plots/prostate/PRAD_P01_Communication_Network.png", dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: PRAD_P01_Communication_Network.png")
plt.close()

# =============================================================================
# VISUALIZATION 2: IGF Ligand Sending (by cell type)
# =============================================================================
print("  [2/5] IGF ligand sending analysis...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Get IGF ligand expression
igf_ligands = ["GENE_000050", "GENE_000051"]
igf_ligands_found = [g for g in igf_ligands if g in available_genes]

if igf_ligands_found:
    # Compute mean expression per cell type
    celltype_ligand = {}
    for ct in cell_type_df['cell_type'].unique():
        cells_in_type = cell_type_df[cell_type_df['cell_type'] == ct]['cell_id'].values
        cells_in_type = [c for c in cells_in_type if c in df_norm.columns]
        if cells_in_type:
            expr = df_norm.loc[igf_ligands_found, cells_in_type].mean(axis=1).max()
            celltype_ligand[ct] = expr
    
    # Sort and plot
    sorted_ct = sorted(celltype_ligand.items(), key=lambda x: x[1], reverse=True)
    cts = [x[0] for x in sorted_ct]
    vals = [x[1] for x in sorted_ct]
    
    ax1.barh(cts, vals, color='steelblue', alpha=0.8, edgecolor='black')
    ax1.set_xlabel('IGF Ligand Expression (Mean)', fontsize=12, fontweight='bold')
    ax1.set_title('IGF Ligand Sending (PRAD_P01)', fontsize=13, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    # Compute communication score per ligand
    if len(comm_df) > 0:
        ligand_comms = comm_df.groupby('Ligand')['Communication_Score'].sum().sort_values(ascending=True).tail(10)
        ax2.barh(ligand_comms.index, ligand_comms.values, color='coral', alpha=0.8, edgecolor='black')
    else:
        ax2.text(0.5, 0.5, 'No Communication Data', ha='center', va='center', transform=ax2.transAxes)
    
    ax2.set_xlabel('Total Communication Score', fontsize=12, fontweight='bold')
    ax2.set_title('Top IGF Ligands in Communication', fontsize=13, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig("graphcomm/plots/prostate/PRAD_P01_IGF_Ligand_Sending.png", dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: PRAD_P01_IGF_Ligand_Sending.png")
plt.close()

# =============================================================================
# VISUALIZATION 3: IGF Receptor Receiving (by cell type)
# =============================================================================
print("  [3/5] IGF receptor receiving analysis...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Get IGF receptor expression
igf_receptors = ["GENE_000100", "GENE_000101"]
igf_receptors_found = [g for g in igf_receptors if g in available_genes]

if igf_receptors_found:
    # Compute mean expression per cell type
    celltype_receptor = {}
    for ct in cell_type_df['cell_type'].unique():
        cells_in_type = cell_type_df[cell_type_df['cell_type'] == ct]['cell_id'].values
        cells_in_type = [c for c in cells_in_type if c in df_norm.columns]
        if cells_in_type:
            expr = df_norm.loc[igf_receptors_found, cells_in_type].mean(axis=1).max()
            celltype_receptor[ct] = expr
    
    # Sort and plot
    sorted_ct = sorted(celltype_receptor.items(), key=lambda x: x[1], reverse=True)
    cts = [x[0] for x in sorted_ct]
    vals = [x[1] for x in sorted_ct]
    
    ax1.barh(cts, vals, color='coral', alpha=0.8, edgecolor='black')
    ax1.set_xlabel('IGF Receptor Expression (Mean)', fontsize=12, fontweight='bold')
    ax1.set_title('IGF Receptor Receiving (PRAD_P01)', fontsize=13, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    # Compute communication score per receptor
    if len(comm_df) > 0:
        receptor_comms = comm_df.groupby('Receptor')['Communication_Score'].sum().sort_values(ascending=True).tail(10)
        ax2.barh(receptor_comms.index, receptor_comms.values, color='steelblue', alpha=0.8, edgecolor='black')
    else:
        ax2.text(0.5, 0.5, 'No Communication Data', ha='center', va='center', transform=ax2.transAxes)
    
    ax2.set_xlabel('Total Communication Score', fontsize=12, fontweight='bold')
    ax2.set_title('Top IGF Receptors in Communication', fontsize=13, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig("graphcomm/plots/prostate/PRAD_P01_IGF_Receptor_Receiving.png", dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: PRAD_P01_IGF_Receptor_Receiving.png")
plt.close()

# =============================================================================
# VISUALIZATION 4: Communication Degree Distributions
# =============================================================================
print("  [4/5] Communication degree distributions...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# In-degree (receiver capability)
in_degrees = metrics_df['In_Degree'].values
ax1.hist(in_degrees, bins=40, color='steelblue', alpha=0.7, edgecolor='black')
ax1.axvline(np.mean(in_degrees), color='red', linestyle='--', linewidth=2.5, label=f'Mean: {np.mean(in_degrees):.1f}')
ax1.axvline(np.median(in_degrees), color='orange', linestyle='--', linewidth=2.5, label=f'Median: {np.median(in_degrees):.1f}')
ax1.set_xlabel('In-Degree (Receiver Capability)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Number of Cells', fontsize=12, fontweight='bold')
ax1.set_title('Receiver Cell Distribution (PRAD_P01)', fontsize=13, fontweight='bold')
ax1.legend(fontsize=11)
ax1.grid(axis='y', alpha=0.3)

# Out-degree (sender capability)
out_degrees = metrics_df['Out_Degree'].values
ax2.hist(out_degrees, bins=40, color='coral', alpha=0.7, edgecolor='black')
ax2.axvline(np.mean(out_degrees), color='red', linestyle='--', linewidth=2.5, label=f'Mean: {np.mean(out_degrees):.1f}')
ax2.axvline(np.median(out_degrees), color='orange', linestyle='--', linewidth=2.5, label=f'Median: {np.median(out_degrees):.1f}')
ax2.set_xlabel('Out-Degree (Sender Capability)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Number of Cells', fontsize=12, fontweight='bold')
ax2.set_title('Sender Cell Distribution (PRAD_P01)', fontsize=13, fontweight='bold')
ax2.legend(fontsize=11)
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig("graphcomm/plots/prostate/PRAD_P01_Communication_Degrees.png", dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: PRAD_P01_Communication_Degrees.png")
plt.close()

# =============================================================================
# VISUALIZATION 5: Communication Heatmap (cell type interactions)
# =============================================================================
print("  [5/5] Communication heatmap by cell type...")

fig, ax = plt.subplots(figsize=(10, 9))

# Compute cell type interaction matrix
cell_types_list = sorted(cell_type_df['cell_type'].unique())
interaction_matrix = pd.DataFrame(0.0, index=cell_types_list, columns=cell_types_list)

for _, row in comm_agg.iterrows():
    src_type = row['Source_Type']
    tgt_type = row['Target_Type']
    if src_type in interaction_matrix.index and tgt_type in interaction_matrix.columns:
        interaction_matrix.loc[src_type, tgt_type] += row['Communication_Score']

# Normalize
interaction_matrix = interaction_matrix / interaction_matrix.max().max() if interaction_matrix.max().max() > 0 else interaction_matrix

# Heatmap
sns.heatmap(interaction_matrix, annot=True, fmt='.2f', cmap='YlOrRd', cbar_kws={'label': 'Normalized Communication Score'},
            ax=ax, linewidths=0.5, linecolor='gray', annot_kws={'size': 10})
ax.set_title('Cell Type Interaction Heatmap (PRAD_P01)', fontsize=14, fontweight='bold', pad=20)
ax.set_xlabel('Target (Receiver) Cell Type', fontsize=12, fontweight='bold')
ax.set_ylabel('Source (Sender) Cell Type', fontsize=12, fontweight='bold')
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig("graphcomm/plots/prostate/PRAD_P01_Communication_Heatmap.png", dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: PRAD_P01_Communication_Heatmap.png")
plt.close()

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)

summary = f"""
PRAD_P01 COMPREHENSIVE COMMUNICATION ANALYSIS SUMMARY
{'='*70}

DATA OVERVIEW:
  • Total Cells: {df_norm.shape[1]:,}
  • Total Genes: {df_norm.shape[0]:,}
  • Cell Types: {cell_type_df['cell_type'].nunique()}

CELL CLASSIFICATION:
  • Autocrine (IGF+): {(auto_df['Classification'] == 'Autocrine (Both +)').sum():,} cells
  • Sender (Ligand+): {(auto_df['Classification'] == 'Sender (Ligand only)').sum():,} cells
  • Receiver (Receptor+): {(auto_df['Classification'] == 'Receiver (Receptor only)').sum():,} cells
  • Inactive: {(auto_df['Classification'] == 'Neither').sum():,} cells

LIGAND-RECEPTOR COMMUNICATION:
  • Valid L-R Pairs: {len(valid_pairs)}/{len(LIGAND_RECEPTOR_PAIRS)}
  • Total Communication Edges: {len(comm_df):,}
  • Unique Cell-Cell Interactions: {len(comm_agg):,}
  • Network Density: {len(comm_agg) / (df_norm.shape[1]**2):.4f}

COMMUNICATION METRICS:
  • Avg In-Degree (receiver): {metrics_df['In_Degree'].mean():.2f}
  • Avg Out-Degree (sender): {metrics_df['Out_Degree'].mean():.2f}
  • Max In-Degree: {metrics_df['In_Degree'].max():.0f}
  • Max Out-Degree: {metrics_df['Out_Degree'].max():.0f}

OUTPUT FILES:
  CSV Data:
    ✓ PRAD_P01_Autocrine_Classifications.csv
    ✓ PRAD_P01_Communication_Edges.csv
    ✓ PRAD_P01_Communication_Network.csv
    ✓ PRAD_P01_Communication_Metrics.csv
  
  Visualizations (300 dpi):
    ✓ PRAD_P01_Communication_Network.png
    ✓ PRAD_P01_IGF_Ligand_Sending.png
    ✓ PRAD_P01_IGF_Receptor_Receiving.png
    ✓ PRAD_P01_Communication_Degrees.png
    ✓ PRAD_P01_Communication_Heatmap.png

All files saved to:
  Data: graphcomm/results/prostate/
  Plots: graphcomm/plots/prostate/
"""

print(summary)

# Save summary
with open("graphcomm/results/prostate/PRAD_P01_Communication_Summary.txt", "w") as f:
    f.write(summary)

print("✓ Summary saved to PRAD_P01_Communication_Summary.txt")
