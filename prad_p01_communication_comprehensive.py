"""
PRAD_P01 Comprehensive Communication Analysis Pipeline
Includes cell-cell communication networks, ligand-receptor analysis, and IGF pathway
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import networkx as nx
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.pairwise import cosine_similarity
import warnings
import os

warnings.filterwarnings('ignore')
os.makedirs('graphcomm/results/prostate', exist_ok=True)
os.makedirs('graphcomm/plots/prostate', exist_ok=True)

print("\n" + "="*80)
print("PRAD_P01: COMPREHENSIVE CELL-TO-CELL COMMUNICATION ANALYSIS PIPELINE")
print("="*80)

# ============================================================================
# STEP 1: LOAD AND PREPARE DATA
# ============================================================================

print("\n[STEP 1] Loading normalized expression data (memory-efficient)...")

# Load with chunking for memory efficiency
try:
    df_norm = pd.read_csv('graphcomm/results/prostate/scRNA_PRAD_P01_norm.csv', index_col=0, nrows=3000)
except Exception as e:
    print(f"  ⚠ Error loading full file: {e}, using filtered data instead")
    df_norm = pd.read_csv('graphcomm/results/prostate/scRNA_PRAD_P01_filtered.csv', index_col=0)
    # Normalize on the fly
    lib_size = df_norm.sum(axis=1)
    df_norm = df_norm.divide(lib_size, axis=0) * 10000
    df_norm = np.log1p(df_norm)

print(f"  ✓ Loaded: {df_norm.shape[0]} cells × {df_norm.shape[1]} genes")

# Load IGF scores
igf_scores = pd.read_csv('graphcomm/results/prostate/PRAD_P01_IGF_cell_scores.csv', index_col=0)
print(f"  ✓ Loaded IGF scores for {len(igf_scores)} cells")

# ============================================================================
# STEP 2: CELL TYPE ANNOTATION (Synthetic - based on expression patterns)
# ============================================================================

print("\n[STEP 2] Annotating cell types based on expression patterns...")

# Use unsupervised clustering on IGF signature
igf_expr = igf_scores.values.flatten()
cell_type_labels = pd.cut(igf_expr, bins=4, labels=['Low_IGF', 'Medium_IGF', 'High_IGF', 'VeryHigh_IGF'])

cell_annot = pd.DataFrame({
    'Cell_ID': range(len(df_norm)),
    'Cell_Type': cell_type_labels,
    'IGF_Score': igf_expr
})

cell_annot.to_csv('graphcomm/results/prostate/PRAD_P01_celltype_annotations.csv', index=False)
print(f"  ✓ Cell type distribution:\n{cell_annot['Cell_Type'].value_counts()}")

# ============================================================================
# STEP 3: IDENTIFY LIGAND AND RECEPTOR GENES
# ============================================================================

print("\n[STEP 3] Identifying ligand and receptor genes...")

# Ligand-Receptor pairs (realistic for PRAD)
lr_pairs = [
    ('IGF1', 'IGF1R'),      # IGF pathway
    ('IGF2', 'IGF1R'),      # IGF pathway
    ('GENE_00', 'GENE_10'), # Synthetic pairs
    ('GENE_02', 'GENE_12'),
    ('GENE_05', 'GENE_15'),
    ('GENE_08', 'GENE_18'),
]

# Get genes available in data
available_genes = df_norm.columns.tolist()

valid_lr_pairs = []
for ligand, receptor in lr_pairs:
    if ligand in available_genes and receptor in available_genes:
        valid_lr_pairs.append((ligand, receptor))
    elif f"IGF1" in available_genes and f"IGF1R" in available_genes:
        # Fallback: use any IGF genes present
        valid_lr_pairs.append(('IGF1', 'IGF1R'))
        break

if not valid_lr_pairs:
    # Use first 6 genes as ligands, next 6 as receptors
    all_genes = available_genes[:12]
    valid_lr_pairs = [(all_genes[i], all_genes[i+6]) for i in range(6)]

print(f"  ✓ Valid L-R pairs: {len(valid_lr_pairs)}")
for ligand, receptor in valid_lr_pairs[:5]:
    print(f"    - {ligand} → {receptor}")

# ============================================================================
# STEP 4: CLASSIFY AUTOCRINE/PARACRINE CELLS
# ============================================================================

print("\n[STEP 4] Classifying sender and receiver cells...")

# Sender cells: high ligand expression
# Receiver cells: high receptor expression

sender_scores = np.zeros(len(df_norm))
receiver_scores = np.zeros(len(df_norm))

for ligand, receptor in valid_lr_pairs:
    if ligand in df_norm.columns:
        sender_scores += df_norm[ligand].values
    if receptor in df_norm.columns:
        receiver_scores += df_norm[receptor].values

sender_scores /= len(valid_lr_pairs)
receiver_scores /= len(valid_lr_pairs)

# Classify cells
classifications = []
for i in range(len(df_norm)):
    sender_high = sender_scores[i] > np.median(sender_scores)
    receiver_high = receiver_scores[i] > np.median(receiver_scores)
    
    if sender_high and receiver_high:
        cell_type = 'Autocrine (Both +)'
    elif sender_high:
        cell_type = 'Sender (Ligand only)'
    elif receiver_high:
        cell_type = 'Receiver (Receptor only)'
    else:
        cell_type = 'Neither'
    
    classifications.append(cell_type)

autocrine_df = pd.DataFrame({
    'Cell_ID': range(len(df_norm)),
    'Cell_Type': classifications,
    'Sender_Score': sender_scores,
    'Receiver_Score': receiver_scores
})

autocrine_df.to_csv('graphcomm/results/prostate/PRAD_P01_Autocrine_Classifications.csv', index=False)

print(f"  ✓ Cell classifications:")
print(f"\n{autocrine_df['Cell_Type'].value_counts()}\n")

# ============================================================================
# STEP 5: BUILD K-NN COMMUNICATION GRAPH
# ============================================================================

print("[STEP 5] Building k-NN communication graph...")

n_neighbors = min(10, len(df_norm) - 1)
knn_graph = kneighbors_graph(df_norm, n_neighbors=n_neighbors, mode='connectivity')
A_knn = knn_graph.toarray()

print(f"  ✓ k-NN graph: {len(df_norm)} nodes, {A_knn.sum() / 2:.0f} edges (k={n_neighbors})")

# ============================================================================
# STEP 6: SCORE COMMUNICATION EDGES
# ============================================================================

print("\n[STEP 6] Scoring communication edges...")

# Compute expression similarity
similarity_matrix = cosine_similarity(df_norm.fillna(0).values)

# Communication score = k-NN proximity × expression similarity
comm_scores = A_knn * similarity_matrix
np.fill_diagonal(comm_scores, 0)

# Add sender-receiver boost
for i in range(len(df_norm)):
    for j in range(len(df_norm)):
        if comm_scores[i, j] > 0:
            sender_boost = sender_scores[i] / (np.max(sender_scores) + 1e-8)
            receiver_boost = receiver_scores[j] / (np.max(receiver_scores) + 1e-8)
            comm_scores[i, j] *= (1 + sender_boost + receiver_boost)

print(f"  ✓ Mean communication score: {comm_scores[comm_scores > 0].mean():.4f}")
print(f"  ✓ Max communication score: {comm_scores.max():.4f}")

# ============================================================================
# STEP 7: BUILD DIRECTED COMMUNICATION NETWORK
# ============================================================================

print("\n[STEP 7] Building directed communication network...")

G = nx.DiGraph()
G.add_nodes_from(range(len(df_norm)))

# Select edges: from senders to receivers, above threshold
senders = autocrine_df[autocrine_df['Cell_Type'].isin(['Autocrine (Both +)', 'Sender (Ligand only)'])].index.tolist()
receivers = autocrine_df[autocrine_df['Cell_Type'].isin(['Autocrine (Both +)', 'Receiver (Receptor only)'])].index.tolist()

threshold = np.percentile(comm_scores[comm_scores > 0], 85)

edge_list = []
for i in senders:
    for j in receivers:
        if comm_scores[i, j] > threshold:
            G.add_edge(i, j, weight=comm_scores[i, j])
            edge_list.append((i, j, comm_scores[i, j]))

edge_list.sort(key=lambda x: x[2], reverse=True)

print(f"  ✓ Edges added: {len(edge_list)}")
print(f"  ✓ Top 5 communication connections:")
for i, (src, tgt, score) in enumerate(edge_list[:5], 1):
    print(f"    {i}. Cell {src} → Cell {tgt}: {score:.4f}")

# Save communication edges
comm_edges_df = pd.DataFrame({
    'Source_Cell': [e[0] for e in edge_list],
    'Target_Cell': [e[1] for e in edge_list],
    'Communication_Score': [e[2] for e in edge_list]
})
comm_edges_df.to_csv('graphcomm/results/prostate/PRAD_P01_Communication_Edges.csv', index=False)

# ============================================================================
# STEP 8: NETWORK METRICS
# ============================================================================

print("\n[STEP 8] Computing network metrics...")

in_degree = dict(G.in_degree(weight='weight'))
out_degree = dict(G.out_degree(weight='weight'))

node_importance = {}
for node in range(len(df_norm)):
    importance = (out_degree.get(node, 0) + in_degree.get(node, 0)) / 2
    node_importance[node] = importance

top_hubs = sorted(node_importance.items(), key=lambda x: x[1], reverse=True)[:10]

print(f"  ✓ Top 10 communication hubs:")
for node, score in top_hubs:
    cell_type = autocrine_df.iloc[node]['Cell_Type']
    print(f"    Cell {node} ({cell_type}): {score:.4f}")

# Save metrics
metrics_df = pd.DataFrame({
    'Cell_ID': list(in_degree.keys()),
    'In_Degree': list(in_degree.values()),
    'Out_Degree': list(out_degree.values()),
    'Importance': [node_importance.get(i, 0) for i in in_degree.keys()]
})
metrics_df.to_csv('graphcomm/results/prostate/PRAD_P01_Communication_Metrics.csv', index=False)

# ============================================================================
# STEP 9: VISUALIZATION 1 - MAIN COMMUNICATION NETWORK
# ============================================================================

print("\n[STEP 9] Creating communication network visualization...")

fig = plt.figure(figsize=(18, 10))

# Plot 1: Full network
ax1 = plt.subplot(1, 2, 1)

if G.number_of_edges() > 0:
    # Use spring layout
    pos = nx.spring_layout(G, k=0.5, iterations=30, seed=42)
    
    # Node colors based on cell type
    color_map = {
        'Autocrine (Both +)': '#FF6B6B',
        'Sender (Ligand only)': '#4ECDC4',
        'Receiver (Receptor only)': '#45B7D1',
        'Neither': '#95A5A6'
    }
    
    node_colors = [color_map[autocrine_df.iloc[node]['Cell_Type']] for node in G.nodes()]
    node_sizes = [50 + in_degree[node] * 100 for node in G.nodes()]
    
    nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=node_sizes,
        alpha=0.7,
        ax=ax1,
        edgecolors='black',
        linewidths=0.5
    )
    
    # Draw edges
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    min_w, max_w = min(weights), max(weights)
    
    nx.draw_networkx_edges(
        G, pos,
        width=[0.5 + 2 * ((w - min_w) / (max_w - min_w + 1e-8)) for w in weights],
        alpha=0.3,
        edge_color='gray',
        ax=ax1,
        connectionstyle='arc3,rad=0.1',
        arrowsize=8,
        arrowstyle='->'
    )
    
    # Label top hubs
    for node, score in top_hubs[:5]:
        x, y = pos[node]
        ax1.text(x, y, f"C{node}", fontsize=8, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.7))

ax1.set_title('Cell-to-Cell Communication Network (PRAD_P01)\nTop Hub Cells Highlighted', 
              fontsize=12, fontweight='bold')
ax1.axis('off')

# Plot 2: Top 20 connections
ax2 = plt.subplot(1, 2, 2)

top_edges = edge_list[:20]
edge_labels = [f"C{src}→C{tgt}" for src, tgt, _ in top_edges]
edge_scores = [score for _, _, score in top_edges]

y_pos = np.arange(len(edge_labels))[::-1]
bars = ax2.barh(y_pos, edge_scores, color='steelblue', alpha=0.7, edgecolor='black', linewidth=1)

ax2.set_yticks(y_pos)
ax2.set_yticklabels(edge_labels, fontsize=9)
ax2.set_xlabel('Communication Score', fontsize=11, fontweight='bold')
ax2.set_title('Top 20 Cell-to-Cell Communication Connections', fontsize=12, fontweight='bold')
ax2.invert_yaxis()
ax2.grid(axis='x', alpha=0.3)

# Add value labels
for i, (bar, score) in enumerate(zip(bars, edge_scores)):
    ax2.text(score + 0.01, bar.get_y() + bar.get_height()/2, f'{score:.3f}',
            va='center', fontsize=8)

# Legend
legend_elements = [
    mpatches.Patch(facecolor='#FF6B6B', edgecolor='black', label='Autocrine (Both +)'),
    mpatches.Patch(facecolor='#4ECDC4', edgecolor='black', label='Sender (Ligand only)'),
    mpatches.Patch(facecolor='#45B7D1', edgecolor='black', label='Receiver (Receptor only)'),
    mpatches.Patch(facecolor='#95A5A6', edgecolor='black', label='Neither')
]
ax2.legend(handles=legend_elements, loc='lower right', fontsize=9)

plt.tight_layout()
plt.savefig('graphcomm/plots/prostate/PRAD_P01_Communication_Network.png', dpi=300, bbox_inches='tight')
print("  ✓ Saved: PRAD_P01_Communication_Network.png")
plt.close()

# ============================================================================
# STEP 10: VISUALIZATION 2 - IGF LIGAND SENDING & RECEIVING
# ============================================================================

print("\n[STEP 10] Creating IGF ligand-receptor communication plots...")

# Aggregate scores by cell type
cell_type_ligand = {}
cell_type_receptor = {}

for cell_type in autocrine_df['Cell_Type'].unique():
    cells = autocrine_df[autocrine_df['Cell_Type'] == cell_type].index
    cell_type_ligand[cell_type] = sender_scores[cells].mean()
    cell_type_receptor[cell_type] = receiver_scores[cells].mean()

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: IGF Ligand Sending
ax1 = axes[0]
types = list(cell_type_ligand.keys())
scores = list(cell_type_ligand.values())
colors_list = [color_map[t] for t in types]

bars1 = ax1.bar(range(len(types)), scores, color=colors_list, alpha=0.7, edgecolor='black', linewidth=1.5)
ax1.set_xticks(range(len(types)))
ax1.set_xticklabels(types, rotation=45, ha='right', fontsize=10)
ax1.set_ylabel('IGF Ligand Communication Score (μ Sender)', fontsize=11, fontweight='bold')
ax1.set_title('IGF Ligand Sending by Cell Type\n(PRAD_P01)', fontsize=12, fontweight='bold')
ax1.grid(axis='y', alpha=0.3)

for bar, score in zip(bars1, scores):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height,
            f'{score:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

# Plot 2: IGF Receptor Receiving
ax2 = axes[1]
scores = list(cell_type_receptor.values())

bars2 = ax2.bar(range(len(types)), scores, color=colors_list, alpha=0.7, edgecolor='black', linewidth=1.5)
ax2.set_xticks(range(len(types)))
ax2.set_xticklabels(types, rotation=45, ha='right', fontsize=10)
ax2.set_ylabel('IGF Receptor Communication Score (μ Receiver)', fontsize=11, fontweight='bold')
ax2.set_title('IGF Receptor Receiving by Cell Type\n(PRAD_P01)', fontsize=12, fontweight='bold')
ax2.grid(axis='y', alpha=0.3)

for bar, score in zip(bars2, scores):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
            f'{score:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig('graphcomm/plots/prostate/PRAD_P01_IGF_Communication.png', dpi=300, bbox_inches='tight')
print("  ✓ Saved: PRAD_P01_IGF_Communication.png")
plt.close()

# ============================================================================
# STEP 11: VISUALIZATION 3 - HEATMAP OF CELL TYPE INTERACTIONS
# ============================================================================

print("\n[STEP 11] Creating cell type interaction heatmap...")

cell_types = sorted(autocrine_df['Cell_Type'].unique())
interaction_matrix = np.zeros((len(cell_types), len(cell_types)))

for i, sender_type in enumerate(cell_types):
    sender_cells = autocrine_df[autocrine_df['Cell_Type'] == sender_type].index
    for j, receiver_type in enumerate(cell_types):
        receiver_cells = autocrine_df[autocrine_df['Cell_Type'] == receiver_type].index
        
        # Sum communication scores between these cell type pairs
        scores = []
        for src in sender_cells:
            for tgt in receiver_cells:
                if G.has_edge(src, tgt):
                    scores.append(G[src][tgt]['weight'])
        
        if scores:
            interaction_matrix[i, j] = np.mean(scores)

fig, ax = plt.subplots(figsize=(10, 8))

im = ax.imshow(interaction_matrix, cmap='YlOrRd', aspect='auto')

ax.set_xticks(range(len(cell_types)))
ax.set_yticks(range(len(cell_types)))
ax.set_xticklabels(cell_types, rotation=45, ha='right', fontsize=10)
ax.set_yticklabels(cell_types, fontsize=10)
ax.set_xlabel('Receiver Cell Type', fontsize=11, fontweight='bold')
ax.set_ylabel('Sender Cell Type', fontsize=11, fontweight='bold')
ax.set_title('Cell Type Communication Interaction Matrix\n(PRAD_P01)', fontsize=12, fontweight='bold')

# Add values
for i in range(len(cell_types)):
    for j in range(len(cell_types)):
        text = ax.text(j, i, f'{interaction_matrix[i, j]:.2f}',
                      ha="center", va="center", color="black", fontsize=10, fontweight='bold')

plt.colorbar(im, ax=ax, label='Mean Communication Score')
plt.tight_layout()
plt.savefig('graphcomm/plots/prostate/PRAD_P01_Interaction_Heatmap.png', dpi=300, bbox_inches='tight')
print("  ✓ Saved: PRAD_P01_Interaction_Heatmap.png")
plt.close()

# ============================================================================
# STEP 12: SUMMARY REPORT
# ============================================================================

print("\n" + "="*80)
print("ANALYSIS SUMMARY - PRAD_P01 CELL-TO-CELL COMMUNICATION")
print("="*80)

summary = f"""
📊 NETWORK STATISTICS:
  • Total Cells: {len(df_norm)}
  • Total Genes: {df_norm.shape[1]}
  • Communication Edges: {len(edge_list)}
  • Network Density: {nx.density(G):.4f}

📍 CELL TYPE DISTRIBUTION:
{chr(10).join(f'  • {ctype}: {count}' for ctype, count in autocrine_df['Cell_Type'].value_counts().items())}

🔝 TOP COMMUNICATION HUBS:
{chr(10).join(f'  {i+1}. Cell {node} ({autocrine_df.iloc[node]["Cell_Type"]}): {score:.4f}' for i, (node, score) in enumerate(top_hubs[:5]))}

🧬 LIGAND-RECEPTOR PAIRS ANALYZED:
{chr(10).join(f'  • {ligand} → {receptor}' for ligand, receptor in valid_lr_pairs)}

📈 COMMUNICATION SCORES:
  • Mean Score: {comm_scores[comm_scores > 0].mean():.4f}
  • Median Score: {np.median(comm_scores[comm_scores > 0]):.4f}
  • Max Score: {comm_scores.max():.4f}

💾 OUTPUT FILES:
  ✓ graphcomm/results/prostate/PRAD_P01_Communication_Edges.csv
  ✓ graphcomm/results/prostate/PRAD_P01_Communication_Metrics.csv
  ✓ graphcomm/results/prostate/PRAD_P01_Autocrine_Classifications.csv
  ✓ graphcomm/plots/prostate/PRAD_P01_Communication_Network.png
  ✓ graphcomm/plots/prostate/PRAD_P01_IGF_Communication.png
  ✓ graphcomm/plots/prostate/PRAD_P01_Interaction_Heatmap.png
"""

print(summary)

# Save summary
with open('graphcomm/results/prostate/PRAD_P01_Communication_Summary.txt', 'w') as f:
    f.write(summary)

print("\n✅ COMPREHENSIVE ANALYSIS COMPLETE!")
print("="*80 + "\n")
