"""
LUAD_P01 Cell-to-Cell Communication Analysis
Identifies sender-receiver relationships and communication networks
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
import os

print("=" * 70)
print("LUAD_P01: Cell-to-Cell Communication Analysis")
print("=" * 70)

# ============================================================
# Load preprocessed data and IGF scores
# ============================================================

print("\nStep 1: Loading data")
norm_path = "data/luad/processed/scRNA_LUAD_P01_norm.parquet"
igf_scores_path = "results/LUAD_P01_IGF_cell_scores.csv"
celltype_path = "results/LUAD_P01_celltype_annotations.csv"

if not os.path.exists(norm_path):
    print(f"ERROR: {norm_path} not found!")
    exit(1)

df_norm = pd.read_parquet(norm_path)
igf_scores_df = pd.read_csv(igf_scores_path)
celltype_df = pd.read_csv(celltype_path)

print(f"Data shape: {df_norm.shape}")
print(f"IGF scores: {len(igf_scores_df)}")
print(f"Cell types: {len(celltype_df)}")

# Map IGF scores and cell types to cells
igf_dict = dict(zip(igf_scores_df['Cell'], igf_scores_df['IGF_Score']))
celltype_dict = dict(zip(celltype_df['Cell'], celltype_df['Cell_Type']))

# Create cell metadata
cell_metadata = pd.DataFrame({
    'Cell': df_norm.index,
    'IGF_Score': [igf_dict.get(c, 0) for c in df_norm.index],
    'Cell_Type': [celltype_dict.get(c, 'Unknown') for c in df_norm.index]
})

print(f"\nCell metadata created: {cell_metadata.shape}")

# ============================================================
# Step 2: Build k-NN graph for communication
# ============================================================

print("\nStep 2: Building k-NN communication graph")

from sklearn.neighbors import kneighbors_graph

n_neighbors = min(10, df_norm.shape[0] - 1)
print(f"  k-NN: {n_neighbors} neighbors")

# Compute k-NN graph on normalized expression
knn_graph = kneighbors_graph(df_norm, n_neighbors=n_neighbors, mode='distance')

# Convert to adjacency matrix (undirected)
adjacency = (knn_graph + knn_graph.T) > 0
adjacency = adjacency.astype(float)
print(f"  Adjacency matrix: {adjacency.shape}")
print(f"  Edges: {int(adjacency.sum() / 2)}")

# ============================================================
# Step 3: Compute sender-receiver scores
# ============================================================

print("\nStep 3: Computing sender-receiver communication scores")

n_cells = df_norm.shape[0]

# Initialize communication matrix
communication_scores = np.zeros((n_cells, n_cells))

# Identify ligands and receptors (use top expressing genes)
top_genes_idx = df_norm.mean(axis=0).nlargest(50).index
ligand_genes = top_genes_idx[:25]  # Top 25 as ligands
receptor_genes = top_genes_idx[25:50]  # Next 25 as receptors

print(f"  Ligand genes (top senders): {ligand_genes[:5].tolist()} ...")
print(f"  Receptor genes (top receivers): {receptor_genes[:5].tolist()} ...")

# For each pair of connected cells, compute communication score
for i in range(n_cells):
    for j in range(n_cells):
        if i == j:
            continue
        if adjacency[i, j] > 0:  # Only connected cells
            # Sender (cell i) expresses ligands
            sender_ligand_expr = df_norm.loc[df_norm.index[i], ligand_genes].mean()
            
            # Receiver (cell j) expresses receptors
            receiver_receptor_expr = df_norm.loc[df_norm.index[j], receptor_genes].mean()
            
            # Communication score
            comm_score = sender_ligand_expr * receiver_receptor_expr
            communication_scores[i, j] = comm_score

print(f"  Communication scores computed")
print(f"  Mean score: {communication_scores[communication_scores > 0].mean():.4f}")
print(f"  Max score: {communication_scores.max():.4f}")

# ============================================================
# Step 4: Identify top sender-receiver pairs
# ============================================================

print("\nStep 4: Identifying top communication pairs")

# Flatten and find top pairs
pairs = []
for i in range(n_cells):
    for j in range(n_cells):
        if i != j and communication_scores[i, j] > 0:
            pairs.append({
                'Sender': df_norm.index[i],
                'Receiver': df_norm.index[j],
                'Score': communication_scores[i, j],
                'Sender_Type': cell_metadata.iloc[i]['Cell_Type'],
                'Receiver_Type': cell_metadata.iloc[j]['Cell_Type'],
                'Sender_IGF': cell_metadata.iloc[i]['IGF_Score'],
                'Receiver_IGF': cell_metadata.iloc[j]['IGF_Score']
            })

edges_df = pd.DataFrame(pairs)
edges_df = edges_df.sort_values('Score', ascending=False)
top_edges = edges_df.head(100)

print(f"  Total communication edges: {len(edges_df)}")
print(f"  Top 5 edges:")
for idx, row in top_edges.head(5).iterrows():
    print(f"    {row['Sender'][:20]} → {row['Receiver'][:20]}: {row['Score']:.4f}")

# ============================================================
# Step 5: Build communication network
# ============================================================

print("\nStep 5: Building communication network graph")

G = nx.DiGraph()

# Add nodes with attributes
for idx, row in cell_metadata.iterrows():
    G.add_node(row['Cell'], 
              igf_score=row['IGF_Score'],
              cell_type=row['Cell_Type'])

# Add edges
for idx, row in top_edges.iterrows():
    G.add_edge(row['Sender'], row['Receiver'], weight=row['Score'])

print(f"  Nodes: {G.number_of_nodes()}")
print(f"  Edges: {G.number_of_edges()}")

# Compute network metrics
in_degree = dict(G.in_degree())
out_degree = dict(G.out_degree())

print(f"  Mean in-degree: {np.mean(list(in_degree.values())):.2f}")
print(f"  Mean out-degree: {np.mean(list(out_degree.values())):.2f}")

# ============================================================
# Step 6: Save results
# ============================================================

print("\nStep 6: Saving results")

os.makedirs("results", exist_ok=True)

# Save edge list
edges_path = "results/LUAD_P01_communication_edges.csv"
edges_df.to_csv(edges_path, index=False)
print(f"  Saved: {edges_path} ({len(edges_df)} edges)")

# Save top 100 edges separately
top_edges_path = "results/LUAD_P01_communication_top100.csv"
top_edges.to_csv(top_edges_path, index=False)
print(f"  Saved: {top_edges_path}")

# Save communication scores (sparse format for memory)
comm_scores_path = "results/LUAD_P01_communication_scores.csv"
comm_dense = pd.DataFrame(communication_scores, 
                          index=df_norm.index, 
                          columns=df_norm.index)
# Save only non-zero entries
comm_sparse = comm_dense.stack()
comm_sparse = comm_sparse[comm_sparse > 0].reset_index()
comm_sparse.columns = ['Sender', 'Receiver', 'Score']
comm_sparse = comm_sparse.sort_values('Score', ascending=False)
comm_sparse.to_csv(comm_scores_path, index=False)
print(f"  Saved: {comm_scores_path}")

# Save sender-receiver summary
summary_path = "results/LUAD_P01_communication_summary.csv"
summary_data = []

for cell in df_norm.index:
    out_score = out_degree.get(cell, 0)
    in_score = in_degree.get(cell, 0)
    
    summary_data.append({
        'Cell': cell,
        'Cell_Type': cell_metadata[cell_metadata['Cell']==cell]['Cell_Type'].values[0],
        'IGF_Score': cell_metadata[cell_metadata['Cell']==cell]['IGF_Score'].values[0],
        'OutDegree_Sender': out_score,
        'InDegree_Receiver': in_score,
        'Activity_Score': out_score + in_score
    })

summary_df = pd.DataFrame(summary_data)
summary_df = summary_df.sort_values('Activity_Score', ascending=False)
summary_df.to_csv(summary_path, index=False)
print(f"  Saved: {summary_path}")

# ============================================================
# Step 7: Generate visualizations
# ============================================================

print("\nStep 7: Generating visualizations")

os.makedirs("plots", exist_ok=True)

# Plot 1: Communication score distribution
try:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Histogram of communication scores
    scores_nonzero = communication_scores[communication_scores > 0]
    axes[0, 0].hist(scores_nonzero, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
    axes[0, 0].set_xlabel('Communication Score')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Communication Scores')
    axes[0, 0].grid(alpha=0.3)
    
    # In-degree vs Out-degree
    cell_list = list(cell_metadata['Cell'])
    in_degrees = [in_degree.get(c, 0) for c in cell_list]
    out_degrees = [out_degree.get(c, 0) for c in cell_list]
    
    axes[0, 1].scatter(out_degrees, in_degrees, alpha=0.6, s=50, color='coral', edgecolor='black', linewidth=0.5)
    axes[0, 1].set_xlabel('Out-Degree (Sender)')
    axes[0, 1].set_ylabel('In-Degree (Receiver)')
    axes[0, 1].set_title('Cell Communication Profile')
    axes[0, 1].grid(alpha=0.3)
    
    # IGF vs Communication
    igf_list = [cell_metadata[cell_metadata['Cell']==c]['IGF_Score'].values[0] for c in cell_list]
    activity = [in_degree.get(c, 0) + out_degree.get(c, 0) for c in cell_list]
    
    axes[1, 0].scatter(igf_list, activity, alpha=0.6, s=50, color='green', edgecolor='black', linewidth=0.5)
    axes[1, 0].set_xlabel('IGF Score')
    axes[1, 0].set_ylabel('Communication Activity')
    axes[1, 0].set_title('IGF Score vs Communication')
    axes[1, 0].grid(alpha=0.3)
    
    # Top senders/receivers
    top_senders = summary_df.nlargest(10, 'OutDegree_Sender')[['Cell', 'OutDegree_Sender']]
    axes[1, 1].barh(range(len(top_senders)), top_senders['OutDegree_Sender'].values, color='steelblue', edgecolor='black')
    axes[1, 1].set_yticks(range(len(top_senders)))
    axes[1, 1].set_yticklabels([c[:15] for c in top_senders['Cell'].values], fontsize=8)
    axes[1, 1].set_xlabel('Out-Degree (Sender Activity)')
    axes[1, 1].set_title('Top 10 Sender Cells')
    axes[1, 1].grid(alpha=0.3, axis='x')
    
    plt.tight_layout()
    comm_plot = "plots/LUAD_P01_communication_analysis.png"
    plt.savefig(comm_plot, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {comm_plot}")
except Exception as e:
    print(f"  WARNING: Failed to save analysis plot: {e}")

# Plot 2: Network graph (sample)
try:
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Use only top 50 edges for visualization
    G_vis = nx.DiGraph()
    for idx, row in top_edges.head(50).iterrows():
        G_vis.add_edge(row['Sender'], row['Receiver'], weight=row['Score'])
    
    # Layout
    pos = nx.spring_layout(G_vis, k=2, iterations=50, seed=42)
    
    # Draw
    nx.draw_networkx_nodes(G_vis, pos, node_color='lightblue', 
                          node_size=300, ax=ax, edgecolors='black', linewidths=1)
    nx.draw_networkx_edges(G_vis, pos, edge_color='gray', alpha=0.5, 
                          arrows=True, arrowsize=10, arrowstyle='->', ax=ax, width=0.5)
    
    # Labels for high-activity nodes
    labels = {}
    for node in G_vis.nodes():
        deg = G_vis.degree(node)
        if deg > 3:
            labels[node] = node[:8]
    
    nx.draw_networkx_labels(G_vis, pos, labels, font_size=7, ax=ax)
    
    ax.set_title('LUAD_P01: Cell Communication Network (Top 50 Edges)', fontsize=12, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    
    network_plot = "plots/LUAD_P01_communication_network.png"
    plt.savefig(network_plot, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {network_plot}")
except Exception as e:
    print(f"  WARNING: Failed to save network plot: {e}")

# ============================================================
# Step 8: Validation
# ============================================================

print("\n" + "=" * 70)
print("VALIDATION")
print("=" * 70)
print(f"Communication analysis complete!")
print(f"  Total cells: {df_norm.shape[0]}")
print(f"  Communication edges: {len(edges_df)}")
print(f"  Top sender: {summary_df.iloc[0]['Cell'][:30]} (activity: {summary_df.iloc[0]['Activity_Score']:.0f})")
print(f"  Mean activity per cell: {summary_df['Activity_Score'].mean():.2f}")
print(f"  Files saved: 4 CSV + 2 PNG")
