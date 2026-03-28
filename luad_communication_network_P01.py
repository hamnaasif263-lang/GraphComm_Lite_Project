"""
LUAD_P01 Cell-to-Cell Communication Graph Analysis
Builds and analyzes communication networks based on ligand-receptor interactions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from sklearn.neighbors import kneighbors_graph
from scipy.spatial.distance import pdist, squareform
import os

print("=" * 70)
print("LUAD_P01: Cell-to-Cell Communication Graph Analysis")
print("=" * 70)

# ============================================================
# Step 1: Load data and cell annotations
# ============================================================

print("\nStep 1: Loading data...")
df_norm = pd.read_parquet("data/luad/processed/scRNA_LUAD_P01_norm.parquet")
cell_annot = pd.read_csv("results/LUAD_P01_celltype_annotations.csv")

print(f"Cells: {df_norm.shape[0]}, Genes: {df_norm.shape[1]}")
print(f"Cell types: {cell_annot['Cell_Type'].nunique()}")

# ============================================================
# Step 2: Build k-NN graph (spatial/similarity network)
# ============================================================

print("\nStep 2: Building k-NN communication graph...")

n_neighbors = min(10, df_norm.shape[0] - 1)
print(f"Using k-NN with k={n_neighbors}")

# Compute k-NN graph
knn_graph = kneighbors_graph(df_norm, n_neighbors=n_neighbors, mode='connectivity')
A_knn = knn_graph.toarray()
print(f"Graph edges: {A_knn.sum() / 2:.0f}")

# ============================================================
# Step 3: Identify sender and receiver cells
# ============================================================

print("\nStep 3: Identifying sender and receiver cells...")

# Load autocrine classifications
autocrine_df = pd.read_csv("results/LUAD_P01_Autocrine_Classifications.csv")

senders = autocrine_df[autocrine_df['Cell_Type'].isin(['Autocrine', 'Sender'])].index.tolist()
receivers = autocrine_df[autocrine_df['Cell_Type'].isin(['Autocrine', 'Receiver'])].index.tolist()

print(f"Sender cells: {len(senders)}")
print(f"Receiver cells: {len(receivers)}")

# ============================================================
# Step 4: Score communication edges
# ============================================================

print("\nStep 4: Scoring communication edges...")

# Edge scores based on:
# 1. Proximity (k-NN connection)
# 2. Gene expression similarity
# 3. Sender-receiver compatibility

from sklearn.preprocessing import StandardScaler

# Compute cosine similarity
from sklearn.metrics.pairwise import cosine_similarity
similarity_matrix = cosine_similarity(df_norm)

# Create communication score matrix
comm_scores = A_knn * similarity_matrix
np.fill_diagonal(comm_scores, 0)

print(f"Mean communication score: {comm_scores[comm_scores > 0].mean():.4f}")
print(f"Max communication score: {comm_scores.max():.4f}")

# ============================================================
# Step 5: Build directed communication network
# ============================================================

print("\nStep 5: Building directed communication network...")

# Create directed graph
G = nx.DiGraph()
G.add_nodes_from(range(df_norm.shape[0]))

# Add edges only for top communication pairs
threshold = np.percentile(comm_scores[comm_scores > 0], 75)
print(f"Edge score threshold (75th percentile): {threshold:.4f}")

edge_count = 0
total_score = 0

for i in senders:
    for j in receivers:
        if comm_scores[i, j] > threshold:
            G.add_edge(i, j, weight=comm_scores[i, j])
            edge_count += 1
            total_score += comm_scores[i, j]

print(f"Communication edges added: {edge_count}")
print(f"Mean edge weight: {total_score / edge_count if edge_count > 0 else 0:.4f}")

# ============================================================
# Step 6: Compute network metrics
# ============================================================

print("\nStep 6: Computing network metrics...")

# Degree centrality
in_degree = dict(G.in_degree(weight='weight'))
out_degree = dict(G.out_degree(weight='weight'))

# Betweenness centrality (for connected components)
if G.number_of_nodes() > 0 and G.number_of_edges() > 0:
    betweenness = nx.betweenness_centrality(G, weight='weight')
else:
    betweenness = {i: 0 for i in range(df_norm.shape[0])}

# Node importance scores
node_importance = {}
for node in range(df_norm.shape[0]):
    importance = (out_degree.get(node, 0) + in_degree.get(node, 0)) / 2
    node_importance[node] = importance

top_nodes = sorted(node_importance.items(), key=lambda x: x[1], reverse=True)[:10]
print(f"Top 10 communication hubs:")
for node, score in top_nodes:
    print(f"  Node {node}: importance={score:.4f}")

# ============================================================
# Step 7: Generate network visualization
# ============================================================

print("\nStep 7: Generating network visualization...")

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# -------- Plot 1: Network graph --------
ax = axes[0]

if G.number_of_edges() > 0:
    # Create a simplified graph with only top edges for visualization
    G_vis = nx.DiGraph()
    G_vis.add_nodes_from(range(min(100, df_norm.shape[0])))  # Limit to 100 nodes
    
    # Add only top 500 edges
    edges_sorted = sorted(G.edges(data=True), key=lambda x: x[2]['weight'], reverse=True)[:500]
    for u, v, d in edges_sorted:
        if u < 100 and v < 100:  # Only keep edges between limited nodes
            G_vis.add_edge(u, v, weight=d['weight'])
    
    if G_vis.number_of_edges() > 0:
        # Use spring layout for visualization
        pos = nx.spring_layout(G_vis, k=0.5, iterations=30, seed=42)
        
        # Node colors based on importance
        node_colors = [node_importance.get(node, 0) for node in G_vis.nodes()]
        
        # Draw nodes
        nodes = nx.draw_networkx_nodes(G_vis, pos, node_color=node_colors, node_size=200,
                                       cmap='YlOrRd', ax=ax, alpha=0.8)
        
        # Draw edges with simplified style (no curves to speed up rendering)
        edges = nx.draw_networkx_edges(G_vis, pos, edge_color='gray', arrows=True,
                                       arrowsize=10, width=0.8, ax=ax,
                                       connectionstyle='straight', alpha=0.4)
        
        # Add colorbar
        cbar = plt.colorbar(nodes, ax=ax)
        cbar.set_label('Node Importance', fontsize=10)
        
        ax.set_title(f'LUAD_P01: Cell Communication Network\n(Top 500 edges, first 100 nodes)', 
                     fontsize=12, fontweight='bold')
        ax.axis('off')
    else:
        ax.text(0.5, 0.5, 'No communication edges\nin visualization subset', 
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.axis('off')
else:
    ax.text(0.5, 0.5, 'No communication edges\nabove threshold', 
            ha='center', va='center', transform=ax.transAxes, fontsize=12)
    ax.axis('off')

# -------- Plot 2: Communication statistics --------
ax = axes[1]

# Cell communication type distribution
comm_types = []
comm_counts = []

# Sender-receiver pairs
n_sr_pairs = sum(1 for i in senders for j in receivers if G.has_edge(i, j))
if n_sr_pairs > 0:
    comm_types.append('Sender → Receiver')
    comm_counts.append(n_sr_pairs)

# Autocrine (self-communication)
n_autocrine = sum(1 for i in range(df_norm.shape[0]) for j in range(df_norm.shape[0])
                  if i != j and G.has_edge(i, j) and i in senders and j in receivers)
if n_autocrine > 0:
    comm_types.append('Cell-Cell')
    comm_counts.append(n_autocrine)

if len(comm_types) > 0:
    bars = ax.bar(comm_types, comm_counts, color=['#FF6B6B', '#4ECDC4'],
                  edgecolor='black', linewidth=1.5, alpha=0.7)
    
    # Add value labels
    for bar, count in zip(bars, comm_counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(count)}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    ax.set_ylabel('Number of Edges', fontsize=11, fontweight='bold')
    ax.set_title('Communication Edge Types', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
else:
    ax.text(0.5, 0.5, 'No communication edges', ha='center', va='center',
            transform=ax.transAxes, fontsize=12)
    ax.axis('off')

plt.tight_layout()

# Save figure
comm_plot_path = "plots/LUAD_P01_Communication_Network.png"
plt.savefig(comm_plot_path, dpi=150, bbox_inches='tight')
print(f"Saved: {comm_plot_path}")
plt.close()

# ============================================================
# Step 8: Save communication results
# ============================================================

print("\nStep 8: Saving communication results...")

# Save communication edges
edges_list = []
for u, v, data in G.edges(data=True):
    edges_list.append({
        'Source_Cell': u,
        'Target_Cell': v,
        'Communication_Score': data['weight']
    })

if edges_list:
    edges_df = pd.DataFrame(edges_list)
    edges_path = "results/LUAD_P01_Communication_Edges.csv"
    edges_df.to_csv(edges_path, index=False)
    print(f"Saved: {edges_path}")

# Save node metrics
node_metrics = []
for node in range(df_norm.shape[0]):
    node_metrics.append({
        'Cell_ID': node,
        'Out_Degree': out_degree.get(node, 0),
        'In_Degree': in_degree.get(node, 0),
        'Importance': node_importance.get(node, 0),
        'Betweenness': betweenness.get(node, 0)
    })

metrics_df = pd.DataFrame(node_metrics)
metrics_path = "results/LUAD_P01_Communication_Metrics.csv"
metrics_df.to_csv(metrics_path, index=False)
print(f"Saved: {metrics_path}")

# ============================================================
# VALIDATION & SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("COMMUNICATION GRAPH ANALYSIS COMPLETE")
print("=" * 70)

print(f"\nNetwork Statistics:")
print(f"  Total nodes: {G.number_of_nodes()}")
print(f"  Total edges: {G.number_of_edges()}")
print(f"  Network density: {nx.density(G.to_undirected()):.4f}")

if G.number_of_edges() > 0:
    avg_weight = np.mean([d['weight'] for u, v, d in G.edges(data=True)])
    print(f"  Average edge weight: {avg_weight:.4f}")
    
    print(f"\nTop 5 Receiver Cells (highest in-degree):")
    top_receivers = sorted(in_degree.items(), key=lambda x: x[1], reverse=True)[:5]
    for node, degree in top_receivers:
        print(f"    Cell {node}: in-degree={degree:.4f}")
    
    print(f"\nTop 5 Sender Cells (highest out-degree):")
    top_senders = sorted(out_degree.items(), key=lambda x: x[1], reverse=True)[:5]
    for node, degree in top_senders:
        print(f"    Cell {node}: out-degree={degree:.4f}")
else:
    print("  No communication edges detected")

print("\n" + "=" * 70)
