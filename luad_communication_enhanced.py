import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
from matplotlib.patches import FancyBboxPatch
import warnings
warnings.filterwarnings('ignore')

# Load data
print("Loading communication network data...")
comm_edges = pd.read_csv("results/LUAD_P01_Communication_Edges.csv")
comm_metrics = pd.read_csv("results/LUAD_P01_Communication_Metrics.csv")
autocrine_class = pd.read_csv("results/LUAD_P01_Autocrine_Classifications.csv")

# Build network
print("Building directed network...")
G = nx.DiGraph()

# Add all nodes
for cell_id in range(1308):
    G.add_node(cell_id)

# Add weighted edges
for _, row in comm_edges.iterrows():
    source = int(row['Source_Cell'])
    target = int(row['Target_Cell'])
    weight = row['Communication_Score']
    G.add_edge(source, target, weight=weight)

print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Create comprehensive visualization
fig = plt.figure(figsize=(20, 12))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

# ============ SUBPLOT 1: MAIN NETWORK ============
ax1 = fig.add_subplot(gs[0, :])

# Use spring layout with larger k for better spacing
print("Computing layout...")
pos = nx.spring_layout(G, k=0.5, iterations=20, seed=42)

# Get node properties
in_degree = dict(G.in_degree())
out_degree = dict(G.out_degree())

# Create color map based on autocrine classification
node_colors = []
color_map = {
    'Autocrine (Both +)': '#FF6B6B',
    'Sender (Ligand only)': '#4ECDC4',
    'Receiver (Receptor only)': '#45B7D1',
    'Neither': '#95A5A6'
}

for node_id in G.nodes():
    cell_row = autocrine_class[autocrine_class['Cell_ID'] == node_id]
    if len(cell_row) > 0:
        cell_type = cell_row.iloc[0]['Cell_Type']
        node_colors.append(color_map.get(cell_type, '#95A5A6'))
    else:
        node_colors.append('#95A5A6')

# Draw network
node_sizes = [in_degree[node] * 50 + 100 for node in G.nodes()]

nx.draw_networkx_nodes(
    G, pos,
    node_color=node_colors,
    node_size=node_sizes,
    alpha=0.7,
    ax=ax1,
    edgecolors='black',
    linewidths=0.5
)

# Draw edges with varying thickness
edges = G.edges()
weights = [G[u][v]['weight'] for u, v in edges]
min_weight, max_weight = min(weights), max(weights)

nx.draw_networkx_edges(
    G, pos,
    width=[1 + 2*((w-min_weight)/(max_weight-min_weight)) for w in weights],
    alpha=0.2,
    edge_color='gray',
    ax=ax1,
    connectionstyle='arc3,rad=0.1',
    arrowsize=10,
    arrowstyle='->'
)

# Label top hub nodes
top_hubs = sorted(in_degree.items(), key=lambda x: x[1], reverse=True)[:5]
for node_id, in_deg in top_hubs:
    x, y = pos[node_id]
    ax1.text(x, y, f"C{node_id}", fontsize=9, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))

ax1.set_title('LUAD_P01 Cell-to-Cell Communication Network\n(1,308 cells, 2,446 communication edges)', 
              fontsize=14, fontweight='bold', pad=20)
ax1.axis('off')

# Legend
legend_elements = [
    mpatches.Patch(facecolor='#FF6B6B', edgecolor='black', label='Autocrine (31%)'),
    mpatches.Patch(facecolor='#4ECDC4', edgecolor='black', label='Sender (13.9%)'),
    mpatches.Patch(facecolor='#45B7D1', edgecolor='black', label='Receiver (13.6%)'),
    mpatches.Patch(facecolor='#95A5A6', edgecolor='black', label='Inactive (41.5%)')
]
ax1.legend(handles=legend_elements, loc='upper left', fontsize=10, framealpha=0.95)

# ============ SUBPLOT 2: IN-DEGREE DISTRIBUTION ============
ax2 = fig.add_subplot(gs[1, 0])

in_degrees = list(in_degree.values())
ax2.hist(in_degrees, bins=30, color='#4ECDC4', alpha=0.7, edgecolor='black')
ax2.axvline(np.mean(in_degrees), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(in_degrees):.2f}')
ax2.axvline(np.median(in_degrees), color='orange', linestyle='--', linewidth=2, label=f'Median: {np.median(in_degrees):.2f}')
ax2.set_xlabel('In-Degree (Receiver Capability)', fontsize=11, fontweight='bold')
ax2.set_ylabel('Number of Cells', fontsize=11, fontweight='bold')
ax2.set_title('Receiver Cell Distribution\n(How many cells target each cell)', fontsize=12, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

# ============ SUBPLOT 3: OUT-DEGREE DISTRIBUTION ============
ax3 = fig.add_subplot(gs[1, 1])

out_degrees = list(out_degree.values())
ax3.hist(out_degrees, bins=30, color='#45B7D1', alpha=0.7, edgecolor='black')
ax3.axvline(np.mean(out_degrees), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(out_degrees):.2f}')
ax3.axvline(np.median(out_degrees), color='orange', linestyle='--', linewidth=2, label=f'Median: {np.median(out_degrees):.2f}')
ax3.set_xlabel('Out-Degree (Sender Capability)', fontsize=11, fontweight='bold')
ax3.set_ylabel('Number of Cells', fontsize=11, fontweight='bold')
ax3.set_title('Sender Cell Distribution\n(How many cells each cell targets)', fontsize=12, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)

# Add network statistics text box
stats_text = (
    f"Network Statistics:\n"
    f"  • Total Nodes: {G.number_of_nodes()}\n"
    f"  • Total Edges: {G.number_of_edges()}\n"
    f"  • Network Density: {nx.density(G):.4f}\n"
    f"  • Avg In-Degree: {np.mean(in_degrees):.2f}\n"
    f"  • Avg Out-Degree: {np.mean(out_degrees):.2f}\n"
    f"  • Max In-Degree: {max(in_degrees):.2f}\n"
    f"  • Max Out-Degree: {max(out_degrees):.2f}"
)

fig.text(0.98, 0.02, stats_text, fontsize=10, family='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
         ha='right', va='bottom')

plt.savefig('results/LUAD_P01_Communication_Network.png', dpi=300, bbox_inches='tight')
print("✓ Network visualization saved: results/LUAD_P01_Communication_Network.png")

# Generate additional heatmap: top 50 hub cells communication matrix
print("\nGenerating hub cell communication heatmap...")

fig2, ax = plt.subplots(figsize=(14, 12))

# Get top 50 hub cells
top_50_hubs = sorted(in_degree.items(), key=lambda x: x[1], reverse=True)[:50]
hub_ids = [x[0] for x in top_50_hubs]

# Create adjacency matrix for top hubs
adj_matrix = np.zeros((50, 50))
for i, source_id in enumerate(hub_ids):
    for j, target_id in enumerate(hub_ids):
        if G.has_edge(source_id, target_id):
            adj_matrix[i][j] = G[source_id][target_id]['weight']

# Plot heatmap
im = ax.imshow(adj_matrix, cmap='YlOrRd', aspect='auto')
ax.set_xticks(range(50))
ax.set_yticks(range(50))
ax.set_xticklabels([f'C{hub_ids[i]}' for i in range(50)], rotation=90, fontsize=8)
ax.set_yticklabels([f'C{hub_ids[i]}' for i in range(50)], fontsize=8)
ax.set_xlabel('Target Cell', fontsize=12, fontweight='bold')
ax.set_ylabel('Source Cell', fontsize=12, fontweight='bold')
ax.set_title('Top 50 Hub Cells - Communication Interaction Matrix\n(Communication scores between hub cells)', 
             fontsize=13, fontweight='bold', pad=15)

# Add colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Communication Score', fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig('results/LUAD_P01_Communication_Hubs_Matrix.png', dpi=300, bbox_inches='tight')
print("✓ Hub cell heatmap saved: results/LUAD_P01_Communication_Hubs_Matrix.png")

plt.close('all')

# Generate summary statistics
print("\n" + "="*60)
print("COMMUNICATION NETWORK ANALYSIS SUMMARY")
print("="*60)

print(f"\nNetwork Size:")
print(f"  Nodes (cells): {G.number_of_nodes()}")
print(f"  Edges (connections): {G.number_of_edges()}")
print(f"  Density: {nx.density(G):.4f}")

print(f"\nNode Degree Statistics:")
print(f"  In-Degree (receivers):")
print(f"    Mean: {np.mean(in_degrees):.2f}")
print(f"    Median: {np.median(in_degrees):.2f}")
print(f"    Max: {max(in_degrees):.2f} (Cell {max(in_degree, key=in_degree.get)})")
print(f"  Out-Degree (senders):")
print(f"    Mean: {np.mean(out_degrees):.2f}")
print(f"    Median: {np.median(out_degrees):.2f}")
print(f"    Max: {max(out_degrees):.2f} (Cell {max(out_degree, key=out_degree.get)})")

print(f"\nTop 10 Receiver Hubs:")
for i, (cell_id, in_deg) in enumerate(top_50_hubs[:10], 1):
    out_deg = out_degree[cell_id]
    print(f"  {i:2d}. Cell {cell_id:4d}: In={in_deg:6.2f}, Out={out_deg:6.2f}")

print(f"\nTop 10 Sender Hubs (by out-degree):")
top_senders = sorted(out_degree.items(), key=lambda x: x[1], reverse=True)[:10]
for i, (cell_id, out_deg) in enumerate(top_senders, 1):
    in_deg = in_degree[cell_id]
    print(f"  {i:2d}. Cell {cell_id:4d}: Out={out_deg:6.2f}, In={in_deg:6.2f}")

print("\n✓ All communication plots generated successfully!")
print("="*60)
