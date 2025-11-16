# ==========================================================
# visualize.py
# ----------------------------------------------------------
# Visualization tools for GraphComm-Lite:
#   - IGF pathway activity scatter plot
#   - Communication score heatmap
#   - Network graph of cell communication
# ==========================================================

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import networkx as nx


# ----------------------------------------------------------
def plot_igf_activity(igf_values, save_path: str = None):
    """
    Scatter plot of IGF pathway activity.
    """
    plt.figure(figsize=(8, 4))
    plt.scatter(range(len(igf_values)), igf_values, s=20)
    plt.title("IGF Pathway Activity per Cell")
    plt.xlabel("Cell Index")
    plt.ylabel("IGF Activity Score")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        plt.close()
    else:
        plt.show()


# ----------------------------------------------------------
def plot_dormancy_heatmap(scores, save_path: str = None):
    """
    Heatmap of GraphComm-Lite risk (dormancy) scores.
    """
    plt.figure(figsize=(10, 1.5))
    sns.heatmap([scores], cmap="viridis", cbar=True)
    plt.title("Dormancy / Communication Scores")
    plt.xlabel("Cell Index")
    plt.yticks([], [])
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        plt.close()
    else:
        plt.show()


# ----------------------------------------------------------
def plot_communication_graph(G, scores, save_path: str = None):
    """
    Visualize the communication graph with node colors based on risk scores.
    """
    plt.figure(figsize=(8, 8))

    # Normalize scores for color mapping
    norm_scores = (scores - np.min(scores)) / (np.max(scores) - np.min(scores))

    pos = nx.spring_layout(G, seed=42)  # nice stable layout

    nx.draw(
        G,
        pos,
        node_color=norm_scores,
        cmap="coolwarm",
        node_size=120,
        edge_color="gray",
        linewidths=0.2,
    )

    plt.title("Cell-to-Cell IGF Communication Network")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        plt.close()
    else:
        plt.show()
