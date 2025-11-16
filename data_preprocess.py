# ==========================================================
# data_preprocess.py
# ----------------------------------------------------------
# Handles scRNA-seq preprocessing for GraphComm-Lite:
#   1. Loads CSV expression data
#   2. Scales it
#   3. Reduces dimensionality with PCA (auto-adjusted)
#   4. Builds a k-NN graph
#   5. Extracts IGF pathway expression
# ==========================================================

import numpy as np
import pandas as pd
import networkx as nx
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import kneighbors_graph

# ----------------------------------------------------------
# Function: preprocess_scRNA
# ----------------------------------------------------------
def preprocess_scRNA(filepath: str, n_components: int = 50, n_neighbors: int = 8):
    """
    Load and preprocess scRNA-seq data.
    Returns reduced PCA features, adjacency matrix, and graph object.
    Automatically adjusts PCA components for small datasets.
    """
    try:
        df = pd.read_csv(filepath, index_col=0)
    except Exception:
        df = pd.read_csv(filepath)  # try without index if above fails

    print(f"🧬 Loaded scRNA-seq data with shape: {df.shape}")

    # If dataset too small for PCA, limit n_components safely
    max_allowed = min(df.shape[0], df.shape[1])
    n_components = min(n_components, max_allowed)
    if n_components < 2:
        print("⚠️ Dataset too small; using no PCA reduction.")
        n_components = 1

    # Normalize and reduce dimensions
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(df)

    if n_components > 1:
        pca = PCA(n_components=n_components)
        data_pca = pca.fit_transform(data_scaled)
        print(f"✅ PCA reduced data to {n_components} components.")
    else:
        data_pca = data_scaled
        print("✅ Skipped PCA reduction (tiny dataset).")

    # Build adjacency matrix using k-NN
    n_neighbors = min(n_neighbors, df.shape[0] - 1)
    A = kneighbors_graph(data_pca, n_neighbors=n_neighbors, mode='connectivity').toarray()
    G = nx.from_numpy_array(A)

    print(f"🔗 Graph built with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    return data_pca, A, G


# ----------------------------------------------------------
# Function: extract_igf_expression
# ----------------------------------------------------------
def extract_igf_expression(df: pd.DataFrame, igf_gene_list: list):
    """
    Extracts average IGF pathway activity per cell.
    Returns a pandas Series.
    """
    found = [g for g in igf_gene_list if g in df.columns]
    if not found:
        print("⚠️ No IGF genes found in dataset columns.")
        return pd.Series(np.zeros(df.shape[0]), index=df.index)

    igf_expr = df[found].mean(axis=1)
    print(f"✅ Extracted IGF expression using {len(found)} genes: {found}")
    return igf_expr
