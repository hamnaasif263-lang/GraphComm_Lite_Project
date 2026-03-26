import os
import glob
import logging
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)

# Existing IGF ligand/receptor helper (kept for compatibility)
ligand_receptor_pairs = [
    ("IGF1", "IGF1R"),
    ("IGF2", "IGF1R"),
    ("IGF1", "INSR")
]


def build_igf_graph(df, threshold=0.2):
    """Build IGF ligand-receptor network.
    
    If no IGF genes are present, falls back to k-NN based on top genes.
    """
    n = len(df)
    A = np.zeros((n, n))
    found_pairs = 0
    
    for ligand, receptor in ligand_receptor_pairs:
        if ligand not in df.columns or receptor not in df.columns:
            continue
        found_pairs += 1
        for i in range(n):
            for j in range(n):
                if df.iloc[i][ligand] > threshold and df.iloc[j][receptor] > threshold:
                    A[i, j] += 1
    
    # If no IGF ligand-receptor pairs found, build graph from k-NN of top genes
    if found_pairs == 0:
        logging.warning("No IGF ligand-receptor pairs found in data. Building from top genes.")
        # Use top 10 expressed genes to build a similarity graph
        top_genes = df.mean(axis=0).nlargest(10).index.tolist()
        if len(top_genes) > 0:
            from sklearn.metrics.pairwise import cosine_similarity
            from sklearn.preprocessing import StandardScaler
            X = StandardScaler().fit_transform(df[top_genes].values)
            sim = cosine_similarity(X)
            # Set diagonal to 0 and use top 8 neighbors per cell
            np.fill_diagonal(sim, 0)
            for i in range(n):
                top_neighbors = np.argsort(sim[i])[-8:]  # 8 nearest neighbors
                A[i, top_neighbors] = 1
    
    return (A > 0).astype(float)


# --- New pathway analysis utilities ---
def load_pathway_genes(pathway_name: str, pathways_dir: str = "data/pathways") -> List[str]:
    path = os.path.join(pathways_dir, f"{pathway_name}.txt")
    if not os.path.exists(path):
        logging.warning("Pathway file not found: %s", path)
        return []
    with open(path, "r", encoding="utf-8") as fh:
        genes = [l.strip() for l in fh if l.strip()]
    return genes


def compute_pathway_score(df: pd.DataFrame, genes: List[str], method: str = "mean") -> pd.Series:
    """Compute per-cell pathway activity score.

    method: 'mean' (default) computes mean expression across pathway genes per cell.
            'zscore' computes mean expression then z-score normalizes the resulting series.
    """
    if len(genes) == 0:
        return pd.Series([], dtype=float)
    genes_present = [g for g in genes if g in df.columns]
    if len(genes_present) == 0:
        logging.warning("None of pathway genes present in expression dataframe: %s", genes[:5])
        # return series of zeros with same index
        return pd.Series(0.0, index=df.index)

    expr = df[genes_present].astype(float)
    score = expr.mean(axis=1)

    if method == "zscore":
        mu = score.mean()
        sigma = score.std()
        if sigma == 0 or np.isnan(sigma):
            logging.info("Score std is zero; returning zeros for zscore normalization")
            return pd.Series(0.0, index=score.index)
        return (score - mu) / sigma

    return score


def normalize_minmax(series: pd.Series, eps: float = 1e-9) -> pd.Series:
    """Min-max normalize a pandas Series to [0,1], safely handling zero-range."""
    s_min = series.min()
    s_max = series.max()
    denom = (s_max - s_min)
    if denom == 0 or np.isnan(denom):
        return pd.Series(0.0, index=series.index)
    return (series - s_min) / (denom + eps)


def plot_pathway_scores(scores: pd.Series, pathway: str, sample_name: str, out_path: str):
    """Save a compact heatmap-like plot and a per-cell scatter plot (saved as a single PNG)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    try:
        fig, axes = plt.subplots(2, 1, figsize=(6, 5), gridspec_kw={"height_ratios": [1, 2]})

        # Heatmap: vertical bar of scores
        im = axes[0].imshow(scores.values.reshape(1, -1), aspect="auto", cmap="viridis")
        axes[0].set_yticks([])
        axes[0].set_xticks([])
        axes[0].set_title(f"{pathway} activity — {sample_name}")
        fig.colorbar(im, ax=axes[0], orientation="horizontal", fraction=0.05)

        # Scatter per cell
        axes[1].scatter(range(len(scores)), scores.values, s=6, alpha=0.6)
        axes[1].set_xlabel("Cell (ordered)")
        axes[1].set_ylabel("Pathway score")

        plt.tight_layout()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
    except Exception as e:
        logging.warning("Failed to plot pathway %s for %s: %s", pathway, sample_name, e)


def analyze_sample_pathways(
    sc_path: str,
    pathways: List[str],
    pathways_dir: str = "data/pathways",
    results_dir: str = "results",
    plots_dir: str = "plots",
    score_method: str = "mean",
):
    """Compute pathway scores for a single scRNA CSV file and save per-pathway CSVs and plots.

    Also produces a summary CSV with all pathway scores per cell: `results/<sample>_pathway_scores.csv`.
    """
    logging.info("Analyzing sample: %s", sc_path)
    df = pd.read_csv(sc_path, index_col=0)
    sample_basename = os.path.splitext(os.path.basename(sc_path))[0]

    all_scores = pd.DataFrame(index=df.index)

    for pathway in pathways:
        genes = load_pathway_genes(pathway, pathways_dir=pathways_dir)
        score = compute_pathway_score(df, genes, method=score_method)
        norm = normalize_minmax(score)

        # Save per-pathway CSV
        os.makedirs(results_dir, exist_ok=True)
        per_path_csv = os.path.join(results_dir, f"{pathway}_{sample_basename}.csv")
        out_df = pd.DataFrame({"cell": df.index, "score": score.values, "score_norm": norm.values})
        out_df.to_csv(per_path_csv, index=False)
        logging.info("Wrote %s", per_path_csv)

        # Save plot
        os.makedirs(plots_dir, exist_ok=True)
        per_path_plot = os.path.join(plots_dir, f"{pathway}_{sample_basename}.png")
        plot_pathway_scores(norm, pathway, sample_basename, per_path_plot)
        logging.info("Wrote plot %s", per_path_plot)

        # Add to summary
        all_scores[pathway] = score

    # Save summary CSV
    summary_csv = os.path.join(results_dir, f"{sample_basename}_pathway_scores.csv")
    all_scores.insert(0, "cell", all_scores.index)
    all_scores.to_csv(summary_csv, index=False)
    logging.info("Wrote pathway summary %s", summary_csv)


def analyze_all_samples(
    sc_glob: str = "data/scRNA_*.csv",
    pathways: List[str] = None,
    **kwargs,
):
    if pathways is None:
        # default list that matches the user request
        pathways = [
            "MAPK",
            "JAK_STAT",
            "TGF_beta",
            "Hypoxia_HIF1A",
            "CellCycle_G1S",
            "CellCycle_G2M",
            "EMT",
            "Apoptosis",
            "p53",
            "Wnt",
            "Notch",
            "Hedgehog",
        ]

    files = sorted(glob.glob(sc_glob))
    logging.info("Found %d samples to analyze", len(files))
    for f in files:
        analyze_sample_pathways(f, pathways, **kwargs)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compute pathway activity per sample and save results/plots")
    parser.add_argument("--sample", type=str, help="Path to a single scRNA CSV (cell x gene). If omitted, runs on all data/scRNA_*.csv", default=None)
    parser.add_argument("--method", type=str, choices=["mean", "zscore"], default="mean")
    parser.add_argument("--results-dir", type=str, default="results")
    parser.add_argument("--plots-dir", type=str, default="plots")
    args = parser.parse_args()

    pathways_default = [
        "MAPK",
        "JAK_STAT",
        "TGF_beta",
        "Hypoxia_HIF1A",
        "CellCycle_G1S",
        "CellCycle_G2M",
        "EMT",
        "Apoptosis",
        "p53",
        "Wnt",
        "Notch",
        "Hedgehog",
    ]

    if args.sample:
        analyze_sample_pathways(args.sample, pathways_default, score_method=args.method, results_dir=args.results_dir, plots_dir=args.plots_dir)
    else:
        analyze_all_samples(pathways=pathways_default, score_method=args.method, results_dir=args.results_dir, plots_dir=args.plots_dir)
