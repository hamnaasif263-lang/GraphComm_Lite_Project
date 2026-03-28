#!/usr/bin/env python3
"""Process all per-sample scRNA files in `data/` one-by-one.

Behavior:
- Converts any 10x HDF5 files found via `scanpy.read_10x_h5` to CSVs (saved as `data/scRNA_<name>.csv`).
- Looks for CSVs named `scRNA_*.csv` and processes each with the pipeline functions (preprocess -> IGF graph -> train GNN -> RF demo).
- Saves per-sample results to `results/<sample>_results.csv` and plots to `plots/` with the sample prefix.

Usage: drop your per-sample H5 or CSV files into `data/` then run this script.
"""
from pathlib import Path
import os
import sys
import subprocess

import numpy as np
import pandas as pd

from data_preprocess import preprocess_scRNA, extract_igf_expression
from pathway_analysis import build_igf_graph
from train_predict import train_graph_model, drug_response_prediction
from visualize import plot_igf_activity, plot_dormancy_heatmap, plot_communication_graph

DATA_DIR = Path("data")
OUT_DIR = Path("results")
PLOTS_DIR = Path("plots")
OUT_DIR.mkdir(exist_ok=True)
PLOTS_DIR.mkdir(exist_ok=True)


def convert_h5_if_needed(path_h5: Path):
    """Convert a 10x H5 to CSV using existing helper if present.
    Returns path to CSV (either existing or newly created).
    """
    csv_out = DATA_DIR / f"scRNA_{path_h5.stem}.csv"
    # attempt to use download_and_convert.ann_to_csv_h5 if available
    try:
        from download_and_convert import ann_to_csv_h5
        print(f"Converting H5 {path_h5} -> {csv_out}")
        ann_to_csv_h5(path_h5, csv_out)
        return csv_out
    except Exception as e:
        print(f"Could not convert {path_h5} automatically: {e}")
        return None


def process_csv(csv_path: Path, sample_name: str):
    print(f"\n--- Processing sample {sample_name}: {csv_path} ---")
    # 1) Preprocess
    features_np, A_knn, G_knn = preprocess_scRNA(str(csv_path))

    # 2) Read expression DF for IGF graph
    try:
        df_expr = pd.read_csv(csv_path, index_col=0)
    except Exception:
        df_expr = pd.read_csv(csv_path)

    # 3) Build IGF graph
    A_igf = build_igf_graph(df_expr, threshold=0.25)

    # 4) Combine with snRNA features (lightweight)
    try:
        from data import snRNA_U1
    except Exception:
        # fallback simple base freq
        snrna_vec = np.array([0.25, 0.25, 0.25, 0.25], dtype=float)
    try:
        snrna_tiled = np.tile(snrna_vec, (features_np.shape[0], 1))
        features_combined = np.hstack([features_np, snrna_tiled])
    except Exception:
        features_combined = features_np

    # 5) Prepare tensors and labels
    import torch
    features_t = torch.tensor(features_combined, dtype=torch.float32)
    adj_t = torch.tensor(A_igf, dtype=torch.float32)

    if "label" in df_expr.columns:
        raw_labels = df_expr["label"].astype(int)
    else:
        # Try to generate labels from IGF pathway genes
        igf_genes = ["IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1"]
        available_genes = [g for g in igf_genes if g in df_expr.columns]
        
        if available_genes:
            # Use available IGF pathway genes
            igf_series = df_expr[available_genes].mean(axis=1)
            thresh = igf_series.median()
            raw_labels = (igf_series > thresh).astype(int)
            print(f"✅ Generated labels from {len(available_genes)} IGF genes: {available_genes}")
        else:
            # Fallback: use highly variable genes or random stratified split
            print("⚠️ No IGF genes found. Using highly expressed genes for label generation.")
            gene_means = df_expr.mean(axis=0).sort_values(ascending=False)
            top_genes = gene_means.head(10).index.tolist()
            if len(top_genes) > 0:
                expr_score = df_expr[top_genes].mean(axis=1)
                thresh = expr_score.median()
                raw_labels = (expr_score > thresh).astype(int)
            else:
                # Final fallback: 50-50 random split
                print("⚠️ Using random 50-50 split for labels.")
                raw_labels = pd.Series(np.random.randint(0, 2, size=df_expr.shape[0]), index=df_expr.index)
        
        # Validate label distribution
        n_pos = (raw_labels == 1).sum()
        n_neg = (raw_labels == 0).sum()
        print(f"📊 Label distribution: {n_pos} positive ({100*n_pos/(n_pos+n_neg):.1f}%), {n_neg} negative")
        
        if n_pos == 0 or n_neg == 0:
            print("⚠️ Imbalanced labels detected! Rebalancing with stratified split...")
            raw_labels = pd.Series(np.random.randint(0, 2, size=df_expr.shape[0]), index=df_expr.index)
            n_pos = (raw_labels == 1).sum()
            n_neg = (raw_labels == 0).sum()
            print(f"📊 Rebalanced: {n_pos} positive, {n_neg} negative")

    labels_t = torch.tensor(raw_labels.values, dtype=torch.long)

    # 6) Train GNN (small epochs for demo)
    print("Training GraphComm-Lite (CPU) for sample", sample_name)
    model = train_graph_model(features_t, adj_t, labels_t, epochs=30)

    # 7) Get embeddings and RF demo
    with torch.no_grad():
        embeddings = model(features_t, adj_t).detach().numpy()

    rf = drug_response_prediction(embeddings, raw_labels.values)
    preds = rf.predict(embeddings)
    try:
        probs = rf.predict_proba(embeddings)[:, 1]
    except Exception:
        probs = None

    acc = (preds == raw_labels.values).mean()
    print(f"{sample_name} - training proxy accuracy: {acc:.4f}")

    # 8) Save results CSV
    out_csv = OUT_DIR / f"{sample_name}_results.csv"
    out_df = pd.DataFrame({'cell': df_expr.index.astype(str), 'label': raw_labels.values, 'predicted': preds})
    if probs is not None:
        out_df['pred_prob'] = probs
    out_df.to_csv(out_csv, index=False)
    print(f"Saved results to {out_csv}")

    # 9) Plots (prefix with sample)
    try:
        igf_vals = df_expr[[g for g in ["IGF1","IGF2","IGF1R"] if g in df_expr.columns]].mean(axis=1)
        plot_igf_activity(igf_vals, save_path=PLOTS_DIR / f"igf_activity_{sample_name}.png")
    except Exception:
        pass
    try:
        scores = embeddings[:, 0] if embeddings.ndim > 1 else embeddings
        plot_dormancy_heatmap(scores, save_path=PLOTS_DIR / f"dormancy_heatmap_{sample_name}.png")
    except Exception:
        pass
    try:
        plot_communication_graph(G_knn, scores, save_path=PLOTS_DIR / f"communication_graph_{sample_name}.png")
    except Exception:
        pass

    return {'name': sample_name, 'accuracy': acc, 'n_cells': df_expr.shape[0]}


def main():
    # convert any .h5 files first
    h5_files = list(DATA_DIR.glob('*.h5'))
    for h5 in h5_files:
        convert_h5_if_needed(h5)

    # find scRNA_*.csv files
    csvs = sorted(DATA_DIR.glob('scRNA_*.csv'))
    if not csvs:
        print('No scRNA_*.csv files found in data/. Please add per-sample CSVs or H5s and re-run.')
        sys.exit(1)

    summary = []
    for csv in csvs:
        sample_name = csv.stem.replace('scRNA_','')
        try:
            res = process_csv(csv, sample_name)
            summary.append(res)
        except Exception as e:
            print(f"Failed processing {csv}: {e}")

    print('\n=== Summary ===')
    for r in summary:
        print(f"{r['name']}: cells={r['n_cells']} accuracy={r['accuracy']:.4f}")


if __name__ == '__main__':
    main()
