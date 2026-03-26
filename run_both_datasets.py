#!/usr/bin/env python3
"""Run preprocessing, graph building and prediction for two datasets.

Assumes the CSV files were created by `download_and_convert.py`:
- data/scRNA_LUAD.csv
- data/scRNA_airway.csv

This script runs the same pipeline (preprocess -> build IGF graph -> train GNN -> RF demo)
for each dataset and saves simple per-cell results to `results/*.csv`.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import torch

from data_preprocess import preprocess_scRNA, extract_igf_expression
from pathway_analysis import build_igf_graph
from train_predict import train_graph_model, drug_response_prediction


def load_snRNA_features(fasta_path="data/snRNA_U1.fa"):
    # lightweight fallback if fasta missing
    import os
    if not os.path.exists(fasta_path):
        # return small zero-vector tiled
        vec = np.array([0.25, 0.25, 0.25, 0.25], dtype=float)
        return vec
    try:
        from Bio import SeqIO
        record = next(SeqIO.parse(fasta_path, "fasta"))
        seq = str(record.seq).upper()
        vec = np.array([seq.count(b) for b in ["A", "C", "G", "T"]], dtype=float)
        vec = vec / vec.sum()
        return vec
    except Exception:
        return np.array([0.25, 0.25, 0.25, 0.25], dtype=float)


def process_dataset(csv_path, result_csv_path, name="dataset"):
    print(f"\n=== Processing {name}: {csv_path} ===")
    csv_path = Path(csv_path)
    if not csv_path.exists():
        print(f"File not found: {csv_path}")
        return None

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
    snrna_vec = load_snRNA_features()
    try:
        snrna_tiled = np.tile(snrna_vec, (features_np.shape[0], 1))
        features_combined = np.hstack([features_np, snrna_tiled])
    except Exception:
        features_combined = features_np

    # 5) Prepare tensors
    features_t = torch.tensor(features_combined, dtype=torch.float32)
    adj_t = torch.tensor(A_igf, dtype=torch.float32)

    # 6) Derive labels (IGF median-split) or fallback
    if "label" in df_expr.columns:
        raw_labels = df_expr["label"].astype(int)
    else:
        igf_series = extract_igf_expression(df_expr, ["IGF1", "IGF2", "IGF1R"])  # pathway genes
        thresh = igf_series.median()
        raw_labels = (igf_series > thresh).astype(int)

    labels_t = torch.tensor(raw_labels.values, dtype=torch.long)

    # 7) Train GNN (small epochs for demo)
    print("Training GraphComm-Lite (CPU)...")
    model = train_graph_model(features_t, adj_t, labels_t, epochs=30)

    # 8) Get embeddings
    with torch.no_grad():
        embeddings = model(features_t, adj_t).detach().numpy()

    # 9) RF demo and compute training accuracy
    rf = drug_response_prediction(embeddings, raw_labels.values)
    preds = rf.predict(embeddings)
    try:
        probs = rf.predict_proba(embeddings)[:, 1]
    except Exception:
        probs = None

    acc = (preds == raw_labels.values).mean()
    print(f"{name} - training proxy accuracy: {acc:.4f}")

    # 10) Save results
    res_dir = Path("results")
    res_dir.mkdir(exist_ok=True)

    out_df = pd.DataFrame({
        'cell': df_expr.index.astype(str),
        'label': raw_labels.values,
        'predicted': preds,
    })
    if probs is not None:
        out_df['pred_prob'] = probs

    out_df.to_csv(result_csv_path, index=False)
    print(f"Saved results to {result_csv_path}")

    return {'name': name, 'accuracy': acc, 'n_cells': df_expr.shape[0]}


def main():
    pairs = [
        ("data/scRNA_LUAD.csv", "results/luad_results.csv", "LUAD"),
        ("data/scRNA_airway.csv", "results/airway_results.csv", "Airway"),
    ]

    summary = []
    for csv_path, result_csv, name in pairs:
        res = process_dataset(csv_path, result_csv, name=name)
        if res is not None:
            summary.append(res)

    print("\n=== Summary ===")
    for r in summary:
        print(f"{r['name']}: cells={r['n_cells']} accuracy={r['accuracy']:.4f}")


if __name__ == '__main__':
    main()
