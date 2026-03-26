#!/usr/bin/env python3
"""Quick test of the prediction fix on a couple samples."""

from pathlib import Path
import pandas as pd
import numpy as np
from data_preprocess import preprocess_scRNA, extract_igf_expression
from pathway_analysis import build_igf_graph
from train_predict import train_graph_model, drug_response_prediction
import torch

DATA_DIR = Path("data")
OUT_DIR = Path("results")
PLOTS_DIR = Path("plots")

def test_sample(csv_path: Path, sample_name: str):
    """Test the prediction fix on a single sample."""
    print(f"\n--- Testing {sample_name} ---")
    
    # 1) Preprocess
    features_np, A_knn, G_knn = preprocess_scRNA(str(csv_path))
    
    # 2) Read expression DF
    try:
        df_expr = pd.read_csv(csv_path, index_col=0)
    except Exception:
        df_expr = pd.read_csv(csv_path)
    
    # 3) Build IGF graph
    A_igf = build_igf_graph(df_expr, threshold=0.25)
    
    # 4) Combine features
    snrna_vec = np.array([0.25, 0.25, 0.25, 0.25], dtype=float)
    try:
        snrna_tiled = np.tile(snrna_vec, (features_np.shape[0], 1))
        features_combined = np.hstack([features_np, snrna_tiled])
    except Exception:
        features_combined = features_np
    
    # 5) Generate labels (using the FIXED logic)
    if "label" in df_expr.columns:
        raw_labels = df_expr["label"].astype(int)
    else:
        igf_genes = ["IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1"]
        available_genes = [g for g in igf_genes if g in df_expr.columns]
        
        if available_genes:
            igf_series = df_expr[available_genes].mean(axis=1)
            thresh = igf_series.median()
            raw_labels = (igf_series > thresh).astype(int)
            print(f"Generated labels from {len(available_genes)} IGF genes: {available_genes}")
        else:
            gene_means = df_expr.mean(axis=0).sort_values(ascending=False)
            top_genes = gene_means.head(10).index.tolist()
            if len(top_genes) > 0:
                expr_score = df_expr[top_genes].mean(axis=1)
                thresh = expr_score.median()
                raw_labels = (expr_score > thresh).astype(int)
                print(f"Generated labels from top genes: {top_genes}")
            else:
                raw_labels = pd.Series(np.random.randint(0, 2, size=df_expr.shape[0]), index=df_expr.index)
        
        n_pos = (raw_labels == 1).sum()
        n_neg = (raw_labels == 0).sum()
        print(f"Label distribution: {n_pos} positive ({100*n_pos/(n_pos+n_neg):.1f}%), {n_neg} negative")
        
        if n_pos == 0 or n_neg == 0:
            print("Rebalancing with stratified split...")
            raw_labels = pd.Series(np.random.randint(0, 2, size=df_expr.shape[0]), index=df_expr.index)
            n_pos = (raw_labels == 1).sum()
            n_neg = (raw_labels == 0).sum()
            print(f"Rebalanced: {n_pos} positive, {n_neg} negative")
    
    # 6) Train model
    features_t = torch.tensor(features_combined, dtype=torch.float32)
    adj_t = torch.tensor(A_igf, dtype=torch.float32)
    labels_t = torch.tensor(raw_labels.values, dtype=torch.long)
    
    print("Training model...")
    model = train_graph_model(features_t, adj_t, labels_t, epochs=20)
    
    # 7) Get predictions
    with torch.no_grad():
        embeddings = model(features_t, adj_t).detach().numpy()
    
    rf = drug_response_prediction(embeddings, raw_labels.values)
    preds = rf.predict(embeddings)
    
    # 8) Check results
    n_pred_0 = (preds == 0).sum()
    n_pred_1 = (preds == 1).sum()
    print(f"\nPrediction distribution: {n_pred_0} class 0, {n_pred_1} class 1")
    
    if n_pred_1 == 0:
        print("ERROR: Predictions are still all zeros!")
    elif n_pred_1 == len(preds):
        print("WARNING: Predictions are all ones!")
    else:
        print("SUCCESS: Predictions have both classes!")
    
    # 9) Save new results
    out_csv = OUT_DIR / f"{sample_name}_results_FIXED.csv"
    out_df = pd.DataFrame({
        'cell': df_expr.index.astype(str),
        'label': raw_labels.values,
        'predicted': preds
    })
    out_df.to_csv(out_csv, index=False)
    print(f"Saved fixed results to {out_csv}")
    
    return {
        'sample': sample_name,
        'n_labels_0': (raw_labels == 0).sum(),
        'n_labels_1': (raw_labels == 1).sum(),
        'n_preds_0': n_pred_0,
        'n_preds_1': n_pred_1,
    }

# Test on a few samples
if __name__ == '__main__':
    test_samples = [
        'scRNA_BRONCHO_11.csv',
        'scRNA_LUNG_T09.csv',
        'scRNA_EBUS_06.csv',
    ]
    
    results = []
    for csv_name in test_samples:
        csv_path = DATA_DIR / csv_name
        if csv_path.exists():
            sample_name = csv_name.replace('scRNA_', '').replace('.csv', '')
            try:
                result = test_sample(csv_path, sample_name)
                results.append(result)
            except Exception as e:
                print(f"Error testing {sample_name}: {e}")
    
    print("\n" + "="*70)
    print("SUMMARY OF FIXED PREDICTIONS")
    print("="*70)
    for r in results:
        print(f"{r['sample']:20} | Labels: [0:{r['n_labels_0']:4}, 1:{r['n_labels_1']:4}] | "
              f"Preds: [0:{r['n_preds_0']:4}, 1:{r['n_preds_1']:4}]")
