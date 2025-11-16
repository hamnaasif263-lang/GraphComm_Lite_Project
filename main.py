"""
GraphComm-Lite  —  CPU-only prototype
AI-driven discovery of dormant-cell communication networks
and FDA-approved drug repurposing (IGF-pathway focus)
"""

import torch
import pandas as pd
import numpy as np
from Bio import SeqIO

from data_preprocess import preprocess_scRNA, extract_igf_expression
from graph_model import GraphCommLite
from train_predict import train_graph_model, drug_response_prediction
from pathway_analysis import build_igf_graph
from drug_integration import query_fda_approved_drugs
from relapse_analysis import get_relapse_labels, evaluate_classifier, top_upregulated_genes
from visualize import plot_igf_activity, plot_dormancy_heatmap, plot_communication_graph

# ============================================================
# 1️⃣  Load and preprocess single-cell RNA-seq data
# ============================================================

#  👇 Replace this path with your downloaded LUAD file
sc_path = "data/scRNA_data.csv"
features_np, A_knn, G_knn = preprocess_scRNA(sc_path)

# ============================================================
# 2️⃣  Load the U1 snRNA (RNU1-1) sequence and encode it
# ============================================================

def load_snRNA_features(fasta_path="data/snRNA_U1.fa"):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    seq = str(record.seq).upper()
    # simple base-frequency encoding
    vec = np.array([seq.count(b) for b in ["A", "C", "G", "T"]], dtype=float)
    vec = vec / vec.sum()
    return vec

snrna_vec = load_snRNA_features()
snrna_tiled = np.tile(snrna_vec, (features_np.shape[0], 1))
features_combined = np.hstack([features_np, snrna_tiled])

# ============================================================
# 3️⃣  Build IGF-pathway communication graph
# ============================================================

# read same scRNA data as DataFrame to check expression
try:
    df_expr = pd.read_csv(sc_path, index_col=0)
except Exception:
    df_expr = pd.read_csv(sc_path)

# build IGF graph using expression dataframe
A_igf = build_igf_graph(df_expr, threshold=0.25)

# ============================================================
# 4️⃣  Train GraphComm-Lite on CPU
# ============================================================

features = torch.tensor(features_combined, dtype=torch.float32)
adj = torch.tensor(A_igf, dtype=torch.float32)

# Derive labels from real data: use IGF-pathway activity per cell
# If an explicit label column exists in the CSV, prefer it. Otherwise use median-split on IGF pathway score.
if "label" in df_expr.columns:
    raw_labels = df_expr["label"].astype(int)
else:
    igf_series = extract_igf_expression(df_expr, ["IGF1", "IGF2", "IGF1R"])  # pathway genes
    thresh = igf_series.median()
    raw_labels = (igf_series > thresh).astype(int)

labels = torch.tensor(raw_labels.values, dtype=torch.long)

model = train_graph_model(features, adj, labels)

with torch.no_grad():
    embeddings = model(features, adj).detach().numpy()

# ============================================================
# 5️⃣  Drug-response prediction (demo)
# ============================================================

# For drug-response demo we don't have measured responses; use IGF-high vs IGF-low as a proxy
drug_labels = raw_labels.values
rf_model = drug_response_prediction(embeddings, drug_labels)

sample_pred = rf_model.predict(embeddings[0].reshape(1, -1))
print("\nPredicted FDA-drug response for sample cell:",
    "Effective" if sample_pred[0] == 1 else "Ineffective")

# ============================================================
# Relapse prediction evaluation
# ============================================================
try:
    # Try to get relapse labels using a relapse-associated signature (example)
    relapse_signature = ["AKT1", "PIK3CA", "PTEN", "MAPK1", "MAPK3"]
    relapse_labels = get_relapse_labels(df_expr, signature=relapse_signature, label_columns=["relapse", "relapse_status"])
    eval_res = evaluate_classifier(embeddings, relapse_labels.values)
    print(f"\nRelapse prediction eval: {eval_res}")
except Exception as e:
    print(f"\nCould not compute relapse labels/evaluation: {e}")

# ============================================================
# Identify top drug targets from differential expression between relapse-high vs relapse-low
# ============================================================
try:
    # choose candidate druggable genes to inspect
    candidate_targets = ["IGF1R", "PIK3CA", "AKT1", "EGFR", "MET"]
    top_targets = top_upregulated_genes(df_expr, relapse_labels.values, top_k=3, candidate_genes=candidate_targets)
    print(f"\nTop upregulated candidate targets in relapse-high cells: {top_targets}")
    for gene in top_targets:
        df_drugs = query_fda_approved_drugs(gene)
        if not df_drugs.empty:
            print(f"\nFDA-approved / mapped drugs for {gene}:")
            # Print available columns safely
            cols = [c for c in ["drug_name", "interaction_source", "interaction_claim_source", "approval_status"] if c in df_drugs.columns]
            if cols:
                print(df_drugs[cols].head())
            else:
                print(df_drugs.head())
        else:
            print(f"\nNo drugs found for {gene} via API/local DB.")
except Exception as e:
    print(f"\nCould not identify drug targets: {e}")

# --------------------------------------------
# VISUALIZATION SECTION
# --------------------------------------------

import os
os.makedirs("plots", exist_ok=True)

# Plot IGF pathway activity (saved to plots/)
try:
    df_expr = pd.read_csv(sc_path, index_col=0)
    igf_genes = open("data/igf_pathway_genes.txt").read().splitlines()
    igf_values = df_expr[[g for g in igf_genes if g in df_expr.columns]].mean(axis=1)
    plot_igf_activity(igf_values, save_path="plots/igf_activity.png")
    print("Saved IGF activity plot to plots/igf_activity.png")
except Exception as e:
    print("Could not generate IGF activity plot:", e)

# Plot heatmap of model scores
try:
    scores = embeddings[:, 0] if embeddings.ndim > 1 else embeddings
    plot_dormancy_heatmap(scores, save_path="plots/dormancy_heatmap.png")
    print("Saved dormancy heatmap to plots/dormancy_heatmap.png")
except Exception as e:
    print("Could not generate dormancy heatmap:", e)

# Plot network graph
try:
    scores = embeddings[:, 0] if embeddings.ndim > 1 else embeddings
    plot_communication_graph(G_knn, scores, save_path="plots/communication_graph.png")
    print("Saved communication graph to plots/communication_graph.png")
except Exception as e:
    print("Could not generate communication graph:", e)

# ============================================================
# 6️⃣  Retrieve FDA-approved drugs for IGF/PI3K/AKT targets
# ============================================================

targets = ["IGF1R", "PIK3CA", "AKT1"]
for gene in targets:
    try:
        df_drugs = query_fda_approved_drugs(gene)
        if not df_drugs.empty:
            print(f"\nFDA-approved drugs for {gene}:")
            print(df_drugs[["drug_name", "interaction_claim_source"]].head())
        else:
            print(f"\nNo approved drugs found for {gene} (check API connectivity).")
    except Exception as e:
        print(f"\nCould not query {gene}: {e}")

print("\n✅ GraphComm-Lite pipeline completed successfully.")
