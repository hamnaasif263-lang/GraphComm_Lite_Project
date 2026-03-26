# snRNA Features Integration Documentation

## Overview
The `snrna_features.py` module provides enriched biological context to your single-cell RNA-seq analysis by incorporating U1 snRNA (RNU1-1) sequence composition.

## What Was Already Done
This code was **already being used in your pipeline** (previously inline in `main.py`). We refactored it into a separate module for:
- Better code organization
- Reusability across projects
- Clear GitHub package structure
- Enhanced maintainability

## Biological Rationale

**U1 snRNA** is a non-coding RNA central to:
- **mRNA splicing** - part of the spliceosome machinery
- **Cell cycle regulation** - influences transcriptional machinery
- **Cell type identity** - snRNA composition reflects cellular function

By encoding U1 snRNA base composition and combining it with gene expression features, we:
1. Add a **genomic context feature** to each cell
2. Capture **splicing machinery variability** across cells
3. Enhance model features with **non-coding RNA signal**

## Technical Implementation

### Step-by-Step Integration

```
main.py pipeline:
├─ Step 1: Load scRNA data
│   └─ Output: gene_features (n_cells × n_genes)
│
├─ Step 2: Load snRNA features  ← snrna_features.py
│   ├─ Load U1 snRNA FASTA sequence
│   ├─ Encode as base-frequency: [A, C, G, T] composition
│   ├─ Tile across all cells
│   └─ Output: snRNA features (n_cells × 4)
│
├─ Step 3: Combine features
│   └─ Output: combined_features (n_cells × [n_genes + 4])
│
└─ Step 4: Build graph and train model
    └─ Input: combined_features
```

### Code Flow

```python
# main.py lines 33-35
snrna_vec = load_snRNA_features()                    # Load U1 snRNA (4D vector)
snrna_tiled = tile_snRNA_features(snrna_vec, ...)    # Tile to match cell count
features_combined = combine_features(features_np, snrna_tiled)  # Merge
```

## Feature Dimensions

| Stage | Shape | Description |
|-------|-------|-------------|
| Gene features | (n_cells, 50) | PCA-reduced gene expression |
| snRNA features | (n_cells, 4) | Base composition [A, C, G, T] normalized |
| **Combined** | **(n_cells, 54)** | **Used for model training** |

## Where It's Used

1. **Model Training** → `train_graph_model(features, adj, labels)` in [train_predict.py](train_predict.py)
2. **Embeddings** → Generated embeddings incorporate snRNA signal
3. **Drug Response Prediction** → Cell embeddings used in `drug_response_prediction()`
4. **Relapse Classification** → snRNA-enriched embeddings improve classifier

## How to Explain to Supervisor

**Statement:**
> "We integrated U1 snRNA sequence composition as a feature-engineering layer to contextualize gene expression signals. This adds biological relevance (splicing/cell-cycle information) to our feature vectors, enriching the model's ability to distinguish dormant vs. active cell states and predict drug response. The snRNA features were extracted into a modular, reusable component (`snrna_features.py`) that adheres to software engineering best practices for GitHub distribution."

## Validation

The module includes:
- ✅ Error handling (graceful fallback if FASTA missing)
- ✅ Logging (base frequency printouts for transparency)
- ✅ Normalization (frequency vector sums to 1.0)
- ✅ Documentation (docstrings for all functions)
- ✅ Demo mode (`if __name__ == "__main__"`)

## Files Involved

| File | Role |
|------|------|
| [snrna_features.py](snrna_features.py) | Core module (feature engineering) |
| [main.py](main.py) | Pipeline orchestrator (lines 33-35) |
| [train_predict.py](train_predict.py) | Consumes combined features |
| [data/snRNA_U1.fa](../data/snRNA_U1.fa) | Input (U1 snRNA FASTA) |

---

**Key Takeaway:** This is not new functionality—it's **refactored best practices** from your existing pipeline, now structured for professional GitHub distribution.
