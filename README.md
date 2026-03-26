# GraphComm-Lite Cancer Analysis

AI-driven discovery of dormant-cell communication networks and FDA-approved drug repurposing for cancer (IGF-pathway focus).

## Overview

GraphComm-Lite is a lightweight, CPU-friendly bioinformatics pipeline for:
- **Single-cell RNA-seq preprocessing** – normalize, reduce dimensionality, build cell-cell communication graphs
- **Graph-based modeling** – neural network for learning cell state embeddings
- **Dormancy prediction** – identify dormant vs. active cell states
- **Drug response prediction** – predict therapeutic response using cell embeddings
- **Relapse analysis** – identify relapse-associated gene signatures
- **Drug targets** – map upregulated genes to FDA-approved drugs via DGIdb API

## Key Features

✅ **Modular architecture** – Each processing step is a separate importable module  
✅ **No GPU required** – Full CPU-only implementation (PyTorch CPU backend)  
✅ **Biological feature engineering** – Integration of U1 snRNA composition  
✅ **Pathway analysis** – Compute scores for 12+ biological pathways  
✅ **API integration** – Query FDA-approved drugs for gene targets  
✅ **Professional visualizations** – IGF activity heatmaps, network graphs, embeddings  

## Pipeline Architecture

```
Step 1: Data Preprocessing (data_preprocess.py)
   └─ Load scRNA-seq → normalize → PCA → k-NN graph

Step 2: snRNA Features (snrna_features.py)
   └─ Load U1 snRNA sequence → encode base composition → tile

Step 3: IGF Pathway Graph (pathway_analysis.py)
   └─ Build cell-cell communication network from IGF expression

Step 4: Model Training (graph_model.py + train_predict.py)
   └─ Train GraphCommLite neural network on combined features

Step 5: Analysis & Visualization
   ├─ Drug response prediction (train_predict.py)
   ├─ Relapse prediction (relapse_analysis.py)
   ├─ FDA drug mapping (drug_integration.py)
   └─ Plots (visualize.py)
```

## Installation

### Prerequisites
- Python 3.8+
- Virtual environment (recommended)

### Setup (Windows PowerShell)

```powershell
# Clone the repository
git clone https://github.com/YOUR_USERNAME/Graphcomm-Lite-Cancer-Analysis.git
cd Graphcomm-Lite-Cancer-Analysis

# Create virtual environment
python -m venv .venv
.venv\Scripts\Activate.ps1

# Install dependencies
pip install -r requirements.txt
```

### Setup (macOS/Linux Bash)

```bash
git clone https://github.com/YOUR_USERNAME/Graphcomm-Lite-Cancer-Analysis.git
cd Graphcomm-Lite-Cancer-Analysis

python -m venv .venv
source .venv/bin/activate

pip install -r requirements.txt
```

## Usage

### Run Full Pipeline

```bash
python main.py
```

This executes all 6 steps and generates:
- **plots/** – Visualizations (IGF activity, dormancy scores, communication network)
- **results/** – Pathway scores, predictions
- Console output with predictions and drug targets

### Override Data Path

```bash
$env:SC_PATH="path/to/your/scRNA_data.csv"
python main.py
```

### Use Individual Modules

```python
from data_preprocess import preprocess_scRNA, extract_igf_expression
from snrna_features import load_snRNA_features, combine_features
from train_predict import train_graph_model, drug_response_prediction
from drug_integration import query_fda_approved_drugs

# Load and preprocess data
features, adj, graph = preprocess_scRNA("data/scRNA_data.csv")

# Add snRNA context
snrna_vec = load_snRNA_features()
combined_features = combine_features(features, snrna_vec)

# Train and predict
model = train_graph_model(features, adj, labels)
drugs = query_fda_approved_drugs("IGF1R")
```

## Module Reference

| Module | Purpose |
|--------|---------|
| [data_preprocess.py](data_preprocess.py) | scRNA-seq loading, scaling, PCA, k-NN graph |
| [snrna_features.py](snrna_features.py) | U1 snRNA feature engineering |
| [pathway_analysis.py](pathway_analysis.py) | IGF pathway graph, pathway scoring |
| [graph_model.py](graph_model.py) | GraphCommLite neural network architecture |
| [train_predict.py](train_predict.py) | Model training, drug response prediction |
| [drug_integration.py](drug_integration.py) | FDA drug queries via DGIdb API |
| [relapse_analysis.py](relapse_analysis.py) | Relapse prediction, differential analysis |
| [visualize.py](visualize.py) | Plotting utilities |
| [main.py](main.py) | Pipeline orchestrator |

## Data Format

### Input: scRNA-seq CSV
```
                  GENE1   GENE2   GENE3   ...
CELL_001          12.3    5.6     8.2
CELL_002          10.1    6.8     7.1
...
```

(Rows = cells, Columns = genes, values = log-normalized expression or counts)

### Output: Embeddings
- Shape: (n_cells, 2) – 2D model outputs for dormancy scoring
- Interpretable scores correlate with cell state predictions

## Biological Context

### IGF Pathway
The **Insulin-like Growth Factor (IGF) pathway** is critical for:
- Cell proliferation & survival
- Tumor progression in lung adenocarcinoma (LUAD)
- Dormancy escape mechanisms

### U1 snRNA Integration
- Captures **splicing machinery activity** across cells
- Non-coding RNA signals correlate with cell cycle state
- Enhances model features with genomic context

### Relapse Signatures
We identify genes (e.g., AKT1, PIK3CA, EGFR) whose upregulation correlates with relapse risk, then map them to FDA-approved drugs.

## Documentation

- [SNRNA_INTEGRATION.md](SNRNA_INTEGRATION.md) – Technical details on snRNA feature engineering
- [main.py](main.py) – Inline comments explaining each pipeline step

## Requirements

- PyTorch 2.0+
- scikit-learn 1.2+
- pandas 1.5+
- numpy 1.23+
- BioPython 1.81+
- networkx 3.0+
- matplotlib 3.5+ (visualization)

See [requirements.txt](requirements.txt) for all dependencies.

## Performance

Tested on:
- **Dataset:** LUAD scRNA-seq (20k cells × 5k genes)
- **Hardware:** CPU-only (Intel i7, 16 GB RAM)
- **Runtime:** ~2-3 minutes for full pipeline

## Citation

If you use this work, please cite:
```
GraphComm-Lite Cancer Analysis Pipeline
[Your Lab/University]
[Publication Year]
GitHub: https://github.com/YOUR_USERNAME/Graphcomm-Lite-Cancer-Analysis
```

## License

[Specify your license: MIT, GPL, proprietary, etc.]

## Contact & Support

For questions or issues:
- **Email:** [Your email]
- **Institution:** [Your university]
- **GitHub Issues:** https://github.com/YOUR_USERNAME/Graphcomm-Lite-Cancer-Analysis/issues

## Troubleshooting

### "FileNotFoundError: data/scRNA_data.csv"
- Place your scRNA-seq CSV in the `data/` folder
- Or set `SC_PATH` environment variable

### "No module named 'torch'"
- Run: `pip install -r requirements.txt`

### "DGIdb API timeout"
- Drug queries fall back to local CSV if API unavailable
- Place `data/drug_db.csv` for offline mode

## Version History

- **v1.0** (2026-03-26) – Initial release
  - Core pipeline with 9 modules
  - snRNA feature integration
  - Drug response & relapse prediction

---

**Made with ❤️ for cancer research**
