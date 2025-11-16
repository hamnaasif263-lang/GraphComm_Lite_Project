# GraphComm-Lite

Lightweight prototype of GraphComm for investigating IGF and related pathways
in LUAD single-cell data, with relapse analysis and drug integration.

This repository contains code to:
- preprocess scRNA-seq from `data/scRNA_data.csv`
- build a simple graph-based model (`GraphCommLite`)
- evaluate relapse-associated signals and identify drug targets
- map targets to drugs via DGIdb or a local fallback CSV

See `main.py` for the pipeline entrypoint.

To run locally (Windows PowerShell):

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python -u main.py
```

If you want to push this project to GitHub, either install and login with `gh` (GitHub CLI) or create a new repository on GitHub and add it as a remote.
