#!/usr/bin/env python3
"""Download two single-person 10x HDF5 matrices and convert them to CSVs.

Files downloaded:
- data/GSE171524_lung_tumor_filtered_feature_bc_matrix.h5
- data/GSE176078_filtered_feature_bc_matrix.h5

Converted CSVs written:
- data/scRNA_LUAD.csv
- data/scRNA_airway.csv

This script attempts a memory-safe conversion: it first tries to convert the AnnData
to a dense DataFrame; if that fails due to memory, it writes the CSV in small cell
chunks to avoid allocating the whole dense matrix at once.
"""
import os
import sys
import requests
import shutil
import gzip
import csv
from pathlib import Path

try:
    import scanpy as sc
except Exception:
    print("scanpy is required. Install with: pip install scanpy anndata numpy scipy")
    raise

from scipy import sparse

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

FILES = [
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171524/suppl/GSE171524_lung_tumor_filtered_feature_bc_matrix.h5", DATA_DIR / "GSE171524_lung_tumor_filtered_feature_bc_matrix.h5", "data/scRNA_LUAD.csv"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_filtered_feature_bc_matrix.h5", DATA_DIR / "GSE176078_filtered_feature_bc_matrix.h5", "data/scRNA_airway.csv"),
]

def download_file(url, out_path, chunk_size=1024*1024):
    out_path = Path(out_path)
    if out_path.exists():
        print(f"Skipped download; file already exists: {out_path}")
        return out_path
    print(f"Downloading {url} -> {out_path}")
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(out_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
    print(f"Downloaded: {out_path} ({out_path.stat().st_size} bytes)")
    return out_path


def ann_to_csv_h5(h5_path, out_csv, chunk_cells=500):
    h5_path = Path(h5_path)
    out_csv = Path(out_csv)
    print(f"Reading H5 (10x) with scanpy: {h5_path}")
    ad = sc.read_10x_h5(str(h5_path))
    # ad: observations=cells, variables=genes
    n_cells = ad.n_obs
    n_genes = ad.n_vars
    print(f"AnnData loaded: {n_cells} cells x {n_genes} genes (sparse={sparse.issparse(ad.X)})")

    # Per-user request: subsample LUAD to exactly 500 cells, keep airway full
    try:
        original_n = int(ad.n_obs)
    except Exception:
        original_n = n_cells

    # determine dataset by target CSV name
    out_name = out_csv.name.lower()
    if 'luad' in out_name:
        print(f"Dataset detected: LUAD (original cells={original_n})")
        if original_n > 500:
            ad = ad[:500, :]
            print(f"Subsampled LUAD -> kept {int(ad.n_obs)} cells (requested 500)")
        else:
            print(f"LUAD has {original_n} cells < 500; keeping all ({original_n})")
    else:
        print(f"Dataset detected: keep ALL cells for {out_name} (original cells={original_n})")

    # Try direct dense conversion first
    try:
        print("Attempting dense conversion via ad.to_df()...")
        df = ad.to_df()
        print("Dense conversion succeeded; writing CSV (this may take a while)...")
        df.to_csv(out_csv)
        print(f"Wrote {out_csv} (shape {df.shape})")
        return out_csv
    except MemoryError as me:
        print("Dense conversion failed (MemoryError). Falling back to chunked streaming conversion.")
    except Exception as e:
        print("Dense conversion raised an exception; will attempt chunked write. Exception:", e)

    # Chunked write: write header (genes) then write each cell row in small batches
    var_names = list(ad.var_names)
    obs_names = list(ad.obs_names)

    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        header = ['cell'] + var_names
        writer.writerow(header)

        # iterate by chunks of cells
        for start in range(0, n_cells, chunk_cells):
            end = min(start + chunk_cells, n_cells)
            block = ad[start:end, :].X
            if sparse.issparse(block):
                block = block.toarray()
            # block: shape (rows, genes)
            for i_row, row in enumerate(block):
                cell_id = obs_names[start + i_row]
                # convert floats to strings to avoid locale issues
                writer.writerow([cell_id] + [str(x) for x in row.tolist()])
            print(f"Wrote cells {start}:{end} -> {out_csv}")

    print(f"Finished chunked write to {out_csv}")
    return out_csv


def main():
    for url, h5_path, out_csv in FILES:
        try:
            downloaded = download_file(url, h5_path)
        except Exception as e:
            print(f"Failed to download {url}: {e}")
            continue

        try:
            ann_to_csv_h5(downloaded, out_csv)
        except Exception as e:
            print(f"Failed to convert {downloaded} to CSV: {e}")
            continue

    print("All done. Generated CSVs (if downloads/conversions succeeded):")
    for _u, _h, out_csv in FILES:
        print(" -", out_csv)


if __name__ == '__main__':
    main()
