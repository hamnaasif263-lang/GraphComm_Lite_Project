#!/usr/bin/env python
"""PRAD_P01 - Pure numpy solution, no pandas"""

import numpy as np
from pathlib import Path
import sys

DATA = Path('data/prostate/raw/GSM4203181_data.raw.matrix.txt')
OUT = Path('graphcomm/results/prostate')
OUT.mkdir(parents=True, exist_ok=True)

# Phase 1: Parse file
print("Phase 1: Reading file...")
genes = []
cells_header = None
data = []

with open(DATA) as f:
    cells_header = f.readline().strip().split('\t')[1:]
    
    for i, line in enumerate(f):
        parts = line.strip().split('\t')
        genes.append(parts[0])
        data.append(np.array(parts[1:], dtype=np.float32))
        
        if (i + 1) % 5000 == 0:
            print(f"  {i+1} genes read")

genes = np.array(genes, dtype='<U50')
cells = np.array(cells_header, dtype='<U50')
counts = np.stack(data, axis=0)  # genes x cells

print(f"Matrix: {counts.shape}")

# Phase 2: Filter
print("Phase 2: Filtering...")
genes_per_cell = (counts > 0).sum(axis=0)
cell_mask = genes_per_cell >= 200
counts = counts[:, cell_mask]
cells = cells[cell_mask]

cells_per_gene = (counts > 0).sum(axis=1)
gene_mask = cells_per_gene >= 3
counts = counts[gene_mask, :]
genes = genes[gene_mask]

print(f"  Filtered: {counts.shape}")

# Phase 3: Normalize & Log
print("Phase 3: Normalizing...")
lib_size = counts.sum(axis=0)
counts_norm = counts * (1e4 / lib_size)[np.newaxis, :]
counts_log = np.log1p(counts_norm)

# Phase 4: Save
print("Phase 4: Saving...")
with open(OUT / 'scRNA_PRAD_P01_norm.csv', 'w') as f:
    f.write('\t' + '\t'.join(cells) + '\n')
    for i, g in enumerate(genes):
        f.write(g + '\t' + '\t'.join(counts_log[i, :].astype(str)) + '\n')
        if (i + 1) % 5000 == 0:
            print(f"    Saved {i+1} genes")

with open(OUT / 'scRNA_PRAD_P01_summary.txt', 'w') as f:
    f.write(f"Cells: {counts.shape[1]}\nGenes: {counts.shape[0]}\nLib Mean: {lib_size.mean():.0f}\n")

# Vis sample
vis_idx = np.random.choice(counts_log.shape[1], min(3000, counts_log.shape[1]), replace=False)
with open(OUT / 'scRNA_PRAD_P01_vis_sample.csv', 'w') as f:
    f.write('\t' + '\t'.join(cells[vis_idx]) + '\n')
    for i, g in enumerate(genes):
        f.write(g + '\t' + '\t'.join(counts_log[i, vis_idx].astype(str)) + '\n')

# Phase 5: IGF
print("Phase 5: IGF Analysis...")
IGF = {"IGF1", "IGF2", "IGF1R", "IRS1", "IRS2", "AKT1", "MTOR", "FOXO1"}
igf_idx = [j for j, g in enumerate(genes) if g in IGF]
print(f"  Found {len(igf_idx)} IGF genes")

if igf_idx:
    igf_scores = counts_log[igf_idx, :].mean(axis=0)
    
    with open(OUT / 'PRAD_P01_IGF_cell_scores.csv', 'w') as f:
        f.write('cell_barcode,igf_activity_score\n')
        for c, s in zip(cells, igf_scores):
            f.write(f'{c},{s}\n')
    
    with open(OUT / 'PRAD_P01_IGF_summary.csv', 'w') as f:
        f.write('metric,value\n')
        f.write(f'n_genes,{len(igf_idx)}\n')
        f.write(f'n_cells,{len(igf_scores)}\n')
        f.write(f'mean,{igf_scores.mean()}\n')
        f.write(f'median,{np.median(igf_scores)}\n')
        f.write(f'std,{igf_scores.std()}\n')
        f.write(f'min,{igf_scores.min()}\n')
        f.write(f'max,{igf_scores.max()}\n')
    
    # Plots
    print("Phase 6: Plots...")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
        ax.hist(igf_scores, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
        ax.axvline(igf_scores.mean(), color='red', linestyle='--', label=f'Mean: {igf_scores.mean():.3f}')
        ax.axvline(np.median(igf_scores), color='green', linestyle='--', label=f'Median: {np.median(igf_scores):.3f}')
        ax.set_xlabel('IGF Score')
        ax.set_ylabel('Cells')
        ax.set_title('IGF Pathway Activity PRAD_P01')
        ax.legend()
        ax.grid(alpha=0.3)
        Path('graphcomm/plots/prostate').mkdir(parents=True, exist_ok=True)
        plt.savefig('graphcomm/plots/prostate/PRAD_P01_IGF_activity_histogram.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("  Histogram done")
        
        top30 = np.argsort(igf_scores)[-30:][::-1]
        fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
        ax.barh(range(30), igf_scores[top30], color='steelblue')
        ax.set_yticklabels([str(c)[:15] for c in cells[top30]], fontsize=8)
        ax.set_xlabel('IGF Score')
        ax.set_title('Top 30 Cells')
        ax.invert_yaxis()
        plt.savefig('graphcomm/plots/prostate/PRAD_P01_IGF_top30_cells.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("  Top 30 done")
    except:
        print("  Plots skipped")

print("\n✅ COMPLETED")
print(f"Output: {OUT}")
