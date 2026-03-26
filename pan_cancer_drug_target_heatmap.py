"""
Pan-Cancer Therapeutic Target Heatmap
Validates drug repurposing candidates across BC_P01, LUAD_P01, and PRAD_P01

This script:
1. Loads expression data (handling multiple formats)
2. Filters for Tumor Epithelial Cells
3. Calculates mean expression for druggable gene panels
4. Z-scores across cancer types
5. Generates heatmap visualization
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

OUTPUT_DIR = Path("results")
OUTPUT_DIR.mkdir(exist_ok=True)
PLOTS_DIR = Path("plots")
PLOTS_DIR.mkdir(exist_ok=True)

# Druggable gene panels (organized by therapeutic pathway)
DRUGGABLE_PANELS = {
    "IGF/Insulin Axis\n(Target: Linsitinib)": ["IGF1R", "INSR", "IRS1", "IGF1", "IGF2"],
    "PI3K/mTOR Axis\n(Target: Everolimus/Alpelisib)": ["PIK3CA", "AKT1", "MTOR", "RPS6KB1", "PI3K"],
    "Cell Cycle/CDK\n(Target: Palbociclib/Ribociclib)": ["CDK4", "CDK6", "CCND1", "CDK2"],
    "Dormancy/Suppression\n(Target: Metformin)": ["CDKN1B", "FOXO1", "FOXO3"],
}

# Flatten gene list for easy access
ALL_GENES = []
GENE_TO_PATHWAY = {}
for pathway, genes in DRUGGABLE_PANELS.items():
    ALL_GENES.extend(genes)
    for gene in genes:
        GENE_TO_PATHWAY[gene] = pathway

print(f"\n📊 PAN-CANCER THERAPEUTIC TARGET HEATMAP")
print(f"{'='*70}")
print(f"Analyzing {len(ALL_GENES)} druggable genes across 3 cancer types")

# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

def load_long_format(filepath):
    """Load long-format expression data (cell, gene, count columns)"""
    df = pd.read_csv(filepath)
    if 'cell' not in df.columns or 'gene' not in df.columns or 'count' not in df.columns:
        return None
    
    # Pivot to matrix format
    matrix = df.pivot_table(index='cell', columns='gene', values='count', fill_value=0)
    return matrix

def load_matrix_format(filepath):
    """Load matrix-format expression data (cells × genes)"""
    try:
        df = pd.read_csv(filepath, index_col=0)
        return df
    except:
        return None

def load_parquet_format(filepath):
    """Load parquet-format expression data"""
    try:
        df = pd.read_parquet(filepath)
        return df
    except:
        return None

def load_patient_data_flexible(patient_id, possible_paths):
    """Try to load data from multiple possible paths"""
    print(f"\n📥 Loading {patient_id}...", end=" ")
    
    for path in possible_paths:
        path_obj = Path(path)
        if not path_obj.exists():
            continue
        
        # Try different loading methods
        if path.endswith('.parquet'):
            expr_df = load_parquet_format(path)
        elif path.endswith('.csv'):
            # Try long format first
            expr_df = load_long_format(path)
            if expr_df is None:
                # Try matrix format
                expr_df = load_matrix_format(path)
        else:
            continue
        
        if expr_df is not None and expr_df.shape[0] > 0 and expr_df.shape[1] > 5:
            print(f"✓ ({expr_df.shape[0]:,} cells × {expr_df.shape[1]} genes)")
            return expr_df
    
    print(f"❌ NOT FOUND")
    return None

# ============================================================================
# CELL TYPE FILTERING
# ============================================================================

def filter_tumor_epithelial_cells(expr_df, patient_id):
    """Filter for tumor epithelial cells using available metadata"""
    
    # Try to find cell type annotations
    celltype_paths = [
        f"results/{patient_id}_celltype_annotations.csv",
        f"data/*/results/{patient_id}_celltype_annotations.csv",
        f"graphcomm/results/*/{patient_id}_celltype_annotations.csv",
    ]
    
    celltype_df = None
    for path in celltype_paths:
        if Path(path).exists():
            celltype_df = pd.read_csv(path)
            break
    
    if celltype_df is None:
        print(f"  ⚠️  No cell type annotations found; using all cells")
        return expr_df
    
    # Find cell ID column and cell type column
    # Different files have different column names
    cell_col = None
    celltype_col = None
    
    for col in celltype_df.columns:
        if 'cell' in col.lower() and 'type' not in col.lower():
            cell_col = col
        elif 'type' in col.lower() or 'classification' in col.lower():
            celltype_col = col
    
    if cell_col is None or celltype_col is None:
        print(f"  ⚠️  Could not identify cell type columns; using all cells")
        return expr_df
    
    # Filter for epithelial-like cell types
    epithelial_keywords = ["epithelial", "epithelium", "tumor", "cancer", "carcinoma", "ep", "t_epi"]
    mask = celltype_df[celltype_col].astype(str).str.lower().str.contains(
        '|'.join(epithelial_keywords), 
        na=False
    )
    
    if mask.sum() == 0:
        print(f"  ⚠️  No epithelial cells found using keywords; using top 50% by cell count")
        # Use top 50% of cells by cell ID length (heuristic)
        return expr_df.head(int(expr_df.shape[0] * 0.5))
    
    epithelial_cells = celltype_df.loc[mask, cell_col].tolist()
    # Filter expression data
    filtered_expr = expr_df.loc[expr_df.index.intersection(epithelial_cells)]
    
    print(f"  ✓ Filtered to {filtered_expr.shape[0]:,} Tumor Epithelial Cells")
    return filtered_expr

# ============================================================================
# EXPRESSION ANALYSIS
# ============================================================================

def calculate_mean_expression(expr_df, genes):
    """Calculate mean expression for genes in the panel"""
    
    # Filter to genes that exist
    available_genes = [g for g in genes if g in expr_df.columns]
    missing_genes = [g for g in genes if g not in expr_df.columns]
    
    if missing_genes:
        print(f"  ⚠️  Missing: {', '.join(missing_genes[:3])}", end="")
        if len(missing_genes) > 3:
            print(f" + {len(missing_genes)-3} more")
        else:
            print()
    
    if not available_genes:
        print(f"  ❌ No genes found!")
        return {}
    
    mean_expr = {}
    for gene in available_genes:
        # Use mean of non-zero values (common in single-cell data)
        values = expr_df[gene]
        values_nonzero = values[values > 0]
        if len(values_nonzero) > 0:
            mean_expr[gene] = values_nonzero.mean()
        else:
            mean_expr[gene] = values.mean()
    
    return mean_expr

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def main():
    """Generate Pan-Cancer Therapeutic Target Heatmap"""
    
    # 1. Load data for all patients
    patient_data_paths = {
        "BC_P01": [
            "data/breast_cancer/raw/scRNA_BC_P01_norm.csv",
            "data/breast_cancer/processed/scRNA_BC_P01_norm.csv",
        ],
        "LUAD_P01": [
            "data/luad/processed/scRNA_LUAD_P01_norm.parquet",
            "data/luad/processed/scRNA_LUAD_P01_norm.csv",
        ],
        "PRAD_P01": [
            "graphcomm/results/prostate/scRNA_PRAD_P01_norm.csv",
            "data/prostate/processed/scRNA_PRAD_P01_norm.csv",
        ],
    }
    
    all_data = {}
    for patient_id, paths in patient_data_paths.items():
        expr_df = load_patient_data_flexible(patient_id, paths)
        
        if expr_df is None:
            print(f"❌ Skipping {patient_id}")
            continue
        
        # Handle transposed data (genes as rows)
        if expr_df.shape[1] > expr_df.shape[0]:
            # Likely cells > genes, correct orientation
            pass
        elif expr_df.shape[0] < 100:
            # Likely transposed (genes as rows, cells as columns)
            print(f"  ℹ️  Transposing data (detected genes as rows)")
            expr_df = expr_df.T
        
        # Filter for tumor epithelial cells
        expr_filtered = filter_tumor_epithelial_cells(expr_df, patient_id)
        
        # Calculate mean expression
        mean_expr = calculate_mean_expression(expr_filtered, ALL_GENES)
        all_data[patient_id] = mean_expr
    
    if not all_data:
        print("❌ No patient data loaded!")
        return
    
    print(f"\n✓ Successfully loaded {len(all_data)} patients")
    
    # 2. Create expression matrix
    print(f"\n📊 Creating expression matrix...")
    
    # Ensure all patients have all genes (fill missing with 0)
    for patient_id in all_data:
        for gene in ALL_GENES:
            if gene not in all_data[patient_id]:
                all_data[patient_id][gene] = 0.0
    
    # Convert to DataFrame
    heatmap_df = pd.DataFrame(all_data).T
    heatmap_df = heatmap_df[ALL_GENES]
    
    print(f"  ✓ Matrix shape: {heatmap_df.shape}")
    print(f"\n  Raw Expression Matrix:")
    print(heatmap_df.round(3))
    
    # 3. Z-score normalization (crucial!)
    print(f"\n📈 Z-score normalization (across cancer types)...")
    
    heatmap_scaled = heatmap_df.copy()
    for gene in heatmap_scaled.columns:
        mean_val = heatmap_scaled[gene].mean()
        std_val = heatmap_scaled[gene].std()
        
        if std_val > 0:
            heatmap_scaled[gene] = (heatmap_scaled[gene] - mean_val) / std_val
        else:
            heatmap_scaled[gene] = 0
    
    print(f"  ✓ Scaled values (z-score):")
    print(heatmap_scaled.round(2))
    
    # 4. Create heatmap
    print(f"\n🎨 Creating visualization...")
    
    fig, ax = plt.subplots(figsize=(15, 6))
    
    # Blue to Red diverging colormap
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    
    # Symmetric limits
    vmax = np.abs(heatmap_scaled.values).max()
    vmin = -vmax
    
    # Plot
    sns.heatmap(
        heatmap_scaled,
        annot=True,
        fmt=".2f",
        cmap=cmap,
        center=0,
        vmin=vmin,
        vmax=vmax,
        cbar_kws={"label": "Z-Score (Standardized Expression)"},
        linewidths=0.5,
        linecolor="gray",
        ax=ax,
        square=False,
    )
    
    # Formatting
    ax.set_title(
        "Pan-Cancer Therapeutic Target Heatmap\nDrug Repurposing Candidate Validation",
        fontsize=14,
        fontweight="bold",
        pad=20
    )
    ax.set_xlabel("Druggable Genes (by Therapeutic Pathway)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Cancer Types", fontsize=12, fontweight="bold")
    
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    
    # Pathway separators
    pathway_boundaries = [0]
    pos = 0
    for pathway, genes in DRUGGABLE_PANELS.items():
        pos += len(genes)
        pathway_boundaries.append(pos)
    
    for boundary in pathway_boundaries[1:-1]:
        ax.axvline(x=boundary, color="black", linewidth=2, linestyle="--")
    
    plt.tight_layout()
    
    # 5. Save figure
    output_path = PLOTS_DIR / "Pan_Cancer_Drug_Target_Heatmap.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"  ✓ Saved: {output_path}")
    plt.close()
    
    # 6. Top therapeutic targets
    print(f"\n🎯 THERAPEUTIC TARGET SUMMARY")
    print(f"{'='*70}")
    
    top_targets = {}
    for patient_id in heatmap_df.index:
        top_gene = heatmap_scaled.loc[patient_id].idxmax()
        z_score = heatmap_scaled.loc[patient_id, top_gene]
        raw_expr = heatmap_df.loc[patient_id, top_gene]
        
        top_targets[patient_id] = {
            "gene": top_gene,
            "z_score": z_score,
            "raw_expr": raw_expr,
        }
        
        pathway = GENE_TO_PATHWAY.get(top_gene, "Unknown").split('\n')[0]
        print(f"\n{patient_id}:")
        print(f"  ⭐ Top Target: {top_gene}")
        print(f"  🔬 Pathway: {pathway}")
        print(f"  📊 Mean Expression: {raw_expr:.4f}")
        print(f"  📈 Vulnerability Score (Z-Score): {z_score:+.3f}")
    
    # 7. Gene-level summary
    print(f"\n📊 GENE-LEVEL RANKING (Mean Expression)")
    print(f"{'='*70}")
    
    gene_summary = heatmap_df.mean().sort_values(ascending=False)
    print("\nTop Druggable Genes (by Expression Level):")
    for i, (gene, value) in enumerate(gene_summary.items()[:10], 1):
        pathway = GENE_TO_PATHWAY.get(gene, "Unknown").split('\n')[0]
        print(f"  {i:2d}. {gene:10s} | Mean: {value:7.4f} | {pathway}")
    
    # 8. Drug repurposing insights
    print(f"\n💊 DRUG REPURPOSING INSIGHTS")
    print(f"{'='*70}")
    
    for patient_id in ["BC_P01", "LUAD_P01", "PRAD_P01"]:
        if patient_id not in heatmap_scaled.index:
            continue
        
        print(f"\n{patient_id} (Priority Therapeutic Pathways):")
        
        patient_data = heatmap_scaled.loc[patient_id]
        
        # Calculate pathway scores
        pathway_scores = {}
        for pathway, genes in DRUGGABLE_PANELS.items():
            pathway_genes = [g for g in genes if g in patient_data.index]
            if pathway_genes:
                avg_score = patient_data[pathway_genes].mean()
                pathway_scores[pathway] = avg_score
        
        # Sort by vulnerability
        sorted_pathways = sorted(pathway_scores.items(), key=lambda x: x[1], reverse=True)
        
        for i, (pathway, score) in enumerate(sorted_pathways[:3], 1):
            drugs_dict = {
                "IGF": "Linsitinib",
                "PI3K": "Everolimus, Alpelisib",
                "Cell Cycle": "Palbociclib, Ribociclib",
                "Dormancy": "Metformin",
            }
            
            drug_name = ""
            for key, drug in drugs_dict.items():
                if key in pathway:
                    drug_name = f"({drug})"
                    break
            
            pathway_name = pathway.split('\n')[0]
            print(f"  {i}. {pathway_name:28s} {drug_name:35s} Score: {score:+.3f}")
    
    # 9. Save results
    results_table = pd.DataFrame({
        "Patient": list(top_targets.keys()),
        "Top_Gene": [v["gene"] for v in top_targets.values()],
        "Raw_Expression": [v["raw_expr"] for v in top_targets.values()],
        "Z_Score": [v["z_score"] for v in top_targets.values()],
    })
    
    results_path = OUTPUT_DIR / "Pan_Cancer_Drug_Targets.csv"
    results_table.to_csv(results_path, index=False)
    print(f"\n✓ Saved results: {results_path}")
    
    heatmap_path = OUTPUT_DIR / "Pan_Cancer_Scaled_Expression.csv"
    heatmap_scaled.to_csv(heatmap_path)
    print(f"✓ Saved scaled expression matrix: {heatmap_path}")
    
    # Final summary
    print(f"\n{'='*70}")
    print(f"✅ ANALYSIS COMPLETE")
    print(f"{'='*70}")
    print(f"\n📌 TOP THERAPEUTIC TARGETS:")
    for patient_id, target in top_targets.items():
        print(f"   {patient_id}: {target['gene']}")


if __name__ == "__main__":
    main()
