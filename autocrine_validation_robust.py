"""
Robust Python Script to Validate Autocrine Signaling in LUAD_P01
Checks for IGF2 (ligand) and IGF1R (receptor) expression
Creates visualization based on gene availability
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def validate_autocrine_signaling():
    """
    Main function to validate autocrine signaling in LUAD_P01 dataset
    """
    
    print("="*70)
    print("AUTOCRINE SIGNALING VALIDATION - LUAD_P01")
    print("="*70)
    
    # ==================== STEP 1: LOAD DATASET ====================
    print("\n[STEP 1] Loading dataset...")
    
    # Try to find the file in multiple locations
    possible_paths = [
        "scRNA_LUAD_P01_norm.csv",
        "data/scRNA_LUAD_P01_norm.csv",
        "results/scRNA_LUAD_P01_norm.csv",
        "data/luad/processed/scRNA_LUAD_P01_norm.csv"
    ]
    
    df = None
    loaded_path = None
    
    for path in possible_paths:
        if os.path.exists(path):
            try:
                print(f"  → Found file: {path}")
                df = pd.read_csv(path, index_col=0)
                loaded_path = path
                print(f"  ✓ Successfully loaded: {path}")
                print(f"  ✓ Dataset shape: {df.shape} (cells × genes)")
                break
            except Exception as e:
                print(f"  ✗ Error loading {path}: {e}")
                continue
    
    if df is None:
        print("\n✗ ERROR: Could not find or load scRNA_LUAD_P01_norm.csv")
        print("  Searched locations:")
        for path in possible_paths:
            print(f"    - {path}")
        return False
    
    # ==================== STEP 2: CHECK FOR GENES ====================
    print("\n[STEP 2] Checking for gene expression data...")
    
    genes_to_check = ['IGF2', 'IGF1R']
    found_genes = {}
    missing_genes = []
    
    print(f"  → Total genes in dataset: {df.shape[1]}")
    print(f"  → Searching for: {genes_to_check}")
    
    for gene in genes_to_check:
        if gene in df.columns:
            found_genes[gene] = df[gene]
            expr_count = (df[gene] > 0).sum()
            expr_percent = (expr_count / len(df)) * 100
            print(f"  ✓ {gene}: FOUND (expressed in {expr_count} cells, {expr_percent:.1f}%)")
        else:
            missing_genes.append(gene)
            print(f"  ✗ {gene}: NOT FOUND")
    
    # ==================== STEP 3: CREATE VISUALIZATION ====================
    print("\n[STEP 3] Creating visualization...")
    
    fig, ax = plt.subplots(figsize=(12, 8), facecolor='white')
    
    if len(found_genes) == 2:
        # CASE 1: BOTH GENES FOUND - Create scatter plot
        print("  ✓ Both genes detected - creating scatter plot...")
        
        igf2_expr = found_genes['IGF2'].values
        igf1r_expr = found_genes['IGF1R'].values
        
        # Create scatter plot
        scatter = ax.scatter(
            igf2_expr, 
            igf1r_expr,
            c=igf2_expr + igf1r_expr,  # Color by combined expression
            s=100,
            alpha=0.6,
            cmap='viridis',
            edgecolors='black',
            linewidth=0.5
        )
        
        # Add labels and title
        ax.set_xlabel('IGF2 Expression (Ligand)', fontsize=14, fontweight='bold')
        ax.set_ylabel('IGF1R Expression (Receptor)', fontsize=14, fontweight='bold')
        ax.set_title('Autocrine Loop Validation\nLUAD_P01 Analysis', 
                     fontsize=16, fontweight='bold', pad=20)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Combined Expression', fontsize=12, fontweight='bold')
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # Calculate and display correlation
        valid_mask = (igf2_expr > 0) & (igf1r_expr > 0)
        if valid_mask.sum() > 0:
            correlation = np.corrcoef(igf2_expr[valid_mask], igf1r_expr[valid_mask])[0, 1]
            n_autocrine = valid_mask.sum()
            text_info = f"Cells with both genes: {n_autocrine}\nCorrelation (IGF2-IGF1R): {correlation:.3f}"
            ax.text(0.98, 0.02, text_info, 
                   transform=ax.transAxes,
                   fontsize=11,
                   verticalalignment='bottom',
                   horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                   family='monospace')
        
        result_status = "✓ POSITIVE"
        
    else:
        # CASE 2: GENES NOT FOUND - Create error message plot
        print("  ✗ Genes not detected - creating negative result visualization...")
        
        ax.set_facecolor('white')
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis('off')
        
        # Add red text message
        ax.text(
            0.5, 0.6,
            'NEGATIVE RESULT',
            ha='center', va='center',
            fontsize=32, fontweight='bold',
            color='red',
            transform=ax.transAxes,
            family='monospace'
        )
        
        ax.text(
            0.5, 0.45,
            'IGF2 Ligand Not Detected',
            ha='center', va='center',
            fontsize=28, fontweight='bold',
            color='darkred',
            transform=ax.transAxes,
            family='monospace'
        )
        
        # Add detailed info
        missing_text = f"Missing genes: {', '.join(missing_genes)}\n\n"
        missing_text += f"Technical Dropout Detected:\n"
        missing_text += f"• Expected genes not found in dataset\n"
        missing_text += f"• Dataset shape: {df.shape}\n"
        missing_text += f"• This indicates genes were filtered\n"
        missing_text += f"  during preprocessing"
        
        ax.text(
            0.5, 0.25,
            missing_text,
            ha='center', va='center',
            fontsize=12,
            color='black',
            transform=ax.transAxes,
            family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8, pad=1)
        )
        
        result_status = "✗ NEGATIVE"
    
    # ==================== STEP 4: SAVE PLOT ====================
    print("\n[STEP 4] Saving visualization...")
    
    output_file = "LUAD_Autocrine_Check.png"
    
    try:
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  ✓ Plot saved: {output_file}")
        file_size = os.path.getsize(output_file) / 1024
        print(f"  ✓ File size: {file_size:.1f} KB")
    except Exception as e:
        print(f"  ✗ Error saving plot: {e}")
        return False
    
    plt.close()
    
    # ==================== SUMMARY ====================
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)
    print(f"\nDataset loaded from: {loaded_path}")
    print(f"Dataset dimensions: {df.shape[0]} cells × {df.shape[1]} genes")
    print(f"\nGenes searched for: {genes_to_check}")
    print(f"Genes found: {list(found_genes.keys()) if found_genes else 'NONE'}")
    print(f"Genes missing: {missing_genes if missing_genes else 'NONE'}")
    print(f"\nResult: {result_status}")
    print(f"Output file: {output_file}")
    print("\n" + "="*70)
    
    return True


if __name__ == "__main__":
    try:
        success = validate_autocrine_signaling()
        if success:
            print("\n✓ Validation completed successfully!")
            sys.exit(0)
        else:
            print("\n✗ Validation failed!")
            sys.exit(1)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
