# ==========================================================
# snrna_features.py
# ----------------------------------------------------------
# Handles loading and encoding of snRNA features (U1 snRNA)
# for GraphComm-Lite:
#   1. Loads FASTA sequence
#   2. Encodes as base-frequency vector
#   3. Tiles across all cells
# ==========================================================

import numpy as np
from Bio import SeqIO


def load_snRNA_features(fasta_path: str = "data/snRNA_U1.fa") -> np.ndarray:
    """
    Load and encode U1 snRNA (RNU1-1) sequence.
    Returns a simple base-frequency encoding (4D vector).
    
    Args:
        fasta_path: Path to FASTA file containing snRNA sequence
        
    Returns:
        np.ndarray: 4D vector with counts of [A, C, G, T] normalized to sum=1
    """
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        seq = str(record.seq).upper()
        # Simple base-frequency encoding
        vec = np.array([seq.count(b) for b in ["A", "C", "G", "T"]], dtype=float)
        vec = vec / vec.sum()
        print(f"✅ Loaded snRNA sequence ({len(seq)} bp) from {fasta_path}")
        print(f"   Base frequencies: A={vec[0]:.3f}, C={vec[1]:.3f}, G={vec[2]:.3f}, T={vec[3]:.3f}")
        return vec
    except FileNotFoundError:
        print(f"⚠️ snRNA file not found at {fasta_path}. Using uniform background.")
        return np.array([0.25, 0.25, 0.25, 0.25])
    except Exception as e:
        print(f"⚠️ Error loading snRNA: {e}. Using uniform background.")
        return np.array([0.25, 0.25, 0.25, 0.25])


def tile_snRNA_features(snrna_vec: np.ndarray, n_cells: int) -> np.ndarray:
    """
    Tile snRNA features across all cells.
    
    Args:
        snrna_vec: 4D snRNA base-frequency vector
        n_cells: Number of cells to tile for
        
    Returns:
        np.ndarray: (n_cells, 4) matrix with replicated features
    """
    snrna_tiled = np.tile(snrna_vec, (n_cells, 1))
    print(f"🔀 Tiled snRNA features to {n_cells} cells: shape {snrna_tiled.shape}")
    return snrna_tiled


def combine_features(gene_features: np.ndarray, snrna_features: np.ndarray) -> np.ndarray:
    """
    Combine gene expression features with snRNA features.
    
    Args:
        gene_features: (n_cells, n_genes) gene expression matrix
        snrna_features: (n_cells, 4) snRNA feature matrix
        
    Returns:
        np.ndarray: (n_cells, n_genes+4) combined feature matrix
    """
    combined = np.hstack([gene_features, snrna_features])
    print(f"🔗 Combined features: gene_features {gene_features.shape} + snRNA {snrna_features.shape} → {combined.shape}")
    return combined


if __name__ == "__main__":
    # Demo usage
    snrna_vec = load_snRNA_features()
    snrna_tiled = tile_snRNA_features(snrna_vec, n_cells=100)
    print(f"\nDemo: Created {snrna_tiled.shape} snRNA feature matrix")
