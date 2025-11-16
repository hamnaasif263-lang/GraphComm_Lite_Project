import numpy as np

ligand_receptor_pairs = [
    ("IGF1", "IGF1R"),
    ("IGF2", "IGF1R"),
    ("IGF1", "INSR")
]

def build_igf_graph(df, threshold=0.2):
    n = len(df)
    A = np.zeros((n, n))
    for ligand, receptor in ligand_receptor_pairs:
        if ligand not in df.columns or receptor not in df.columns:
            continue
        for i in range(n):
            for j in range(n):
                if df.iloc[i][ligand] > threshold and df.iloc[j][receptor] > threshold:
                    A[i, j] += 1
    return (A > 0).astype(float)
