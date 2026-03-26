"""
IGF PATHWAY ANALYSIS MODULE
Analyze IGF ligands, receptors, and downstream signaling per cancer type
Identify dormancy signatures and pathway activation
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Tuple
import warnings
warnings.filterwarnings('ignore')

# IGF PATHWAY GENES
IGF_PATHWAY_GENES = {
    'ligands': ['IGF1', 'IGF2', 'INSULIN'],
    'receptors': ['IGF1R', 'INSR', 'IR', 'IGF2R'],
    'pi3k_akt': ['PIK3CA', 'PIK3CB', 'AKT1', 'AKT2', 'AKT3', 'PKB'],
    'mtor': ['MTOR', 'RPTOR', 'RICTOR', 'PTEN', 'TSC1', 'TSC2'],
    'mapk': ['RAF1', 'MAP2K1', 'ERK1', 'ERK2', 'MAPK1', 'MAPK3'],
    'myc': ['MYC', 'MYCN', 'MYCL', 'CCND1', 'CDKs'],
    'dormancy': ['CDKN1B', 'p27', 'CDKN1A', 'p21', 'CDKN2A', 'p16'],
    'proliferation': ['PCNA', 'KI67', 'MKI67', 'CCNA2', 'CCNB1'],
    'esr1': ['ESR1', 'ER-alpha'],  # Estrogen receptor (breast)
    'ar': ['AR', 'Androgen-receptor'],  # Androgen receptor (prostate)
}

class IGFPathwayAnalyzer:
    """
    Comprehensive IGF pathway analysis including:
    - Gene expression quantification
    - Pathway activation scoring
    - Cancer-type specific patterns
    - Dormancy marker identification
    """
    
    def __init__(self, cancer_type: str = "breast_cancer"):
        """
        Initialize IGF pathway analyzer.
        Args:
            cancer_type: "breast_cancer", "lung_cancer", or "prostate_cancer"
        """
        self.cancer_type = cancer_type
        self.pathway_genes = IGF_PATHWAY_GENES
        self.gene_expression = None
        self.pathway_scores = None
        self.dormancy_signature = None
        
    def generate_cancer_specific_expression(
        self,
        single_cell_data: pd.DataFrame,
        n_genes_per_category: int = 100
    ) -> pd.DataFrame:
        """
        Generate realistic cancer-type specific IGF pathway expression.
        
        Args:
            single_cell_data: Input single-cell data
            n_genes_per_category: Number of pathway genes per category
        
        Returns:
            Expression matrix with pathway genes
        """
        np.random.seed(hash(self.cancer_type) % 2**32)
        
        n_cells = len(single_cell_data)
        expression_data = pd.DataFrame({'cell_index': range(n_cells)})
        
        # Generate all pathway gene expressions with cancer-specific patterns
        for category, genes in self.pathway_genes.items():
            for gene in genes:
                # Base expression
                base_expr = np.random.gamma(2, 2, n_cells)
                
                # Cancer-type specific modulation
                if self.cancer_type == "breast_cancer":
                    if category == 'ligands':  # High IGF ligands
                        base_expr *= 1.8
                    elif gene in ['ESR1', 'ER-alpha']:
                        base_expr *= 2.0  # High ER in breast
                    elif category == 'pi3k_akt':
                        base_expr *= 1.5
                
                elif self.cancer_type == "lung_cancer":
                    if category == 'mapk':  # High MAPK in lung
                        base_expr *= 1.6
                    elif category == 'ligands':
                        base_expr *= 1.3
                    elif category == 'myc':
                        base_expr *= 1.7  # High MYC in lung
                
                elif self.cancer_type == "prostate_cancer":
                    if gene in ['AR', 'Androgen-receptor']:
                        base_expr *= 2.2  # Very high AR in prostate
                    elif category == 'myc':
                        base_expr *= 1.2
                    elif category == 'dormancy':
                        base_expr *= 0.8  # Lower dormancy initially
                
                # Add noise
                expression_data[gene] = base_expr + np.random.normal(0, 0.5, n_cells)
                expression_data[gene] = expression_data[gene].clip(lower=0)
        
        self.gene_expression = expression_data
        return expression_data
    
    def calculate_pathway_scores(self) -> pd.DataFrame:
        """
        Calculate pathway activation scores for each cell.
        
        Returns:
            DataFrame with pathway scores per cell
        """
        if self.gene_expression is None:
            raise ValueError("Must call generate_cancer_specific_expression first")
        
        expr = self.gene_expression.copy()
        scores = pd.DataFrame({'cell_index': expr['cell_index']})
        
        # IGF LIGAND SCORE (IGF1, IGF2 normalized)
        igf_ligands = ['IGF1', 'IGF2']
        available_ligands = [g for g in igf_ligands if g in expr.columns]
        if available_ligands:
            scores['igf_ligand_score'] = expr[available_ligands].mean(axis=1) / expr[available_ligands].std(axis=1).max()
        
        # IGF RECEPTOR SCORE
        igf_receptors = ['IGF1R', 'INSR']
        available_receptors = [g for g in igf_receptors if g in expr.columns]
        if available_receptors:
            scores['igf_receptor_score'] = expr[available_receptors].mean(axis=1) / expr[available_receptors].std(axis=1).max()
        
        # PI3K-AKT PATHWAY SCORE
        pi3k_akt = [g for g in self.pathway_genes['pi3k_akt'] if g in expr.columns]
        if pi3k_akt:
            scores['pi3k_akt_score'] = expr[pi3k_akt].mean(axis=1) / expr[pi3k_akt].std(axis=1).max()
        
        # mTOR PATHWAY SCORE
        mtor = [g for g in self.pathway_genes['mtor'] if g in expr.columns]
        if mtor:
            scores['mtor_score'] = expr[mtor].mean(axis=1) / expr[mtor].std(axis=1).max()
        
        # MAPK PATHWAY SCORE
        mapk = [g for g in self.pathway_genes['mapk'] if g in expr.columns]
        if mapk:
            scores['mapk_score'] = expr[mapk].mean(axis=1) / expr[mapk].std(axis=1).max()
        
        # MYC SCORE
        myc_genes = [g for g in self.pathway_genes['myc'] if g in expr.columns]
        if myc_genes:
            scores['myc_score'] = expr[myc_genes].mean(axis=1) / expr[myc_genes].std(axis=1).max()
        
        # DORMANCY SCORE (NEGATIVE = DORMANT)
        dormancy_genes = [g for g in self.pathway_genes['dormancy'] if g in expr.columns]
        if dormancy_genes:
            scores['dormancy_score'] = expr[dormancy_genes].mean(axis=1) / expr[dormancy_genes].std(axis=1).max()
        
        # PROLIFERATION SCORE
        prolif_genes = [g for g in self.pathway_genes['proliferation'] if g in expr.columns]
        if prolif_genes:
            scores['proliferation_score'] = expr[prolif_genes].mean(axis=1) / expr[prolif_genes].std(axis=1).max()
        
        # CANCER-SPECIFIC SCORES
        if self.cancer_type == "breast_cancer" and 'ESR1' in expr.columns:
            scores['esr1_score'] = expr['ESR1'] / expr['ESR1'].std().max() if expr['ESR1'].std() > 0 else 0
        
        if self.cancer_type == "prostate_cancer" and 'AR' in expr.columns:
            scores['ar_score'] = expr['AR'] / expr['AR'].std().max() if expr['AR'].std() > 0 else 0
        
        # OVERALL IGF PATHWAY ACTIVATION
        pathway_cols = [c for c in scores.columns if c.endswith('_score') and c != 'dormancy_score']
        if pathway_cols:
            scores['igf_pathway_activation'] = scores[pathway_cols].mean(axis=1)
        
        self.pathway_scores = scores
        return scores
    
    def identify_dormancy_signature(self, threshold: float = 0.7) -> Dict:
        """
        Identify cells with dormancy signature.
        
        Dormant cells: High dormancy markers (CDKN1B, p27, p16, p21)
                      + Low proliferation
                      + Variable IGF signaling
        
        Args:
            threshold: Score threshold for classification
        
        Returns:
            Dictionary with dormancy signature details
        """
        if self.pathway_scores is None:
            raise ValueError("Must call calculate_pathway_scores first")
        
        scores = self.pathway_scores.copy()
        
        # Normalize scores to 0-1 range
        for col in scores.columns:
            if col.endswith('_score'):
                min_val = scores[col].min()
                max_val = scores[col].max()
                if max_val > min_val:
                    scores[col] = (scores[col] - min_val) / (max_val - min_val)
        
        # DORMANT: High dormancy + Low proliferation + Low IGF
        dormant_mask = (
            (scores.get('dormancy_score', 0) > threshold) &
            (scores.get('proliferation_score', 1) < (1 - threshold)) &
            (scores.get('igf_pathway_activation', 0) < 0.5)
        )
        
        # PROLIFERATIVE: Low dormancy + High proliferation
        prolif_mask = (
            (scores.get('dormancy_score', 0) < (1 - threshold)) &
            (scores.get('proliferation_score', 0) > threshold)
        )
        
        # TRANSITIONAL (awakening): Moderate dormancy + Starting proliferation + Rising IGF
        transit_mask = (
            (scores.get('dormancy_score', 0).between(0.3, 0.7)) &
            (scores.get('proliferation_score', 0).between(0.3, 0.7)) &
            (scores.get('igf_pathway_activation', 0) > 0.4)
        )
        
        signature = {
            'dormant_cells': dormant_mask.sum(),
            'proliferative_cells': prolif_mask.sum(),
            'transitional_cells': transit_mask.sum(),
            'dormant_pct': (dormant_mask.sum() / len(scores)) * 100,
            'proliferative_pct': (prolif_mask.sum() / len(scores)) * 100,
            'transitional_pct': (transit_mask.sum() / len(scores)) * 100,
            'cancer_type': self.cancer_type,
            'mean_igf_activation': scores.get('igf_pathway_activation', pd.Series([0])).mean(),
            'mean_dormancy': scores.get('dormancy_score', pd.Series([0])).mean(),
            'mean_proliferation': scores.get('proliferation_score', pd.Series([0])).mean(),
        }
        
        self.dormancy_signature = signature
        return signature
    
    def get_summary_statistics(self) -> pd.DataFrame:
        """
        Get summary statistics of pathway analysis.
        """
        if self.pathway_scores is None:
            raise ValueError("Must call calculate_pathway_scores first")
        
        scores = self.pathway_scores.copy()
        
        summary = pd.DataFrame({
            'Cancer_Type': [self.cancer_type],
            'IGF_Ligand_Mean': [scores.get('igf_ligand_score', pd.Series([0])).mean()],
            'IGF_Ligand_Std': [scores.get('igf_ligand_score', pd.Series([0])).std()],
            'IGF_Receptor_Mean': [scores.get('igf_receptor_score', pd.Series([0])).mean()],
            'IGF_Receptor_Std': [scores.get('igf_receptor_score', pd.Series([0])).std()],
            'PI3K_AKT_Mean': [scores.get('pi3k_akt_score', pd.Series([0])).mean()],
            'mTOR_Mean': [scores.get('mtor_score', pd.Series([0])).mean()],
            'MAPK_Mean': [scores.get('mapk_score', pd.Series([0])).mean()],
            'MYC_Mean': [scores.get('myc_score', pd.Series([0])).mean()],
            'Dormancy_Mean': [scores.get('dormancy_score', pd.Series([0])).mean()],
            'Proliferation_Mean': [scores.get('proliferation_score', pd.Series([0])).mean()],
            'IGF_Pathway_Activation': [scores.get('igf_pathway_activation', pd.Series([0])).mean()],
        })
        
        return summary


if __name__ == "__main__":
    # Test the analyzer
    print("Testing IGF Pathway Analyzer...")
    
    for cancer in ["breast_cancer", "lung_cancer", "prostate_cancer"]:
        print(f"\n{cancer.upper()}:")
        
        # Create dummy single-cell data
        dummy_data = pd.DataFrame({
            'cell_index': range(1000),
            'patient_id': np.repeat(range(10), 100),
        })
        
        analyzer = IGFPathwayAnalyzer(cancer_type=cancer)
        expr = analyzer.generate_cancer_specific_expression(dummy_data)
        scores = analyzer.calculate_pathway_scores()
        sig = analyzer.identify_dormancy_signature()
        
        print(f"  Dormant cells: {sig['dormant_pct']:.1f}%")
        print(f"  Proliferative cells: {sig['proliferative_pct']:.1f}%")
        print(f"  Transitional cells: {sig['transitional_pct']:.1f}%")
        print(f"  IGF Pathway Activation: {sig['mean_igf_activation']:.2f}")
