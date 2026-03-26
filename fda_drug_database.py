"""
FDA-APPROVED DRUG DATABASE & RECOMMENDATION ENGINE
Suggests FDA-approved drugs targeting IGF pathway and dormancy reversal
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple

class FDADrugDatabase:
    """
    FDA-approved drug database with IGF pathway, dormancy, and cancer-specific targeting
    """
    
    def __init__(self):
        """Initialize FDA drug database"""
        self.drugs = self._build_database()
    
    def _build_database(self) -> pd.DataFrame:
        """
        Build comprehensive FDA drug database.
        
        Returns:
            DataFrame with drug information
        """
        drugs = [
            # ========== IGF1R INHIBITORS ==========
            {
                'drug_name': 'Linsitinib (OSI-906)',
                'target': 'IGF1R',
                'mechanism': 'Selective IGF1R tyrosine kinase inhibitor',
                'fda_status': 'Clinical Trial',
                'cancer_types': ['breast_cancer', 'lung_cancer', 'prostate_cancer'],
                'dormancy_effect': 'HIGH',
                'igf_inhibition': 'STRONG',
                'clinical_evidence': 'Phase 2 trials ongoing',
                'notes': 'Reverses dormancy by blocking IGF signaling'
            },
            {
                'drug_name': 'BMS-536924',
                'target': 'IGF1R',
                'mechanism': 'IGF1R tyrosine kinase inhibitor',
                'fda_status': 'Experimental',
                'cancer_types': ['breast_cancer', 'lung_cancer', 'prostate_cancer'],
                'dormancy_effect': 'HIGH',
                'igf_inhibition': 'STRONG',
                'clinical_evidence': 'Preclinical + Phase 1',
                'notes': 'Potent IGF pathway suppression'
            },
            {
                'drug_name': 'Figitumumab (CP-751,871)',
                'target': 'IGF1',
                'mechanism': 'Monoclonal antibody against IGF1 ligand',
                'fda_status': 'Clinical Trial',
                'cancer_types': ['lung_cancer', 'prostate_cancer'],
                'dormancy_effect': 'MODERATE',
                'igf_inhibition': 'STRONG',
                'clinical_evidence': 'Phase 2 trials',
                'notes': 'Blocks IGF1 ligand'
            },
            {
                'drug_name': 'Ganitumab (AMG 479)',
                'target': 'IGF1',
                'mechanism': 'Human monoclonal antibody targeting IGF1',
                'fda_status': 'Clinical Trial',
                'cancer_types': ['breast_cancer', 'lung_cancer'],
                'dormancy_effect': 'MODERATE',
                'igf_inhibition': 'STRONG',
                'clinical_evidence': 'Phase 2/3 trials',
                'notes': 'IGF1 sequestration'
            },
            
            # ========== PI3K/AKT/mTOR INHIBITORS ==========
            {
                'drug_name': 'Alpelisib (BYL719)',
                'target': 'PI3K-AKT-mTOR',
                'mechanism': 'Selective PI3K inhibitor (Class IA)',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer'],
                'dormancy_effect': 'MODERATE-HIGH',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for HR+/HER2- breast cancer',
                'notes': 'Synergizes with hormonal therapy'
            },
            {
                'drug_name': 'GSK2636771',
                'target': 'PI3K-beta',
                'mechanism': 'Selective PI3K-beta inhibitor',
                'fda_status': 'Clinical Trial',
                'cancer_types': ['breast_cancer', 'prostate_cancer'],
                'dormancy_effect': 'MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Phase 1/2',
                'notes': 'Selective for PI3K-beta isoform'
            },
            {
                'drug_name': 'Everolimus (Afinitor)',
                'target': 'mTOR (Complex 1)',
                'mechanism': 'mTORC1 inhibitor (rapamycin analog)',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer', 'lung_cancer', 'prostate_cancer'],
                'dormancy_effect': 'MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for ER+ breast, RCC, PNET',
                'notes': 'Blocks mTOR downstream of IGF/AKT'
            },
            {
                'drug_name': 'Temsirolimus (Torisel)',
                'target': 'mTOR',
                'mechanism': 'mTORC1 inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['lung_cancer', 'prostate_cancer'],
                'dormancy_effect': 'MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for advanced RCC',
                'notes': 'Reduces proliferation signaling'
            },
            {
                'drug_name': 'BEZ235',
                'target': 'PI3K + mTOR (dual)',
                'mechanism': 'Dual PI3K/mTOR inhibitor',
                'fda_status': 'Clinical Trial',
                'cancer_types': ['breast_cancer', 'lung_cancer'],
                'dormancy_effect': 'HIGH',
                'igf_inhibition': 'INDIRECT-STRONG',
                'clinical_evidence': 'Phase 1/2',
                'notes': 'Blocks both PI3K and mTOR simultaneously'
            },
            {
                'drug_name': 'AZD8055',
                'target': 'mTORC1 + mTORC2',
                'mechanism': 'ATP-competitive mTOR kinase inhibitor',
                'fda_status': 'Clinical Trial',
                'cancer_types': ['breast_cancer', 'lung_cancer'],
                'dormancy_effect': 'HIGH',
                'igf_inhibition': 'INDIRECT-STRONG',
                'clinical_evidence': 'Phase 1/2',
                'notes': 'Inhibits both mTOR complexes'
            },
            
            # ========== MEK/ERK (MAPK) INHIBITORS ==========
            {
                'drug_name': 'Trametinib (Mekinist)',
                'target': 'MEK1/2',
                'mechanism': 'Selective MEK1/2 inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['lung_cancer'],
                'dormancy_effect': 'LOW-MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for BRAF-mutant melanoma',
                'notes': 'Used with BRAF inhibitors'
            },
            {
                'drug_name': 'Selumetinib (AZD6244)',
                'target': 'MEK1/2',
                'mechanism': 'Allosteric MEK inhibitor',
                'fda_status': 'Clinical Trial',
                'cancer_types': ['lung_cancer', 'breast_cancer'],
                'dormancy_effect': 'LOW',
                'igf_inhibition': 'NONE',
                'clinical_evidence': 'Phase 2',
                'notes': 'MAPK pathway inhibition'
            },
            {
                'drug_name': 'Cobimetinib (Cotellic)',
                'target': 'MEK1/2',
                'mechanism': 'Selective MEK inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['lung_cancer'],
                'dormancy_effect': 'LOW-MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for BRAF-mutant melanoma',
                'notes': 'Used in combination therapy'
            },
            
            # ========== MYC PATHWAY INHIBITORS (Indirect) ==========
            {
                'drug_name': 'Palbociclib (Ibrance)',
                'target': 'CDK4/6 → MYC suppression',
                'mechanism': 'CDK4/6 inhibitor (affects cell cycle)',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer'],
                'dormancy_effect': 'HIGH',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for HR+/HER2- breast cancer',
                'notes': 'Blocks proliferation; maintains dormancy'
            },
            {
                'drug_name': 'Ribociclib (Kisqali)',
                'target': 'CDK4/6 → MYC',
                'mechanism': 'CDK4/6 inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer'],
                'dormancy_effect': 'HIGH',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for HR+/HER2- breast cancer',
                'notes': 'Cytostatic; maintains quiescence'
            },
            {
                'drug_name': 'Abemaciclib (Verzenio)',
                'target': 'CDK4/6',
                'mechanism': 'CDK4/6 inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer'],
                'dormancy_effect': 'HIGH',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for HR+/HER2- breast cancer',
                'notes': 'Can be monotherapy; cell cycle arrest'
            },
            
            # ========== HORMONE THERAPY (Breast Cancer) ==========
            {
                'drug_name': 'Tamoxifen',
                'target': 'ESR1 (ER-alpha)',
                'mechanism': 'Selective estrogen receptor modulator',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer'],
                'dormancy_effect': 'MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Standard of care ER+ breast cancer',
                'notes': 'May induce dormancy via ER antagonism'
            },
            {
                'drug_name': 'Fulvestrant (Faslodex)',
                'target': 'ESR1 degradation',
                'mechanism': 'Selective estrogen receptor degrader',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer'],
                'dormancy_effect': 'MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for ER+ breast cancer',
                'notes': 'Complete ER antagonism'
            },
            {
                'drug_name': 'Letrozole (Femara)',
                'target': 'Estrogen synthesis',
                'mechanism': 'Aromatase inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer'],
                'dormancy_effect': 'LOW-MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Standard for postmenopausal ER+ breast',
                'notes': 'Estrogen deprivation causes dormancy'
            },
            {
                'drug_name': 'Anastrozole (Arimidex)',
                'target': 'Estrogen synthesis',
                'mechanism': 'Aromatase inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer'],
                'dormancy_effect': 'LOW-MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Standard for ER+ postmenopausal breast',
                'notes': 'Reduces ESR1 ligand availability'
            },
            
            # ========== ANDROGEN THERAPY (Prostate Cancer) ==========
            {
                'drug_name': 'Bicalutamide (Casodex)',
                'target': 'AR (Androgen Receptor)',
                'mechanism': 'Androgen receptor antagonist',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['prostate_cancer'],
                'dormancy_effect': 'MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Standard for advanced prostate cancer',
                'notes': 'AR antagonism; may induce dormancy'
            },
            {
                'drug_name': 'Enzalutamide (Xtandi)',
                'target': 'AR signaling',
                'mechanism': 'Next-generation androgen receptor antagonist',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['prostate_cancer'],
                'dormancy_effect': 'MODERATE-HIGH',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for CRPC; improved survival',
                'notes': 'Potent AR inhibition; suppresses MAPK'
            },
            {
                'drug_name': 'Abiraterone (Zytiga)',
                'target': 'CYP17A1 (androgen synthesis)',
                'mechanism': 'CYP17A1 inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['prostate_cancer'],
                'dormancy_effect': 'MODERATE',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for mCRPC',
                'notes': 'Androgen deprivation therapy'
            },
            {
                'drug_name': 'Darolutamide (Nubeqa)',
                'target': 'AR',
                'mechanism': 'Next-generation androgen receptor antagonist',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['prostate_cancer'],
                'dormancy_effect': 'MODERATE-HIGH',
                'igf_inhibition': 'INDIRECT',
                'clinical_evidence': 'Approved for non-mCRPC',
                'notes': 'BBB-penetrant AR antagonist'
            },
            
            # ========== IMMUNOTHERAPY (Pan-cancer) ==========
            {
                'drug_name': 'Pembrolizumab (Keytruda)',
                'target': 'PD-1 checkpoint',
                'mechanism': 'PD-1 inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer', 'lung_cancer', 'prostate_cancer'],
                'dormancy_effect': 'LOW',
                'igf_inhibition': 'NONE',
                'clinical_evidence': 'Approved for multiple cancer types',
                'notes': 'May reduce dormancy by immune activation'
            },
            {
                'drug_name': 'Nivolumab (Opdivo)',
                'target': 'PD-1 checkpoint',
                'mechanism': 'PD-1 inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer', 'lung_cancer'],
                'dormancy_effect': 'LOW',
                'igf_inhibition': 'NONE',
                'clinical_evidence': 'Approved for NSCLC and other cancers',
                'notes': 'Immune checkpoint inhibition'
            },
            {
                'drug_name': 'Tremelimumab (Imjudo)',
                'target': 'CTLA-4 checkpoint',
                'mechanism': 'CTLA-4 inhibitor',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['lung_cancer'],
                'dormancy_effect': 'LOW',
                'igf_inhibition': 'NONE',
                'clinical_evidence': 'Approved for NSCLC',
                'notes': 'Combined with atezolizumab'
            },
            
            # ========== CHEMOTHERAPY (Reference) ==========
            {
                'drug_name': 'Paclitaxel (Taxol)',
                'target': 'Microtubules (non-specific)',
                'mechanism': 'Taxane chemotherapy',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer', 'lung_cancer'],
                'dormancy_effect': 'NONE',
                'igf_inhibition': 'NONE',
                'clinical_evidence': 'Standard chemotherapy',
                'notes': 'Cytotoxic; kills proliferating cells'
            },
            {
                'drug_name': 'Docetaxel (Taxotere)',
                'target': 'Microtubules',
                'mechanism': 'Taxane chemotherapy',
                'fda_status': 'FDA APPROVED',
                'cancer_types': ['breast_cancer', 'lung_cancer', 'prostate_cancer'],
                'dormancy_effect': 'NONE',
                'igf_inhibition': 'NONE',
                'clinical_evidence': 'Standard chemotherapy',
                'notes': 'Does not target dormant cells'
            },
        ]
        
        return pd.DataFrame(drugs)
    
    def recommend_drugs(
        self,
        cancer_type: str,
        igf_activation_level: float,
        target_dormancy: bool = True
    ) -> pd.DataFrame:
        """
        Recommend FDA-approved drugs for specific cancer type and IGF activation.
        
        Args:
            cancer_type: "breast_cancer", "lung_cancer", or "prostate_cancer"
            igf_activation_level: Mean IGF pathway activation (0-1)
            target_dormancy: Whether primary goal is to induce/maintain dormancy
        
        Returns:
            Ranked DataFrame of recommended drugs
        """
        # Filter for cancer type
        candidates = self.drugs[
            self.drugs['cancer_types'].apply(lambda x: cancer_type in x)
        ].copy()
        
        if len(candidates) == 0:
            return pd.DataFrame()
        
        # Scoring based on criteria
        candidates['score'] = 0
        
        # IGF pathway inhibition score
        igf_score_map = {'STRONG': 3, 'INDIRECT-STRONG': 3, 'INDIRECT': 1, 'NONE': 0}
        igf_strength = candidates['igf_inhibition'].map(igf_score_map)
        
        # Weight by IGF activation level
        candidates['score'] += igf_strength * (igf_activation_level + 0.5)
        
        # Dormancy effect score
        if target_dormancy:
            dormancy_score_map = {'HIGH': 3, 'MODERATE-HIGH': 2.5, 'MODERATE': 2, 'LOW-MODERATE': 1, 'LOW': 0}
            dormancy_strength = candidates['dormancy_effect'].map(dormancy_score_map)
            candidates['score'] += dormancy_strength * 2  # Double weight for dormancy
        
        # FDA approval bonus
        candidates['score'] = candidates.apply(
            lambda row: row['score'] + 2 if row['fda_status'] == 'FDA APPROVED' else row['score'],
            axis=1
        )
        
        # Sort by score
        candidates = candidates.sort_values('score', ascending=False)
        
        return candidates[['drug_name', 'target', 'mechanism', 'cancer_types', 
                          'dormancy_effect', 'igf_inhibition', 'fda_status', 'score']]


if __name__ == "__main__":
    db = FDADrugDatabase()
    
    print("\n" + "="*80)
    print("FDA DRUG DATABASE TEST")
    print("="*80)
    
    for cancer_type in ["breast_cancer", "lung_cancer", "prostate_cancer"]:
        print(f"\n{cancer_type.upper()} - HIGH IGF ACTIVATION (0.8)")
        recs = db.recommend_drugs(cancer_type, igf_activation_level=0.8, target_dormancy=True)
        print(recs.head(5).to_string(index=False))
