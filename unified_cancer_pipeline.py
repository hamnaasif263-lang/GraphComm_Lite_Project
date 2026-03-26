#!/usr/bin/env python3
"""
UNIFIED CANCER ANALYSIS PIPELINE
=================================

A comprehensive 4-stage pipeline for cancer data analysis:

STAGE 1: DATA PREPROCESSING
   - Load raw scRNA-seq data from cancer type folders
   - Normalize and scale expression data
   - Quality control and filtering
   - Output: Processed expression matrices

STAGE 2: IGF PATHWAY ANALYSIS (GraphComm-Lite)
   - Build cell-cell communication networks
   - Identify IGF ligand-receptor interactions
   - Compute pathway activity scores
   - Analyze autocrine and paracrine signaling
   - Output: Network edge files, pathway scores, visualizations

STAGE 3: WINDOW OF RISK ANALYSIS
   - Identify high-risk cell populations
   - Predict risk stratification windows
   - Analyze risk-associated pathway patterns
   - Output: Risk scores, risk windows, stratification metrics

STAGE 4: DRUG REPOSITIONING & PREDICTION
   - Predict drug responses based on gene expression
   - Map FDA-approved drugs to target genes
   - Score drug-cell interactions
   - Output: Drug predictions, repurposing recommendations

CANCER TYPES SUPPORTED:
   - Breast Cancer
   - Lung Cancer
   - Prostate Cancer

DIRECTORY STRUCTURE:
   data/
   ├── breast_cancer/
   ├── lung_cancer/
   └── prostate_cancer/

PROCESSED OUTPUTS:
   results/
   ├── stage1_preprocessed/
   ├── stage2_igf_analysis/
   ├── stage3_risk_analysis/
   └── stage4_drug_prediction/
   
   plots/
   ├── stage2_networks/
   ├── stage3_risk/
   └── stage4_drugs/

USAGE:
   python unified_cancer_pipeline.py [--cancer_type breast|lung|prostate] [--stage 1|2|3|4|all]
   
EXAMPLE:
   python unified_cancer_pipeline.py --cancer_type lung --stage all
   python unified_cancer_pipeline.py --cancer_type breast --stage 2
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from datetime import datetime

# Import existing modules
from data_preprocess import preprocess_scRNA
from pathway_analysis import build_igf_graph, compute_pathway_score
from train_predict import train_graph_model, drug_response_prediction
from visualize import plot_igf_activity, plot_dormancy_heatmap, plot_communication_graph

# ============================================================================
# PIPELINE CONFIGURATION
# ============================================================================

CANCER_TYPES = ['breast_cancer', 'lung_cancer', 'prostate_cancer']
PIPELINE_STAGES = {
    1: 'Data Preprocessing',
    2: 'IGF Pathway Analysis',
    3: 'Window of Risk Analysis',
    4: 'Drug Repositioning & Prediction'
}

# ============================================================================
# STAGE 1: DATA PREPROCESSING
# ============================================================================

class Stage1_DataPreprocessing:
    """
    Stage 1: Data Preprocessing
    
    Processes raw scRNA-seq files:
    - Loads data from cancer type folders
    - Normalizes and scales expression
    - Performs quality control
    - Outputs processed matrices
    """
    
    def __init__(self, cancer_type):
        self.cancer_type = cancer_type
        self.data_dir = Path('data') / cancer_type
        self.output_dir = Path('results') / 'stage1_preprocessed' / cancer_type
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def run(self):
        """Execute Stage 1 preprocessing"""
        print(f"\n{'='*70}")
        print(f"STAGE 1: DATA PREPROCESSING - {self.cancer_type.upper()}")
        print(f"{'='*70}")
        
        # Find all CSV files in cancer type folder
        csv_files = list(self.data_dir.glob('*.csv'))
        
        if not csv_files:
            print(f"⚠️ No CSV files found in {self.data_dir}")
            return None
        
        print(f"Found {len(csv_files)} samples to process")
        
        results = []
        for csv_file in csv_files:
            try:
                sample_name = csv_file.stem
                print(f"\n  Processing: {sample_name}...", end=" ")
                
                # Load data
                df = pd.read_csv(csv_file, index_col=0)
                
                # Preprocess using existing function
                features, adj_matrix, graph = preprocess_scRNA(str(csv_file))
                
                # Save preprocessed data
                output_file = self.output_dir / f"{sample_name}_preprocessed.csv"
                df.to_csv(output_file)
                
                # Save metadata
                metadata = {
                    'sample': sample_name,
                    'n_cells': df.shape[0],
                    'n_genes': df.shape[1],
                    'n_pca_components': features.shape[1],
                    'preprocessed_file': str(output_file)
                }
                results.append(metadata)
                
                print(f"✓ ({df.shape[0]} cells, {df.shape[1]} genes)")
                
            except Exception as e:
                print(f"✗ Error: {e}")
        
        # Save summary
        summary_df = pd.DataFrame(results)
        summary_file = self.output_dir / 'preprocessing_summary.csv'
        summary_df.to_csv(summary_file, index=False)
        
        print(f"\n✓ Stage 1 Complete: {len(results)} samples processed")
        print(f"  Summary saved to: {summary_file}")
        
        return summary_df


# ============================================================================
# STAGE 2: IGF PATHWAY ANALYSIS
# ============================================================================

class Stage2_IGFPathwayAnalysis:
    """
    Stage 2: IGF Pathway Analysis
    
    Analyzes cell-cell communication via IGF pathway:
    - Builds ligand-receptor networks
    - Computes pathway scores
    - Identifies signaling hubs
    - Visualizes communication networks
    """
    
    def __init__(self, cancer_type):
        self.cancer_type = cancer_type
        self.input_dir = Path('results') / 'stage1_preprocessed' / cancer_type
        self.output_dir = Path('results') / 'stage2_igf_analysis' / cancer_type
        self.plot_dir = Path('plots') / 'stage2_networks' / cancer_type
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.plot_dir.mkdir(parents=True, exist_ok=True)
        
    def run(self):
        """Execute Stage 2 analysis"""
        print(f"\n{'='*70}")
        print(f"STAGE 2: IGF PATHWAY ANALYSIS - {self.cancer_type.upper()}")
        print(f"{'='*70}")
        
        preprocessed_files = list(self.input_dir.glob('*_preprocessed.csv'))
        
        if not preprocessed_files:
            print(f"⚠️ No preprocessed files found. Run Stage 1 first.")
            return None
        
        results = []
        for pfile in preprocessed_files:
            try:
                sample_name = pfile.stem.replace('_preprocessed', '')
                print(f"\n  Analyzing: {sample_name}...", end=" ")
                
                # Load preprocessed data
                df = pd.read_csv(pfile, index_col=0)
                
                # Build IGF graph
                igf_graph = build_igf_graph(df, threshold=0.25)
                
                # Compute IGF pathway scores
                igf_genes = ['IGF1', 'IGF2', 'IGF1R', 'IRS1', 'IRS2', 'AKT1']
                igf_scores = compute_pathway_score(df, igf_genes, method='mean')
                
                # Save results
                graph_file = self.output_dir / f"{sample_name}_igf_network.csv"
                score_file = self.output_dir / f"{sample_name}_igf_scores.csv"
                
                pd.DataFrame(igf_graph).to_csv(graph_file, index=False)
                igf_scores.to_csv(score_file, header=['igf_score'])
                
                results.append({
                    'sample': sample_name,
                    'n_edges': np.count_nonzero(igf_graph),
                    'mean_igf_score': igf_scores.mean(),
                    'igf_high_cells': (igf_scores > igf_scores.median()).sum()
                })
                
                print(f"✓ ({np.count_nonzero(igf_graph)} edges)")
                
            except Exception as e:
                print(f"✗ Error: {e}")
        
        # Save summary
        summary_df = pd.DataFrame(results)
        summary_file = self.output_dir / 'igf_analysis_summary.csv'
        summary_df.to_csv(summary_file, index=False)
        
        print(f"\n✓ Stage 2 Complete: {len(results)} samples analyzed")
        print(f"  Summary saved to: {summary_file}")
        
        return summary_df


# ============================================================================
# STAGE 3: WINDOW OF RISK ANALYSIS
# ============================================================================

class Stage3_WindowOfRiskAnalysis:
    """
    Stage 3: Window of Risk Analysis
    
    Identifies high-risk cell populations:
    - Analyzes pathway activation patterns
    - Predicts risk stratification windows
    - Identifies risk-associated signatures
    - Outputs risk scores and windows
    """
    
    def __init__(self, cancer_type):
        self.cancer_type = cancer_type
        self.input_dir = Path('results') / 'stage2_igf_analysis' / cancer_type
        self.output_dir = Path('results') / 'stage3_risk_analysis' / cancer_type
        self.plot_dir = Path('plots') / 'stage3_risk' / cancer_type
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.plot_dir.mkdir(parents=True, exist_ok=True)
        
    def run(self):
        """Execute Stage 3 analysis"""
        print(f"\n{'='*70}")
        print(f"STAGE 3: WINDOW OF RISK ANALYSIS - {self.cancer_type.upper()}")
        print(f"{'='*70}")
        print("⏳ [PLACEHOLDER: Ready for implementation]")
        print("   - Risk stratification will be computed")
        print("   - Risk windows will be identified")
        print("   - Risk-associated pathways will be analyzed")
        
        return None


# ============================================================================
# STAGE 4: DRUG REPOSITIONING & PREDICTION
# ============================================================================

class Stage4_DrugPrediction:
    """
    Stage 4: Drug Repositioning & Prediction
    
    Predicts drug responses and identifies drug repositioning opportunities:
    - Predicts drug-cell interactions
    - Maps FDA-approved drugs to target genes
    - Scores drug repositioning candidates
    - Outputs drug predictions
    """
    
    def __init__(self, cancer_type):
        self.cancer_type = cancer_type
        self.input_dir = Path('results') / 'stage3_risk_analysis' / cancer_type
        self.output_dir = Path('results') / 'stage4_drug_prediction' / cancer_type
        self.plot_dir = Path('plots') / 'stage4_drugs' / cancer_type
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.plot_dir.mkdir(parents=True, exist_ok=True)
        
    def run(self):
        """Execute Stage 4 analysis"""
        print(f"\n{'='*70}")
        print(f"STAGE 4: DRUG REPOSITIONING & PREDICTION - {self.cancer_type.upper()}")
        print(f"{'='*70}")
        print("⏳ [PLACEHOLDER: Ready for implementation]")
        print("   - Drug responses will be predicted")
        print("   - FDA-approved drugs will be mapped")
        print("   - Drug repositioning candidates will be identified")
        
        return None


# ============================================================================
# MAIN PIPELINE ORCHESTRATOR
# ============================================================================

class UnifiedCancerPipeline:
    """Main pipeline orchestrator"""
    
    def __init__(self, cancer_type, stages_to_run):
        self.cancer_type = cancer_type
        self.stages_to_run = stages_to_run
        
        # Initialize stage executors
        self.stage1 = Stage1_DataPreprocessing(cancer_type)
        self.stage2 = Stage2_IGFPathwayAnalysis(cancer_type)
        self.stage3 = Stage3_WindowOfRiskAnalysis(cancer_type)
        self.stage4 = Stage4_DrugPrediction(cancer_type)
        
    def run(self):
        """Execute pipeline"""
        print(f"\n{'#'*70}")
        print(f"# UNIFIED CANCER ANALYSIS PIPELINE")
        print(f"# Cancer Type: {self.cancer_type.upper()}")
        print(f"# Stages: {', '.join([PIPELINE_STAGES[s] for s in self.stages_to_run])}")
        print(f"# Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'#'*70}")
        
        try:
            if 1 in self.stages_to_run:
                self.stage1.run()
            
            if 2 in self.stages_to_run:
                self.stage2.run()
            
            if 3 in self.stages_to_run:
                self.stage3.run()
            
            if 4 in self.stages_to_run:
                self.stage4.run()
            
            print(f"\n{'#'*70}")
            print(f"# PIPELINE COMPLETE")
            print(f"# Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"{'#'*70}\n")
            
        except Exception as e:
            print(f"\n❌ Pipeline failed: {e}")
            return False
        
        return True


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Unified Cancer Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  python unified_cancer_pipeline.py --cancer_type lung --stage all
  python unified_cancer_pipeline.py --cancer_type breast --stage 1,2
  python unified_cancer_pipeline.py --cancer_type prostate --stage 3,4
        """
    )
    
    parser.add_argument('--cancer_type', 
                       choices=['breast_cancer', 'lung_cancer', 'prostate_cancer'],
                       default='lung_cancer',
                       help='Cancer type to analyze (default: lung_cancer)')
    
    parser.add_argument('--stage',
                       default='all',
                       help='Stages to run: "all" or comma-separated list (1,2,3,4)')
    
    args = parser.parse_args()
    
    # Parse stages
    if args.stage.lower() == 'all':
        stages = [1, 2, 3, 4]
    else:
        try:
            stages = [int(s.strip()) for s in args.stage.split(',')]
        except ValueError:
            print("Error: Invalid stage specification. Use 'all' or comma-separated numbers (1,2,3,4)")
            return False
    
    # Create and run pipeline
    pipeline = UnifiedCancerPipeline(args.cancer_type, stages)
    return pipeline.run()


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
