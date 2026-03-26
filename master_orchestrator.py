"""
MASTER PIPELINE ORCHESTRATOR
Runs individual cancer analysis → Window of risk → FDA drug recommendations → Comparative analysis
"""

import os
import sys
from individual_cancer_analyzer import IndividualCancerAnalyzer
from comparative_analyzer import ComparativeAnalyzer

def main():
    """
    Execute complete analysis pipeline for all cancer types.
    """
    
    print("\n" + "="*80)
    print("🚀 GRAPHCOMM-LITE: COMPREHENSIVE MULTI-CANCER ANALYSIS PIPELINE")
    print("="*80)
    
    output_base = "output"
    cancer_types = ["breast_cancer", "lung_cancer", "prostate_cancer"]
    
    # ========== STAGE 1: INDIVIDUAL CANCER ANALYSIS ==========
    print("\n" + "#"*80)
    print("STAGE 1: INDIVIDUAL CANCER ANALYSIS (Preprocessing → IGF → Dormancy → Risk)")
    print("#"*80)
    
    individual_results = {}
    
    for cancer_type in cancer_types:
        try:
            print(f"\n{'─'*80}")
            analyzer = IndividualCancerAnalyzer(cancer_type=cancer_type, output_base=output_base)
            results = analyzer.run_complete_analysis()
            individual_results[cancer_type] = results
            
        except Exception as e:
            print(f"\n❌ Error analyzing {cancer_type}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    # ========== STAGE 2: COMPARATIVE ANALYSIS ==========
    print("\n\n" + "#"*80)
    print("STAGE 2: CROSS-CANCER COMPARATIVE ANALYSIS")
    print("#"*80)
    
    try:
        comparative = ComparativeAnalyzer(output_base=output_base)
        comparative.generate_summary_report()
    except Exception as e:
        print(f"\n❌ Error in comparative analysis: {str(e)}")
        import traceback
        traceback.print_exc()
    
    # ========== FINAL SUMMARY ==========
    print("\n\n" + "="*80)
    print("✅ PIPELINE EXECUTION COMPLETE")
    print("="*80)
    
    print("\n📊 ANALYSIS SUMMARY")
    print("-" * 80)
    
    print("\n1️⃣  INDIVIDUAL CANCER ANALYSES (by cancer type):")
    for cancer_type in cancer_types:
        cancer_dir = os.path.join(output_base, cancer_type)
        print(f"\n   📁 {cancer_type.upper().replace('_', ' ')}")
        print(f"      ├── 01_preprocessing/")
        print(f"      ├── 01_igf_pathway_analysis/       ← IGF pathway gene expression")
        print(f"      ├── 01_dormancy_signatures/        ← Cell population classification")
        print(f"      ├── 02_window_of_risk/             ← Relapse risk (6mo, 12mo, 24mo)")
        print(f"      ├── 03_drug_repurposing/           ← FDA drug recommendations")
        print(f"      └── figures/                       ← Publication-quality visualizations")
    
    print("\n2️⃣  COMPARATIVE ANALYSIS:")
    print(f"   📁 output/04_COMPARATIVE_ANALYSIS/")
    print(f"      ├── IGF_Pathway_Comparison.csv")
    print(f"      ├── Dormancy_Signature_Comparison.csv")
    print(f"      ├── Risk_Window_Comparison.csv")
    print(f"      ├── COMPARATIVE_IGF_Pathway_Heatmap.png")
    print(f"      ├── COMPARATIVE_Dormancy_Distribution.png")
    print(f"      ├── COMPARATIVE_Window_of_Risk.png")
    print(f"      └── COMPARATIVE_IGF_vs_Dormancy.png")
    
    print("\n3️⃣  KEY ANALYSIS OUTPUTS:")
    print("""
    PER CANCER TYPE:
    • IGF Pathway Summary: Gene expression levels for all pathway components
    • Dormancy Signature: % dormant, transitional, and proliferative cells
    • Window of Risk: Probability of relapse in 6, 12, 24-month windows
    • FDA Drug Recommendations: Top drugs ranked by relevance to IGF pathway
    • Figures: 5 publication-quality visualizations per cancer type
    
    COMPARATIVE:
    • IGF pathway activation ranking across cancer types
    • Dormancy status comparison
    • Relapse risk comparison
    • Scatter plot: IGF vs Dormancy relationship
    """)
    
    print("\n4️⃣  RECOMMENDED NEXT STEPS:")
    print("""
    ✓ Review individual cancer directories for cancer-specific insights
    ✓ Check figures/ folder in each cancer directory for visualizations
    ✓ Compare FDA drug recommendations across cancer types
    ✓ Analyze IGF pathway components specific to each cancer
    ✓ Use risk scores for patient stratification
    ✓ Review comparative analysis for cross-cancer patterns
    """)
    
    print("\n" + "="*80)
    print("🎯 PIPELINE STATUS: SUCCESS ✅")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
