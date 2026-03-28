"""
Master Professional Report Generator - ASCII Safe Version
Generates comprehensive professional analysis reports
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

class MasterProfessionalReportGenerator:
    """Generates professional analysis reports"""
    
    def __init__(self, output_base="output"):
        self.output_base = Path(output_base)
        self.report_dir = self.output_base
        self.cancer_types = ['breast_cancer', 'lung_cancer', 'prostate_cancer']
        self.cancer_display = {
            'breast_cancer': 'Breast Cancer',
            'lung_cancer': 'Lung Cancer',
            'prostate_cancer': 'Prostate Cancer'
        }
    
    def load_cancer_analysis_data(self):
        """Load all results from individual cancer analyses"""
        results = {}
        
        for cancer in self.cancer_types:
            cancer_dir = self.output_base / cancer
            cancer_results = {
                'cancer_type': cancer,
                'display_name': self.cancer_display[cancer],
                'igf_summary': None,
                'dormancy_signature': None,
                'risk_scores': None,
                'drug_recommendations': None,
                'figures': []
            }
            
            # Load IGF summary
            igf_file = cancer_dir / "01_igf_pathway_analysis" / f"{cancer}_igf_summary.csv"
            if igf_file.exists():
                cancer_results['igf_summary'] = pd.read_csv(igf_file)
            
            # Load dormancy signature
            dorm_file = cancer_dir / "01_dormancy_signatures" / f"{cancer}_dormancy_signature.csv"
            if dorm_file.exists():
                cancer_results['dormancy_signature'] = pd.read_csv(dorm_file)
            
            # Load risk scores
            risk_file = cancer_dir / "02_window_of_risk" / f"{cancer}_risk_scores.csv"
            if risk_file.exists():
                cancer_results['risk_scores'] = pd.read_csv(risk_file)
            
            # Load drug recommendations
            drug_file = cancer_dir / "03_drug_repurposing" / f"{cancer}_drug_recommendations.csv"
            if drug_file.exists():
                cancer_results['drug_recommendations'] = pd.read_csv(drug_file)
            
            # Collect figures
            fig_dir = cancer_dir / "figures"
            if fig_dir.exists():
                cancer_results['figures'] = list(fig_dir.glob('*.png'))
            
            results[cancer] = cancer_results
        
        return results
    
    def create_html_report(self, analyses):
        """Create professional HTML report (ASCII safe)"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Comprehensive Cancer Analysis Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            background: #f5f7fa;
            color: #333;
            line-height: 1.6;
            padding: 20px;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            padding: 40px;
        }}
        header {{
            text-align: center;
            border-bottom: 4px solid #3498DB;
            padding-bottom: 30px;
            margin-bottom: 40px;
        }}
        h1 {{ font-size: 2.5em; color: #2C3E50; }}
        .subtitle {{ color: #7F8C8D; font-size: 1.1em; }}
        .timestamp {{ color: #95A5A6; font-style: italic; margin-top: 10px; }}
        .section {{ margin-bottom: 50px; }}
        .section-title {{ 
            font-size: 1.8em; 
            color: #2C3E50;
            border-left: 5px solid #3498DB;
            padding-left: 15px;
            margin-bottom: 25px;
            margin-top: 35px;
        }}
        .cancer-analysis {{
            background: #F8F9FA;
            border-radius: 8px;
            padding: 25px;
            margin-bottom: 30px;
            border-left: 5px solid #E74C3C;
        }}
        .cancer-title {{
            font-size: 1.5em;
            color: #E74C3C;
            margin-bottom: 15px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .metric-card {{
            background: white;
            padding: 15px;
            border-radius: 6px;
            text-align: center;
            border-top: 3px solid #3498DB;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric-label {{ font-size: 0.9em; color: #7F8C8D; text-transform: uppercase; }}
        .metric-value {{ font-size: 1.5em; font-weight: bold; color: #2C3E50; }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 6px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        th {{
            background-color: #3498DB;
            color: white;
            padding: 12px;
            text-align: left;
        }}
        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #ECF0F1;
        }}
        tr:hover {{ background-color: #F8F9FA; }}
        .figure-gallery {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(350px, 1fr));
            gap: 20px;
            margin: 25px 0;
        }}
        .figure-item {{
            background: white;
            border-radius: 6px;
            overflow: hidden;
            box-shadow: 0 4px 8px rgba(0,0,0,0.15);
        }}
        .figure-image {{ width: 100%; height: auto; display: block; }}
        .figure-caption {{
            padding: 12px;
            background: #F8F9FA;
            font-weight: 600;
            color: #2C3E50;
            font-size: 0.9em;
        }}
        .key-findings {{
            background: #E8F8F5;
            border-left: 5px solid #27AE60;
            padding: 15px;
            margin: 15px 0;
            border-radius: 4px;
        }}
        .key-findings h4 {{ color: #27AE60; margin-bottom: 10px; }}
        footer {{
            border-top: 2px solid #ECF0F1;
            padding-top: 25px;
            text-align: center;
            color: #95A5A6;
            font-size: 0.9em;
            margin-top: 50px;
        }}
        .comparative-section {{
            background: #667eea;
            color: white;
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
        }}
        .stat {{ background: #667eea; color: white; padding: 15px; border-radius: 6px; text-align: center; }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>COMPREHENSIVE MULTI-CANCER ANALYSIS REPORT</h1>
            <p class="subtitle">IGF Pathway Analysis - Clinical Risk Assessment - FDA Drug Recommendations</p>
            <p class="timestamp">Generated: {timestamp}</p>
        </header>
        
        <section class="section">
            <h2 class="section-title">EXECUTIVE SUMMARY</h2>
            <p>Integrated analysis of three major cancer types examining:</p>
            <ul>
                <li>IGF Pathway Activation: Cancer-specific signaling patterns</li>
                <li>Dormancy Signatures: Cell population composition analysis</li>
                <li>Clinical Risk Assessment: Prognostic modeling</li>
                <li>Precision Drug Recommendations: FDA-approved treatments matched to pathway activation</li>
            </ul>
        </section>
        
        <section class="section">
            <h2 class="section-title">INDIVIDUAL CANCER ANALYSES</h2>
"""
        
        # Add individual cancer analyses
        for cancer in self.cancer_types:
            analysis = analyses[cancer]
            cancer_display = self.cancer_display[cancer]
            
            html_content += f"""
            <div class="cancer-analysis">
                <div class="cancer-title">[{cancer_display}]</div>
"""
            
            # Add metrics
            if analysis['igf_summary'] is not None:
                igf = analysis['igf_summary'].iloc[0]
                html_content += """
                <div class="metrics-grid">
"""
                metrics = [
                    ('IGF Ligand', float(igf.get('IGF_Ligand_Mean', 0))),
                    ('IGF Receptor', float(igf.get('IGF_Receptor_Mean', 0))),
                    ('PI3K-AKT', float(igf.get('PI3K_AKT_Mean', 0))),
                    ('mTOR', float(igf.get('mTOR_Mean', 0))),
                    ('MAPK', float(igf.get('MAPK_Mean', 0))),
                    ('MYC', float(igf.get('MYC_Mean', 0))),
                ]
                
                for label, value in metrics:
                    html_content += f"""
                    <div class="metric-card">
                        <div class="metric-label">{label}</div>
                        <div class="metric-value">{value:.3f}</div>
                    </div>
"""
                
                html_content += """
                </div>
"""
            
            # Add key findings
            if analysis['dormancy_signature'] is not None:
                dorm_df = analysis['dormancy_signature']
                
                dormancy_pct = float(dorm_df.get('dormant_pct', [0]).values[0]) if 'dormant_pct' in dorm_df.columns else 0
                transitional_pct = float(dorm_df.get('transitional_pct', [0]).values[0]) if 'transitional_pct' in dorm_df.columns else 0
                proliferative_pct = float(dorm_df.get('proliferative_pct', [0]).values[0]) if 'proliferative_pct' in dorm_df.columns else 0
                
                html_content += f"""
                <div class="key-findings">
                    <h4>Dormancy Signature Results</h4>
                    <ul>
                        <li>Dormant Cells: {dormancy_pct:.1f}%</li>
                        <li>Transitional Cells: {transitional_pct:.1f}%</li>
                        <li>Proliferative Cells: {proliferative_pct:.1f}%</li>
                    </ul>
                </div>
"""
            
            # Add top drugs
            if analysis['drug_recommendations'] is not None:
                drugs_df = analysis['drug_recommendations'].head(5)
                html_content += """
                <div style="margin-top: 20px;">
                    <h4>Top 5 FDA-Approved Drug Recommendations</h4>
                    <table>
                        <tr>
                            <th>Drug Name</th>
                            <th>Mechanism</th>
                            <th>Score</th>
                        </tr>
"""
                for _, drug in drugs_df.iterrows():
                    drug_name = str(drug.get('drug_name', 'Unknown'))
                    mechanism = str(drug.get('mechanism', 'N/A'))
                    score = float(drug.get('score', 0))
                    html_content += f"""
                        <tr>
                            <td>{drug_name}</td>
                            <td>{mechanism}</td>
                            <td>{score:.2f}</td>
                        </tr>
"""
                html_content += """
                    </table>
                </div>
"""
            
            # Add figures
            if analysis['figures']:
                html_content += """
                <div class="figure-gallery">
"""
                for fig_path in analysis['figures']:
                    fig_name = fig_path.name
                    fig_caption = fig_name.replace('.png', '').replace('_', ' ')
                    rel_path = fig_path.relative_to(self.output_base)
                    html_content += f"""
                    <div class="figure-item">
                        <img src="{rel_path.as_posix()}" alt="{fig_caption}" class="figure-image">
                        <div class="figure-caption">{fig_caption}</div>
                    </div>
"""
                html_content += """
                </div>
"""
            
            html_content += """
            </div>
"""
        
        # Add comparative section
        html_content += """
        </section>
        
        <section class="section">
            <h2 class="section-title">COMPARATIVE ANALYSIS</h2>
            <div class="comparative-section">
                <h3>Cross-Cancer Comparison</h3>
                <p>Professional visualizations comparing IGF pathway activation, dormancy signatures, 
                and clinical outcomes across the three cancer types.</p>
            </div>
            
            <div class="figure-gallery">
                <div class="figure-item">
                    <img src="04_COMPARATIVE_ANALYSIS/IGF_Pathway_Comparison_Heatmap.png" 
                         alt="IGF Pathway Comparison" class="figure-image">
                    <div class="figure-caption">IGF Pathway Component Comparison</div>
                </div>
                <div class="figure-item">
                    <img src="04_COMPARATIVE_ANALYSIS/IGF_Activation_Comparison.png" 
                         alt="IGF Activation" class="figure-image">
                    <div class="figure-caption">Overall IGF Pathway Activation</div>
                </div>
                <div class="figure-item">
                    <img src="04_COMPARATIVE_ANALYSIS/Dormancy_Signature_Comparison.png" 
                         alt="Dormancy Comparison" class="figure-image">
                    <div class="figure-caption">Cell Population Composition</div>
                </div>
                <div class="figure-item">
                    <img src="04_COMPARATIVE_ANALYSIS/IGF_Dormancy_Correlation.png" 
                         alt="IGF-Dormancy Correlation" class="figure-image">
                    <div class="figure-caption">IGF-Dormancy Relationship</div>
                </div>
                <div class="figure-item">
                    <img src="04_COMPARATIVE_ANALYSIS/Summary_Statistics_Table.png" 
                         alt="Summary Statistics" class="figure-image">
                    <div class="figure-caption">Summary Statistics Table</div>
                </div>
            </div>
        </section>
        
        <section class="section">
            <h2 class="section-title">CLINICAL IMPLICATIONS</h2>
            <div class="key-findings">
                <h4>Key Discovery</h4>
                <p>Each cancer type demonstrates distinct IGF pathway activation correlating with 
                dormancy signatures and specific FDA-approved drug responses. This supports precision 
                medicine approaches using cancer-type and IGF-state specific combinations.</p>
            </div>
        </section>
        
        <footer>
            <p>Professional multi-cancer analysis using advanced bioinformatics pipelines incorporating 
            IGF pathway analysis, clinical risk modeling, and precision drug repurposing.</p>
        </footer>
    </div>
</body>
</html>
"""
        
        # Save HTML report
        report_file = self.report_dir / "COMPREHENSIVE_ANALYSIS_REPORT.html"
        
        # Use UTF-8 encoding
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"[OK] Professional HTML Report: {report_file.name}")
    
    def generate_reports(self):
        """Generate all professional reports"""
        print("\n" + "="*60)
        print("GENERATING PROFESSIONAL REPORTS")
        print("="*60)
        
        analyses = self.load_cancer_analysis_data()
        self.create_html_report(analyses)
        
        print("\n" + "="*60)
        print("[SUCCESS] PROFESSIONAL REPORTS GENERATION COMPLETE")
        print("="*60 + "\n")


if __name__ == "__main__":
    generator = MasterProfessionalReportGenerator()
    generator.generate_reports()
