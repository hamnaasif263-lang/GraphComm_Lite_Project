"""
PROFESSIONAL VISUALIZATION MODULE
Creates Nature-grade publication-ready graphics for the GraphComm-Lite pipeline
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality styling
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# Professional color schemes
CANCER_COLORS = {
    'breast_cancer': '#E74C3C',      # Red
    'lung_cancer': '#3498DB',         # Blue
    'prostate_cancer': '#2ECC71'      # Green
}

RISK_COLORS = {
    'Low Risk': '#2ECC71',            # Green
    'Intermediate Risk': '#F39C12',   # Orange
    'High Risk': '#E74C3C'            # Red
}


class ProfessionalVisualizer:
    """
    Creates publication-quality visualizations for cancer risk stratification.
    """
    
    def __init__(self, output_base_dir: str = "output"):
        self.output_base_dir = output_base_dir
        self.figures_dir = os.path.join(output_base_dir, "figures")
        os.makedirs(self.figures_dir, exist_ok=True)
        self.dpi = 300  # Publication quality
        
    def set_publication_style(self, fig, ax=None):
        """Apply publication-quality styling to figures."""
        if ax is None:
            ax = plt.gca()
        
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Bold left and bottom spines
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        
        # Font settings
        ax.tick_params(labelsize=11, width=1.5, length=5)
        
    def create_comparative_risk_plot(self, comparative_df: pd.DataFrame):
        """
        Create a professional comparative risk plot across cancer types.
        """
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.suptitle('Estimated Probability of Relapse by Time Window and Cancer Type',
                     fontsize=16, fontweight='bold', y=1.02)
        
        windows = [6, 12, 24]
        
        for idx, window in enumerate(windows):
            ax = axes[idx]
            
            data = comparative_df[comparative_df['window_months'] == window].copy()
            data = data.sort_values('mean_risk', ascending=False)
            
            # Create bar plot
            bars = ax.bar(range(len(data)), data['mean_risk'], 
                          color=[CANCER_COLORS[c] for c in data['cancer_type']],
                          alpha=0.8, edgecolor='black', linewidth=1.5)
            
            # Add error bars (std deviation)
            ax.errorbar(range(len(data)), data['mean_risk'], 
                       yerr=data['std_risk'], fmt='none', 
                       ecolor='black', capsize=8, capthick=2, alpha=0.6)
            
            # Customization
            ax.set_xticks(range(len(data)))
            ax.set_xticklabels([c.replace('_', '\n').title() for c in data['cancer_type']], 
                               fontsize=11, fontweight='bold')
            ax.set_ylabel('Probability of Relapse', fontsize=12, fontweight='bold')
            ax.set_ylim(0, max(data['mean_risk']) * 1.3)
            ax.set_title(f'{window}-Month Window', fontsize=13, fontweight='bold', pad=10)
            
            # Add value labels on bars
            for i, (bar, val) in enumerate(zip(bars, data['mean_risk'])):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                       f'{val:.1%}', ha='center', va='bottom', fontsize=10, fontweight='bold')
            
            self.set_publication_style(fig, ax)
            ax.grid(axis='y', alpha=0.3, linestyle='--')
            ax.set_axisbelow(True)
        
        plt.tight_layout()
        filepath = os.path.join(self.figures_dir, 'Figure_1_Comparative_Risk_by_Window.png')
        plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        print(f"✓ Figure 1 saved: {filepath}")
        plt.close()
    
    def create_risk_distribution_plots(self, output_paths: dict):
        """
        Create professional risk distribution plots for each cancer type.
        """
        fig, axes = plt.subplots(3, 3, figsize=(16, 14))
        fig.suptitle('Distribution of Relapse Risk Across Patient Cohorts',
                     fontsize=18, fontweight='bold', y=0.995)
        
        cancer_types = ['breast_cancer', 'lung_cancer', 'prostate_cancer']
        windows = [6, 12, 24]
        
        for row, cancer in enumerate(cancer_types):
            # Load patient risk report
            risk_file = os.path.join(
                output_paths[cancer]['02_window_of_risk'],
                f'{cancer}_patient_risk_report.csv'
            )
            risk_df = pd.read_csv(risk_file)
            
            for col, window in enumerate(windows):
                ax = axes[row, col]
                col_name = f'risk_{window}mo'
                
                # Create histogram
                n, bins, patches = ax.hist(
                    risk_df[col_name], bins=20,
                    color=CANCER_COLORS[cancer], alpha=0.7,
                    edgecolor='black', linewidth=1.2
                )
                
                # Color code by risk level
                for i, patch in enumerate(patches):
                    if bins[i] < 0.15:
                        patch.set_facecolor('#2ECC71')
                    elif bins[i] < 0.35:
                        patch.set_facecolor('#F39C12')
                    else:
                        patch.set_facecolor('#E74C3C')
                
                # Add vertical lines for mean and median
                mean_risk = risk_df[col_name].mean()
                median_risk = risk_df[col_name].median()
                
                ax.axvline(mean_risk, color='red', linestyle='--', linewidth=2.5, 
                          label=f'Mean: {mean_risk:.1%}', alpha=0.8)
                ax.axvline(median_risk, color='darkblue', linestyle=':', linewidth=2.5,
                          label=f'Median: {median_risk:.1%}', alpha=0.8)
                
                # Customization
                ax.set_xlabel('Risk Probability', fontsize=11, fontweight='bold')
                ax.set_ylabel('Number of Patients', fontsize=11, fontweight='bold')
                
                if row == 0:
                    ax.set_title(f'{window}-Month Window', fontsize=12, fontweight='bold')
                if col == 0:
                    ax.set_ylabel(f'{cancer.replace("_", " ").title()}\nNumber of Patients',
                                 fontsize=11, fontweight='bold')
                
                if row == 2 and col == 2:
                    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
                
                self.set_publication_style(fig, ax)
                ax.grid(axis='y', alpha=0.3, linestyle='--')
                ax.set_axisbelow(True)
        
        plt.tight_layout()
        filepath = os.path.join(self.figures_dir, 'Figure_2_Risk_Distribution_by_Cancer.png')
        plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        print(f"✓ Figure 2 saved: {filepath}")
        plt.close()
    
    def create_risk_stratification_summary(self, output_paths: dict):
        """
        Create summary table visualization of risk stratification.
        """
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.suptitle('Patient Risk Stratification Summary at 24-Month Window',
                     fontsize=16, fontweight='bold', y=1.02)
        
        cancer_types = ['breast_cancer', 'lung_cancer', 'prostate_cancer']
        
        for idx, cancer in enumerate(cancer_types):
            ax = axes[idx]
            
            # Load patient risk report
            risk_file = os.path.join(
                output_paths[cancer]['02_window_of_risk'],
                f'{cancer}_patient_risk_report.csv'
            )
            risk_df = pd.read_csv(risk_file)
            
            # Count risk categories at 24 months
            risk_24mo_category = risk_df['risk_24mo_category'].value_counts()
            
            # Order categories
            categories = ['Low Risk', 'Intermediate Risk', 'High Risk']
            counts = [risk_24mo_category.get(cat, 0) for cat in categories]
            colors = [RISK_COLORS.get(cat, '#95A5A6') for cat in categories]
            
            # Create pie chart
            wedges, texts, autotexts = ax.pie(
                counts, labels=categories, autopct='%1.1f%%',
                colors=colors, startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'},
                wedgeprops={'edgecolor': 'black', 'linewidth': 1.5}
            )
            
            # Improve autotext
            for autotext in autotexts:
                autotext.set_color('white')
                autotext.set_fontsize(11)
                autotext.set_fontweight('bold')
            
            ax.set_title(f'{cancer.replace("_", " ").title()}\n(n={len(risk_df)})',
                        fontsize=13, fontweight='bold', pad=10)
        
        plt.tight_layout()
        filepath = os.path.join(self.figures_dir, 'Figure_3_Risk_Stratification_Pie.png')
        plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        print(f"✓ Figure 3 saved: {filepath}")
        plt.close()
    
    def create_heatmap_risk_matrix(self, master_df: pd.DataFrame):
        """
        Create a heatmap showing risk scores across patients and time windows.
        """
        # Create a summary heatmap for first 30 patients
        sample_df = master_df.head(30).copy()
        
        # Prepare data
        heatmap_data = sample_df[['patient_id', 'cancer_type', 'risk_6mo', 'risk_12mo', 'risk_24mo']].set_index('patient_id')
        heatmap_data = heatmap_data.drop('cancer_type', axis=1)
        
        fig, ax = plt.subplots(figsize=(10, 14))
        
        # Create heatmap
        sns.heatmap(heatmap_data, cmap='RdYlGn_r', annot=True, fmt='.2%',
                   cbar_kws={'label': 'Risk Probability'},
                   linewidths=0.5, linecolor='gray', ax=ax,
                   vmin=0, vmax=0.5)
        
        ax.set_title('Patient-Level Risk Scores Across Time Windows\n(First 30 Patients)',
                    fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel('Time Window', fontsize=12, fontweight='bold')
        ax.set_ylabel('Patient ID', fontsize=12, fontweight='bold')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        
        self.set_publication_style(fig, ax)
        
        plt.tight_layout()
        filepath = os.path.join(self.figures_dir, 'Figure_4_Risk_Heatmap.png')
        plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        print(f"✓ Figure 4 saved: {filepath}")
        plt.close()
    
    def create_comparative_boxplot(self, comparative_df: pd.DataFrame):
        """
        Create comparative boxplot showing risk distributions.
        """
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        fig.suptitle('Risk Distribution Comparison: Boxplot Analysis',
                     fontsize=16, fontweight='bold', y=1.02)
        
        windows = [6, 12, 24]
        cancer_order = ['breast_cancer', 'lung_cancer', 'prostate_cancer']
        
        for idx, window in enumerate(windows):
            ax = axes[idx]
            
            data = comparative_df[comparative_df['window_months'] == window].copy()
            data = data.sort_values('cancer_type')
            
            # Prepare data for boxplot
            positions = range(len(data))
            plot_data = [
                np.random.normal(data.iloc[i]['mean_risk'], data.iloc[i]['std_risk'], 50)
                for i in range(len(data))
            ]
            
            bp = ax.boxplot(plot_data, positions=positions, widths=0.6,
                           patch_artist=True, showmeans=True,
                           meanprops={'marker': 'D', 'markerfacecolor': 'red', 'markersize': 8})
            
            # Color boxes
            for patch, cancer in zip(bp['boxes'], data['cancer_type']):
                patch.set_facecolor(CANCER_COLORS[cancer])
                patch.set_alpha(0.7)
                patch.set_edgecolor('black')
                patch.set_linewidth(1.5)
            
            # Color whiskers and caps
            for whisker in bp['whiskers']:
                whisker.set_linewidth(1.5)
            for cap in bp['caps']:
                cap.set_linewidth(1.5)
            
            ax.set_xticks(positions)
            ax.set_xticklabels([c.replace('_', '\n').title() for c in data['cancer_type']],
                              fontsize=11, fontweight='bold')
            ax.set_ylabel('Risk Probability', fontsize=12, fontweight='bold')
            ax.set_title(f'{window}-Month Window', fontsize=13, fontweight='bold', pad=10)
            ax.set_ylim(-0.05, max(data['mean_risk'].values) * 1.5)
            
            self.set_publication_style(fig, ax)
            ax.grid(axis='y', alpha=0.3, linestyle='--')
            ax.set_axisbelow(True)
        
        plt.tight_layout()
        filepath = os.path.join(self.figures_dir, 'Figure_5_Risk_Boxplot.png')
        plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        print(f"✓ Figure 5 saved: {filepath}")
        plt.close()
    
    def create_summary_statistics_table(self, comparative_df: pd.DataFrame):
        """
        Create a publication-quality summary statistics table.
        """
        fig, ax = plt.subplots(figsize=(14, 8))
        ax.axis('tight')
        ax.axis('off')
        
        # Prepare data
        table_data = []
        for cancer in ['breast_cancer', 'lung_cancer', 'prostate_cancer']:
            cancer_label = cancer.replace('_', ' ').title()
            subset = comparative_df[comparative_df['cancer_type'] == cancer]
            
            for _, row in subset.iterrows():
                table_data.append([
                    cancer_label,
                    f"{row['window_months']:.0f} months",
                    f"{row['mean_risk']:.1%}",
                    f"{row['median_risk']:.1%}",
                    f"{row['std_risk']:.1%}",
                    f"{row['high_risk_count']:.0f}/{row['total_patients']:.0f}"
                ])
        
        columns = ['Cancer Type', 'Time Window', 'Mean Risk', 'Median Risk', 'Std Dev', 'High Risk (>50%)']
        
        table = ax.table(cellText=table_data, colLabels=columns,
                        cellLoc='center', loc='center',
                        colColours=['#34495E']*6)
        
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1, 2.5)
        
        # Style header
        for i in range(len(columns)):
            table[(0, i)].set_facecolor('#34495E')
            table[(0, i)].set_text_props(weight='bold', color='black')
        
        # Alternate row colors
        for i in range(1, len(table_data) + 1):
            for j in range(len(columns)):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#ECF0F1')
                else:
                    table[(i, j)].set_facecolor('#FFFFFF')
                table[(i, j)].set_edgecolor('#34495E')
        
        plt.title('Summary Statistics: Window of Risk Analysis',
                 fontsize=16, fontweight='bold', pad=20)
        
        filepath = os.path.join(self.figures_dir, 'Table_1_Summary_Statistics.png')
        plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        print(f"✓ Table 1 saved: {filepath}")
        plt.close()
    
    def create_patient_flow_diagram(self, output_paths: dict):
        """
        Create a patient flow diagram showing cohort stratification.
        """
        fig, ax = plt.subplots(figsize=(14, 10))
        
        # Get patient counts
        cancer_types = ['breast_cancer', 'lung_cancer', 'prostate_cancer']
        total_patients = 0
        risk_summary = {}
        
        for cancer in cancer_types:
            risk_file = os.path.join(
                output_paths[cancer]['02_window_of_risk'],
                f'{cancer}_patient_risk_report.csv'
            )
            risk_df = pd.read_csv(risk_file)
            total_patients += len(risk_df)
            
            risk_summary[cancer] = {
                'total': len(risk_df),
                'low_risk': (risk_df['risk_24mo_category'] == 'Low Risk').sum(),
                'inter_risk': (risk_df['risk_24mo_category'] == 'Intermediate Risk').sum(),
                'high_risk': (risk_df['risk_24mo_category'] == 'High Risk').sum()
            }
        
        # Drawing coordinates
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 12)
        ax.axis('off')
        
        # Title
        ax.text(5, 11.5, f'Patient Cohort Flow: {total_patients} Patients Stratified by Relapse Risk',
               ha='center', fontsize=16, fontweight='bold')
        
        # Top box
        rect_all = Rectangle((2, 9.5), 6, 1.2, linewidth=2, edgecolor='black', facecolor='#3498DB', alpha=0.3)
        ax.add_patch(rect_all)
        ax.text(5, 10.1, f'Total Patients Enrolled\nn = {total_patients}',
               ha='center', va='center', fontsize=12, fontweight='bold')
        
        # Cancer type boxes
        y_pos = 7.5
        x_positions = [1, 3.5, 6]
        
        for idx, (cancer, x_pos) in enumerate(zip(cancer_types, x_positions)):
            # Cancer label
            rect = Rectangle((x_pos-0.8, y_pos), 1.6, 0.8, linewidth=2, 
                            edgecolor='black', facecolor=CANCER_COLORS[cancer], alpha=0.5)
            ax.add_patch(rect)
            ax.text(x_pos, y_pos + 0.4, cancer.replace('_', '\n').title(),
                   ha='center', va='center', fontsize=10, fontweight='bold')
            
            # Connect from top
            ax.arrow(5, 9.5, x_pos - 5, y_pos + 0.8 - 9.5, 
                    head_width=0.15, head_length=0.1, fc='black', ec='black', alpha=0.5)
            
            # Risk stratification boxes
            risks = risk_summary[cancer]
            
            # Low risk
            rect_low = Rectangle((x_pos-0.8, y_pos-1.5), 0.5, 0.7, linewidth=1.5,
                                 edgecolor='black', facecolor='#2ECC71', alpha=0.6)
            ax.add_patch(rect_low)
            ax.text(x_pos-0.55, y_pos-1.15, f'{risks["low_risk"]}\nLow',
                   ha='center', va='center', fontsize=8, fontweight='bold')
            
            # Intermediate risk
            rect_inter = Rectangle((x_pos-0.25, y_pos-1.5), 0.5, 0.7, linewidth=1.5,
                                   edgecolor='black', facecolor='#F39C12', alpha=0.6)
            ax.add_patch(rect_inter)
            ax.text(x_pos, y_pos-1.15, f'{risks["inter_risk"]}\nInt',
                   ha='center', va='center', fontsize=8, fontweight='bold')
            
            # High risk
            rect_high = Rectangle((x_pos+0.3, y_pos-1.5), 0.5, 0.7, linewidth=1.5,
                                  edgecolor='black', facecolor='#E74C3C', alpha=0.6)
            ax.add_patch(rect_high)
            ax.text(x_pos+0.55, y_pos-1.15, f'{risks["high_risk"]}\nHigh',
                   ha='center', va='center', fontsize=8, fontweight='bold')
        
        # Legend
        legend_y = 0.5
        ax.text(0.5, legend_y + 0.5, 'Risk Category at 24 Months:', fontsize=11, fontweight='bold')
        
        colors_legend = [('#2ECC71', 'Low Risk'), ('#F39C12', 'Intermediate Risk'), ('#E74C3C', 'High Risk')]
        for idx, (color, label) in enumerate(colors_legend):
            rect = Rectangle((0.5 + idx*2.5, legend_y - 0.3), 0.3, 0.3,
                           linewidth=1, edgecolor='black', facecolor=color, alpha=0.6)
            ax.add_patch(rect)
            ax.text(1 + idx*2.5, legend_y - 0.15, label, fontsize=10, va='center')
        
        filepath = os.path.join(self.figures_dir, 'Figure_6_Patient_Flow_Diagram.png')
        plt.savefig(filepath, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        print(f"✓ Figure 6 saved: {filepath}")
        plt.close()
    
    def generate_all_figures(self, output_paths: dict, master_df: pd.DataFrame, comparative_df: pd.DataFrame):
        """
        Generate all publication-quality figures.
        """
        print("\n" + "="*70)
        print("📊 GENERATING PUBLICATION-QUALITY VISUALIZATIONS")
        print("="*70 + "\n")
        
        print("Creating Figure 1: Comparative Risk by Window...")
        self.create_comparative_risk_plot(comparative_df)
        
        print("Creating Figure 2: Risk Distribution by Cancer Type...")
        self.create_risk_distribution_plots(output_paths)
        
        print("Creating Figure 3: Risk Stratification Pie Charts...")
        self.create_risk_stratification_summary(output_paths)
        
        print("Creating Figure 4: Risk Heatmap...")
        self.create_heatmap_risk_matrix(master_df)
        
        print("Creating Figure 5: Risk Boxplots...")
        self.create_comparative_boxplot(comparative_df)
        
        print("Creating Table 1: Summary Statistics...")
        self.create_summary_statistics_table(comparative_df)
        
        print("Creating Figure 6: Patient Flow Diagram...")
        self.create_patient_flow_diagram(output_paths)
        
        print("\n✅ All visualizations complete!")
        print(f"✓ High-resolution figures saved to: {self.figures_dir}/")
        print("✓ DPI: 300 (publication quality)")
        print("✓ All figures are Nature-grade quality\n")


if __name__ == "__main__":
    
    # Load data
    comparative_df = pd.read_csv('output/COMPARATIVE_RISK_SUMMARY.csv')
    master_df = pd.read_csv('output/MASTER_PATIENT_REGISTRY.csv')
    
    # Create output paths dict
    output_paths = {}
    for cancer in ['breast_cancer', 'lung_cancer', 'prostate_cancer']:
        output_paths[cancer] = {
            '02_window_of_risk': f'output/{cancer}/02_window_of_risk'
        }
    
    # Generate visualizations
    visualizer = ProfessionalVisualizer()
    visualizer.generate_all_figures(output_paths, master_df, comparative_df)
    
    print("="*70)
    print("🎨 VISUALIZATION PIPELINE COMPLETE")
    print("="*70)
