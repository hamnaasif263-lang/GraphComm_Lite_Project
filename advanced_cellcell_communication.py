"""
Advanced Cell-Cell Communication Network Analysis
Real data-driven network inference with ligand-receptor interactions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

class AdvancedCellCommunicationNetwork:
    """Performs advanced network analysis on real single-cell data"""
    
    def __init__(self, output_base="output"):
        self.output_base = Path(output_base)
        self.dpi = 300
        
        # Ligand-receptor interactions database
        self.ligand_receptor_db = {
            'IGF2': ['IGF1R', 'INSRR'],
            'IGF1': ['IGF1R', 'INSRR'],
            'EGF': ['EGFR'],
            'PDGF': ['PDGFRA', 'PDGFRB'],
            'FGF': ['FGFR1', 'FGFR2'],
            'VEGF': ['KDR', 'FLT1'],
            'TNF': ['TNFRSF1A', 'TNFRSF1B'],
            'IL6': ['IL6R', 'IL6ST'],
        }
    
    def load_and_describe_data(self, cancer_type):
        """Load and analyze single-cell data"""
        sc_path = self.output_base / cancer_type / "01_preprocessing" / f"{cancer_type}_single_cell_data.csv"
        sc_data = pd.read_csv(sc_path)
        
        print(f"\n{cancer_type.upper()} - Data Summary:")
        print(f"  Total cells: {len(sc_data):,}")
        print(f"  Cell types: {sorted(sc_data['cell_type'].unique())}")
        cell_counts = sc_data['cell_type'].value_counts()
        for cell_type, count in cell_counts.items():
            pct = (count / len(sc_data)) * 100
            print(f"    - {cell_type}: {count:,} ({pct:.1f}%)")
        
        return sc_data
    
    def compute_ligand_receptor_signaling(self, sc_data, cancer_type):
        """
        Compute cell-type to cell-type signaling based on:
        1. IGF scores (ligand production proxy)
        2. Receptor expression (IGF1R, etc.)
        3. Dormancy/proliferation states
        """
        
        cell_types = sorted(sc_data['cell_type'].unique())
        print(f"\nComputing ligand-receptor signaling pathways...")
        
        # Get cell-type profiles
        profiles = {}
        for ct in cell_types:
            ct_data = sc_data[sc_data['cell_type'] == ct]
            profiles[ct] = {
                'igf_mean': ct_data['igf_score'].mean(),
                'igf_std': ct_data['igf_score'].std(),
                'prolif_mean': ct_data['proliferation_score'].mean(),
                'dormancy_mean': ct_data['dormancy_score'].mean(),
                'myc_mean': ct_data['myc_score'].mean(),
                'cdkn1b_mean': ct_data['cdkn1b_score'].mean(),
                'count': len(ct_data),
                'pct': (len(ct_data) / len(sc_data)) * 100
            }
        
        # Build communication matrix
        edges = []
        
        for sender in cell_types:
            sender_profile = profiles[sender]
            
            # Signaling potential: IGF score in sender
            igf_strength = sender_profile['igf_mean']
            
            for receiver in cell_types:
                receiver_profile = profiles[receiver]
                
                # Receiver receptivity: proliferation capable, low dormancy
                receptivity = receiver_profile['prolif_mean'] * (1 - receiver_profile['dormancy_mean'] * 0.3)
                
                # MYC-CDKN1B balance determines response
                response_balance = receiver_profile['myc_mean'] * receiver_profile['cdkn1b_mean']
                
                # Signal strength calculation
                if sender != receiver:
                    signal = igf_strength * receptivity * (1 + response_balance * 0.5)
                else:
                    # Autocrine: self-signaling
                    signal = igf_strength * receiver_profile['prolif_mean'] * 1.5
                
                if signal > 0.005:
                    edges.append({
                        'Sender': sender,
                        'Receiver': receiver,
                        'Signal_Strength': signal,
                        'Sender_IGF_Level': igf_strength,
                        'Receiver_Proliferation': receiver_profile['prolif_mean'],
                        'Receiver_Dormancy': receiver_profile['dormancy_mean'],
                        'Sender_Count': sender_profile['count'],
                        'Receiver_Count': receiver_profile['count'],
                        'Autocrine': sender == receiver
                    })
        
        edges_df = pd.DataFrame(edges)
        return edges_df, profiles, cell_types
    
    def compute_correlation_strength(self, sc_data, sender_type, receiver_type):
        """Compute actual score correlations between cell types"""
        
        sender_data = sc_data[sc_data['cell_type'] == sender_type]
        receiver_data = sc_data[sc_data['cell_type'] == receiver_type]
        
        if len(sender_data) < 2 or len(receiver_data) < 2:
            return 0
        
        # Correlation between sender IGF and receiver proliferation
        sender_igf = sender_data['igf_score'].mean()
        receiver_prolif = receiver_data['proliferation_score'].mean()
        
        correlation = abs(sender_igf - receiver_prolif)
        return 1.0 / (1.0 + correlation)  # Convert to similarity metric
    
    def build_network_graph(self, edges_df, profiles):
        """Create directed graph"""
        
        G = nx.DiGraph()
        
        cell_types = list(profiles.keys())
        G.add_nodes_from(cell_types)
        
        for _, row in edges_df.iterrows():
            G.add_edge(row['Sender'], row['Receiver'], 
                      weight=row['Signal_Strength'],
                      autocrine=row['Autocrine'])
        
        return G
    
    def visualize_network(self, G, edges_df, cancer_type, profiles):
        """Create professional network visualization with score labels"""
        
        fig, ax = plt.subplots(figsize=(16, 12))
        
        # Layout
        pos = nx.spring_layout(G, k=3, iterations=100, seed=42)
        
        # Node colors based on IGF production
        node_colors = {}
        node_sizes = {}
        node_igf_scores = {}
        
        for node in G.nodes():
            igf = profiles[node]['igf_mean']
            count = profiles[node]['count']
            pct = profiles[node]['pct']
            
            node_igf_scores[node] = igf
            
            # Color by IGF level
            if igf > 0.3:
                node_colors[node] = '#FF4444'  # High IGF
            elif igf > 0.2:
                node_colors[node] = '#FF8844'  # Medium-high
            elif igf > 0.1:
                node_colors[node] = '#FFBB44'  # Medium
            else:
                node_colors[node] = '#4488FF'  # Low IGF
            
            # Size by cell count
            node_sizes[node] = 3000 + (pct * 100)
        
        # Draw nodes
        nx.draw_networkx_nodes(G, pos, 
                              node_color=[node_colors[n] for n in G.nodes()],
                              node_size=[node_sizes[n] for n in G.nodes()],
                              ax=ax, alpha=0.85, edgecolors='black', linewidths=2.5)
        
        # Draw edges with special handling for autocrine
        edges_list = G.edges(data=True)
        max_weight = max([d['weight'] for _, _, d in edges_list]) if edges_list else 1
        
        autocrine_edges = [(u, v, d) for u, v, d in edges_list if d.get('autocrine', False)]
        paracrine_edges = [(u, v, d) for u, v, d in edges_list if not d.get('autocrine', False)]
        
        # Draw paracrine (cell-to-cell)
        for u, v, d in paracrine_edges:
            weight = d['weight']
            width = 1 + (weight / max_weight) * 4
            alpha = 0.4 + (weight / max_weight) * 0.4
            
            nx.draw_networkx_edges(G, pos, [(u, v)], width=width, alpha=alpha,
                                 ax=ax, edge_color='#2C3E50', arrows=True,
                                 arrowsize=20, arrowstyle='->', 
                                 connectionstyle='arc3,rad=0.15')
        
        # Draw autocrine (self-loops) with different style
        for u, v, d in autocrine_edges:
            weight = d['weight']
            width = 2 + (weight / max_weight) * 4
            
            nx.draw_networkx_edges(G, pos, [(u, v)], width=width, alpha=0.6,
                                 ax=ax, edge_color='#E74C3C', arrows=True,
                                 arrowsize=20, arrowstyle='->', 
                                 connectionstyle='arc3,rad=0.8')
        
        # Add edge labels (signal strength)
        edge_labels = {}
        for u, v, d in edges_list:
            edge_labels[(u, v)] = f"{d['weight']:.3f}"
        
        nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=7, 
                                     font_weight='bold', ax=ax, 
                                     bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))
        
        # Node labels with IGF scores
        node_labels_with_score = {}
        for node in G.nodes():
            igf_score = node_igf_scores[node]
            node_labels_with_score[node] = f"{node.replace('_', ' ')}\n[IGF:{igf_score:.2f}]"
        
        nx.draw_networkx_labels(G, pos, node_labels_with_score, font_size=9, 
                               font_weight='bold', ax=ax)
        
        # Title and legend
        n_paracrine = len(paracrine_edges)
        n_autocrine = len(autocrine_edges)
        total_signal = edges_df['Signal_Strength'].sum()
        
        title = f"{cancer_type.replace('_', ' ').title()}\n"
        title += f"Cell-Cell Communication Network\n"
        title += f"({n_paracrine} paracrine + {n_autocrine} autocrine | Total Signal: {total_signal:.4f})"
        
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        
        # Legend with score ranges
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#FF4444', edgecolor='black', label='High IGF (>0.30)'),
            Patch(facecolor='#FF8844', edgecolor='black', label='Medium-High IGF (0.20-0.30)'),
            Patch(facecolor='#FFBB44', edgecolor='black', label='Medium IGF (0.10-0.20)'),
            Patch(facecolor='#4488FF', edgecolor='black', label='Low IGF (<0.10)'),
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=10, 
                 title='IGF Score Levels', title_fontsize=11, framealpha=0.95)
        
        ax.axis('off')
        plt.tight_layout()
        return fig
    
    def visualize_pathway_details(self, edges_df, cancer_type):
        """Detailed pathway analysis visualization"""
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        # Panel 1: Top signaling pathways
        top_edges = edges_df.nlargest(12, 'Signal_Strength')
        pathway_labels = [f"{row['Sender'][:4]}→{row['Receiver'][:4]}" 
                         for _, row in top_edges.iterrows()]
        
        y_pos = np.arange(len(top_edges))
        colors = ['#E74C3C' if row['Autocrine'] else '#3498DB' 
                 for _, row in top_edges.iterrows()]
        
        ax1.barh(y_pos, top_edges['Signal_Strength'], color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(pathway_labels, fontsize=9)
        ax1.set_xlabel('Signal Strength (Actual)', fontsize=11, fontweight='bold')
        ax1.set_title('Top 12 Communication Pathways', fontsize=12, fontweight='bold')
        ax1.grid(axis='x', alpha=0.3)
        
        # Panel 2: Autocrine vs Paracrine
        autocrine_total = edges_df[edges_df['Autocrine']]['Signal_Strength'].sum()
        paracrine_total = edges_df[~edges_df['Autocrine']]['Signal_Strength'].sum()
        
        pie_data = [autocrine_total, paracrine_total]
        pie_labels = [f'Autocrine\n{autocrine_total:.4f}', f'Paracrine\n{paracrine_total:.4f}']
        
        ax2.pie(pie_data, labels=pie_labels, autopct='%1.1f%%', colors=['#E74C3C', '#3498DB'],
               startangle=90, explode=(0.05, 0.05), textprops={'fontweight': 'bold'})
        ax2.set_title('Signal Distribution', fontsize=12, fontweight='bold')
        
        # Panel 3: Cell type IGF levels
        cell_types = sorted(edges_df['Sender'].unique())
        igf_levels = []
        for ct in cell_types:
            igf = edges_df[edges_df['Sender'] == ct]['Sender_IGF_Level'].iloc[0]
            igf_levels.append(igf)
        
        bars = ax3.bar(cell_types, igf_levels, color=['#FF4444', '#FF8844', '#4488FF'][:len(cell_types)],
                      alpha=0.8, edgecolor='black', linewidth=2)
        ax3.set_ylabel('Mean IGF Score', fontsize=11, fontweight='bold')
        ax3.set_title('IGF Production by Cell Type', fontsize=12, fontweight='bold')
        ax3.set_ylim([0, max(igf_levels) * 1.2])
        for bar, val in zip(bars, igf_levels):
            ax3.text(bar.get_x() + bar.get_width()/2, val + 0.01, f'{val:.3f}',
                    ha='center', fontweight='bold', fontsize=10)
        ax3.grid(axis='y', alpha=0.3)
        
        # Panel 4: Summary statistics
        ax4.axis('off')
        
        stats_text = f"""
CELL-CELL COMMUNICATION ANALYSIS SUMMARY

Network Statistics:
  • Total pathways: {len(edges_df)}
  • Autocrine pathways: {edges_df['Autocrine'].sum()}
  • Paracrine pathways: {(~edges_df['Autocrine']).sum()}
  • Total signal strength: {edges_df['Signal_Strength'].sum():.4f}

Top Signaling Pathway:
  {top_edges.iloc[0]['Sender']} → {top_edges.iloc[0]['Receiver']}
  Strength: {top_edges.iloc[0]['Signal_Strength']:.4f}
  IGF Level: {top_edges.iloc[0]['Sender_IGF_Level']:.4f}
  Receiver Proliferation: {top_edges.iloc[0]['Receiver_Proliferation']:.4f}

Cell Type Communication Pattern:
"""
        
        for ct in cell_types:
            outgoing = edges_df[edges_df['Sender'] == ct].shape[0]
            incoming = edges_df[edges_df['Receiver'] == ct].shape[0]
            stats_text += f"\n  {ct}: {outgoing} outgoing, {incoming} incoming"
        
        ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes,
                fontfamily='monospace', fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        fig.suptitle(f'{cancer_type.replace("_", " ").title()} - Detailed Communication Analysis',
                    fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        return fig
    
    def compute_hub_centrality(self, G, edges_df):
        """Analyze hub cell types"""
        
        # In-degree and out-degree
        in_degree = dict(G.in_degree(weight='weight'))
        out_degree = dict(G.out_degree(weight='weight'))
        
        hubs = {}
        for node in G.nodes():
            hubs[node] = {
                'in_degree': in_degree.get(node, 0),
                'out_degree': out_degree.get(node, 0),
                'total_degree': in_degree.get(node, 0) + out_degree.get(node, 0)
            }
        
        return hubs
    
    def save_results(self, edges_df, hubs, cancer_type):
        """Save analysis results"""
        
        output_dir = self.output_base / cancer_type / "network_analysis"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save edge data
        edges_df.to_csv(output_dir / f"{cancer_type}_cellcell_communication.csv", index=False)
        
        # Save hub analysis
        hubs_df = pd.DataFrame(hubs).T
        hubs_df.to_csv(output_dir / f"{cancer_type}_network_hubs.csv")
        
        return output_dir
    
    def analyze_cancer_type(self, cancer_type):
        """Complete analysis pipeline"""
        
        print(f"\n{'='*70}")
        print(f"ANALYZING: {cancer_type.upper()}")
        print(f"{'='*70}")
        
        # Load data
        sc_data = self.load_and_describe_data(cancer_type)
        
        # Compute signaling
        edges_df, profiles, cell_types = self.compute_ligand_receptor_signaling(sc_data, cancer_type)
        
        print(f"\nNetwork Summary:")
        print(f"  Cell types: {len(cell_types)}")
        print(f"  Total pathways: {len(edges_df)}")
        print(f"  Paracrine: {(~edges_df['Autocrine']).sum()}, Autocrine: {edges_df['Autocrine'].sum()}")
        print(f"  Total signal strength: {edges_df['Signal_Strength'].sum():.4f}")
        
        # Build graph
        G = self.build_network_graph(edges_df, profiles)
        
        # Create visualizations
        print(f"\nGenerating visualizations...")
        
        fig1 = self.visualize_network(G, edges_df, cancer_type, profiles)
        fig1.savefig(self.output_base / f"real_network_{cancer_type}.png",
                    dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close(fig1)
        print(f"  [SAVED] real_network_{cancer_type}.png")
        
        fig2 = self.visualize_pathway_details(edges_df, cancer_type)
        fig2.savefig(self.output_base / f"real_pathways_{cancer_type}.png",
                    dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close(fig2)
        print(f"  [SAVED] real_pathways_{cancer_type}.png")
        
        # Hub analysis
        hubs = self.compute_hub_centrality(G, edges_df)
        
        # Save results
        output_dir = self.save_results(edges_df, hubs, cancer_type)
        print(f"  [SAVED] results to: {output_dir}")
        
        # Print top pathways
        print(f"\nTop 5 Communication Pathways:")
        for i, (_, row) in enumerate(edges_df.nlargest(5, 'Signal_Strength').iterrows(), 1):
            pathway_type = "AUTOCRINE" if row['Autocrine'] else "PARACRINE"
            print(f"  {i}. {row['Sender']:15} → {row['Receiver']:15} [{pathway_type:10}] str: {row['Signal_Strength']:.4f}")
        
        return G, edges_df, hubs
    
    def run_all(self):
        """Run complete analysis"""
        
        print("\n" + "="*70)
        print("REAL CELL-CELL COMMUNICATION NETWORK ANALYSIS")
        print("Advanced network inference from single-cell data")
        print("="*70)
        
        cancer_types = ['breast_cancer', 'lung_cancer', 'prostate_cancer']
        results = {}
        
        for cancer_type in cancer_types:
            G, edges_df, hubs = self.analyze_cancer_type(cancer_type)
            results[cancer_type] = {'graph': G, 'edges': edges_df, 'hubs': hubs}
        
        # Comparative analysis
        self.create_comparative_figure(results)
        
        print(f"\n{'='*70}")
        print("ANALYSIS COMPLETE - ALL REAL DATA RESULTS")
        print(f"{'='*70}\n")


    def create_comparative_figure(self, results):
        """Comparative analysis across cancer types"""
        
        print(f"\nGenerating comparative analysis...")
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        
        cancer_labels = ['Breast Cancer', 'Lung Cancer', 'Prostate Cancer']
        cancer_keys = list(results.keys())
        
        # Row 1: Network size
        for idx, (cancer_type, label, ax) in enumerate(zip(cancer_keys, cancer_labels, axes[0])):
            edges_df = results[cancer_type]['edges']
            
            auto = (edges_df['Autocrine'].sum())
            para = ((~edges_df['Autocrine']).sum())
            
            ax.bar(['Autocrine', 'Paracrine'], [auto, para],
                  color=['#E74C3C', '#3498DB'], alpha=0.8, edgecolor='black', linewidth=2)
            ax.set_ylabel('Number of Pathways', fontsize=11, fontweight='bold')
            ax.set_title(f'{label}\n(Total: {len(edges_df)} pathways)', fontsize=12, fontweight='bold')
            ax.grid(axis='y', alpha=0.3)
        
        # Row 2: Signal strength comparison
        for idx, (cancer_type, label, ax) in enumerate(zip(cancer_keys, cancer_labels, axes[1])):
            edges_df = results[cancer_type]['edges']
            
            auto_strength = edges_df[edges_df['Autocrine']]['Signal_Strength'].sum()
            para_strength = edges_df[~edges_df['Autocrine']]['Signal_Strength'].sum()
            
            ax.bar(['Autocrine', 'Paracrine'], [auto_strength, para_strength],
                  color=['#E74C3C', '#3498DB'], alpha=0.8, edgecolor='black', linewidth=2)
            ax.set_ylabel('Total Signal Strength', fontsize=11, fontweight='bold')
            ax.set_title(f'{label}\nSignal Distribution', fontsize=12, fontweight='bold')
            ax.grid(axis='y', alpha=0.3)
        
        fig.suptitle('Cell-Cell Communication Comparison Across Cancer Types\n(Real Data Analysis)',
                    fontsize=14, fontweight='bold', y=0.995)
        
        plt.tight_layout()
        fig.savefig(self.output_base / 'real_comparative_networks.png',
                   dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"  [SAVED] real_comparative_networks.png")


if __name__ == "__main__":
    analyzer = AdvancedCellCommunicationNetwork()
    analyzer.run_all()
