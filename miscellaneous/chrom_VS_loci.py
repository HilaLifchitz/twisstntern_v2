#!/usr/bin/env python3
"""
FINAL CORRECTED Experiment: Find parameter combinations for ~10,000 trees
Fixes all plotting and table issues
"""

import msprime
import tskit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import time
import logging
from typing import List, Tuple, Dict
import itertools

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FinalCorrectedExperiment:
    """
    FINAL corrected experiment with proper plotting and complete tables.
    """
    
    def __init__(self, output_dir: str = "final_corrected_experiment"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.target_trees = 10000
        self.results = []
        
    def create_demography(self, Ne: float = 1000) -> msprime.Demography:
        """Create the CORRECT demography model matching config_template.yaml."""
        demography = msprime.Demography()
        
        # Add populations (matching config_template.yaml)
        demography.add_population(name="O", initial_size=Ne)      # Outgroup
        demography.add_population(name="p1", initial_size=Ne)     # Population 1
        demography.add_population(name="p2", initial_size=Ne)     # Population 2
        demography.add_population(name="p3", initial_size=Ne)     # Population 3
        demography.add_population(name="p12", initial_size=Ne)    # Ancestor of p1 and p2
        demography.add_population(name="p123", initial_size=Ne)   # Ancestor of p12 and p3
        demography.add_population(name="ANC", initial_size=Ne)    # Most ancestral population
        
        # Add population splits (matching config_template.yaml exactly)
        # Split 1: p1 and p2 split from p12 at time 100
        demography.add_population_split(time=100, derived=["p1", "p2"], ancestral="p12")
        # Split 2: p12 and p3 split from p123 at time 200
        demography.add_population_split(time=200, derived=["p12", "p3"], ancestral="p123")
        # Split 3: p123 and O split from ANC at time 300
        demography.add_population_split(time=300, derived=["p123", "O"], ancestral="ANC")
        
        return demography
    
    def get_targeted_combinations(self) -> List[Tuple[float, float, float]]:
        """
        Get parameter combinations that should actually reach the target lines.
        For 10,000 trees: rec_rate * chromosome_length ≈ 10 (since Ne=1000)
        """
        combinations = []
        
        # Define parameter ranges that will actually hit our target lines
        rec_rates = [1e-8,1e-7,1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3]
        chromosome_lengths = [1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6,1e7,1e8]
        
        # Target products: 1, 2, and 10 (for 10k trees)
        target_products = [1.0, 2.0, 10.0]
        
        for target_product in target_products:
            for rec_rate in rec_rates:
                for chrom_length in chromosome_lengths:
                    product = rec_rate * chrom_length
                    # Check if this combination is close to our target product
                    if abs(product - target_product) / target_product < 0.15:  # Within 15%
                        combinations.append((rec_rate, chrom_length, target_product))
        
        return combinations
    
    def run_single_simulation(self, rec_rate: float, chromosome_length: float, seed: int = None) -> Dict:
        """Run a single msprime simulation and return statistics."""
        demography = self.create_demography()
        samples = {"O": 10, "p1": 10, "p2": 10, "p3": 10}
        
        start_time = time.time()
        
        ts = msprime.sim_ancestry(
            samples=samples,
            demography=demography,
            sequence_length=chromosome_length,
            recombination_rate=rec_rate,
            ploidy=1,
            random_seed=seed
        )
        
        simulation_time = time.time() - start_time
        num_trees = ts.num_trees
        
        # Calculate how close we are to target
        distance_from_target = abs(num_trees - self.target_trees)
        relative_error = distance_from_target / self.target_trees
        
        return {
            'rec_rate': rec_rate,
            'chromosome_length': chromosome_length,
            'num_trees': num_trees,
            'simulation_time': simulation_time,
            'distance_from_target': distance_from_target,
            'relative_error': relative_error,
            'seed': seed
        }
    
    def run_experiment(self, replicates: int = 2) -> pd.DataFrame:
        """Run the experiment."""
        
        combinations = self.get_targeted_combinations()
        
        logger.info(f"Found {len(combinations)} parameter combinations to test")
        logger.info(f"Target: {self.target_trees} trees")
        
        results = []
        simulation_count = 0
        total_simulations = len(combinations) * replicates
        
        for rec_rate, chrom_length, target_product in combinations:
            for replicate in range(replicates):
                simulation_count += 1
                seed = 42 + simulation_count
                
                logger.info(f"Simulation {simulation_count}/{total_simulations}: "
                           f"rec_rate={rec_rate:.2e}, length={chrom_length:.2e}, "
                           f"product={rec_rate*chrom_length:.2e}, replicate={replicate+1}")
                
                try:
                    result = self.run_single_simulation(rec_rate, chrom_length, seed)
                    result['target_product'] = target_product
                    result['actual_product'] = rec_rate * chrom_length
                    results.append(result)
                    
                except Exception as e:
                    logger.error(f"Simulation failed: {e}")
                    continue
        
        df = pd.DataFrame(results)
        df.to_csv(self.output_dir / "final_results.csv", index=False)
        
        return df
    
    def create_summary_table(self, df: pd.DataFrame) -> pd.DataFrame:
        """Create a summary table with ALL tested combinations."""
        
        # Calculate average results for each parameter combination
        summary = df.groupby(['rec_rate', 'chromosome_length']).agg({
            'num_trees': ['mean', 'std'],
            'relative_error': 'mean',
            'simulation_time': 'mean',
            'actual_product': 'mean'
        }).round(2)
        
        # Flatten column names
        summary.columns = ['_'.join(col).strip() for col in summary.columns]
        summary = summary.reset_index()
        
        # Sort by actual product for logical ordering
        summary = summary.sort_values('actual_product_mean')
        
        # Add formatted columns
        summary['rec_rate_formatted'] = summary['rec_rate'].apply(lambda x: f"{x:.2e}")
        summary['chrom_length_formatted'] = summary['chromosome_length'].apply(lambda x: f"{x:.0e}")
        summary['product_formatted'] = summary['actual_product_mean'].apply(lambda x: f"{x:.2e}")
        
        return summary
    
    def visualize_results(self, df: pd.DataFrame, summary_table: pd.DataFrame) -> None:
        """Create corrected visualizations with proper target lines and complete tables."""
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f'FINAL Corrected Experiment: Finding ~{self.target_trees:,} Trees', fontsize=16)
        
        # 1. Scatter plot: rec_rate vs chromosome_length colored by tree count
        ax1 = axes[0, 0]
        scatter = ax1.scatter(df['rec_rate'], df['chromosome_length'], 
                            c=df['num_trees'], cmap='viridis', alpha=0.7, s=60)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel('Recombination Rate')
        ax1.set_ylabel('Chromosome Length')
        ax1.set_title('Tree Count by Parameter Combination')
        plt.colorbar(scatter, ax=ax1, label='Number of Trees')
        
        # Add target lines that will actually be reached
        x_range = np.logspace(np.log10(df['rec_rate'].min()), np.log10(df['rec_rate'].max()), 100)
        
        # Target lines for products 1, 2, and 10
        y_target_1 = 1.0 / x_range
        y_target_2 = 2.0 / x_range
        y_target_10 = 10.0 / x_range
        
        ax1.plot(x_range, y_target_1, 'b--', alpha=0.8, linewidth=2, label='Product = 1')
        ax1.plot(x_range, y_target_2, 'g--', alpha=0.8, linewidth=2, label='Product = 2')
        ax1.plot(x_range, y_target_10, 'r--', alpha=0.8, linewidth=2, 
                label=f'Product = 10 (~{self.target_trees} trees)')
        
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Tree count distribution
        ax2 = axes[0, 1]
        ax2.hist(df['num_trees'], bins=20, alpha=0.7, edgecolor='black')
        ax2.axvline(self.target_trees, color='red', linestyle='--', linewidth=2,
                   label=f'Target: {self.target_trees}')
        ax2.set_xlabel('Number of Trees')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Distribution of Tree Counts')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. COMPLETE table of ALL tested combinations
        ax3 = axes[1, 0]
        ax3.axis('tight')
        ax3.axis('off')
        
        # Show ALL combinations, not just top 12
        table_data = summary_table[['rec_rate_formatted', 'chrom_length_formatted', 
                                  'product_formatted', 'num_trees_mean', 'relative_error_mean']]
        table_data.columns = ['Rec Rate', 'Chrom Length', 'Product', 'Avg Trees', 'Error %']
        table_data['Error %'] = (table_data['Error %'] * 100).round(1)
        
        # Create table with all data
        table = ax3.table(cellText=table_data.values,
                         colLabels=table_data.columns,
                         cellLoc='center',
                         loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1.1, 1.0)
        ax3.set_title(f'ALL Tested Parameter Combinations ({len(table_data)} total)')
        
        # 4. Product vs Tree Count with actual data points
        ax4 = axes[1, 1]
        ax4.scatter(df['actual_product'], df['num_trees'], alpha=0.6, s=30)
        ax4.set_xlabel('Product (rec_rate × chromosome_length)')
        ax4.set_ylabel('Number of Trees')
        ax4.set_title('Product vs Tree Count')
        ax4.grid(True, alpha=0.3)
        
        # Add theoretical line
        x_theory = np.linspace(df['actual_product'].min(), df['actual_product'].max(), 100)
        y_theory = 1 + x_theory * 1000  # Ne = 1000
        ax4.plot(x_theory, y_theory, 'r--', alpha=0.7, label='Theory: 1 + product × Ne')
        
        # Add target lines
        ax4.axvline(1.0, color='blue', linestyle=':', alpha=0.7, label='Product = 1')
        ax4.axvline(2.0, color='green', linestyle=':', alpha=0.7, label='Product = 2')
        ax4.axvline(10.0, color='red', linestyle=':', alpha=0.7, label='Product = 10')
        
        ax4.legend()
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'final_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    def print_complete_results(self, summary_table: pd.DataFrame) -> None:
        """Print ALL tested combinations."""
        
        print("\n" + "="*100)
        print(f" COMPLETE RESULTS: ALL TESTED PARAMETER COMBINATIONS")
        print("="*100)
        
        print(f"\nAll {len(summary_table)} tested combinations:")
        print("-" * 100)
        print(f"{'Rank':<4} {'Rec Rate':<12} {'Chrom Length':<15} {'Product':<12} {'Avg Trees':<10} {'Error %':<8}")
        print("-" * 100)
        
        for i, (_, row) in enumerate(summary_table.iterrows(), 1):
            error_pct = row['relative_error_mean'] * 100
            product = row['actual_product_mean']
            print(f"{i:<4} {row['rec_rate']:<12.2e} {row['chromosome_length']:<15.0e} "
                  f"{product:<12.2e} {row['num_trees_mean']:<10.0f} {error_pct:<8.1f}")
        
        # Find best combination for 10k trees
        best_for_10k = summary_table.loc[summary_table['relative_error_mean'].idxmin()]
        
        print("\n" + "="*100)
        print("💡 BEST COMBINATION FOR ~10,000 TREES:")
        print("="*100)
        print(f"   rec_rate: {best_for_10k['rec_rate']:.2e}")
        print(f"   chromosome_length: {best_for_10k['chromosome_length']:.0e}")
        print(f"   product: {best_for_10k['actual_product_mean']:.2e}")
        print(f"   Expected trees: ~{best_for_10k['num_trees_mean']:.0f}")
        print(f"   Error from target: {best_for_10k['relative_error_mean']*100:.1f}%")
        
        print(f"\n For your config_template.yaml:")
        print(f"   rec_rate: {best_for_10k['rec_rate']:.2e}")
        print(f"   chromosome_length: {best_for_10k['chromosome_length']:.0e}")
        
        # Save complete results
        with open(self.output_dir / 'complete_results.txt', 'w') as f:
            f.write(f"Complete results for ~{self.target_trees} trees:\n\n")
            f.write("Best combination:\n")
            f.write(f"  rec_rate: {best_for_10k['rec_rate']:.2e}\n")
            f.write(f"  chromosome_length: {best_for_10k['chromosome_length']:.0e}\n")
            f.write(f"  product: {best_for_10k['actual_product_mean']:.2e}\n")
            f.write(f"  Expected trees: ~{best_for_10k['num_trees_mean']:.0f}\n")
            f.write(f"  Error from target: {best_for_10k['relative_error_mean']*100:.1f}%\n\n")
            
            f.write("All tested combinations:\n")
            for i, (_, row) in enumerate(summary_table.iterrows(), 1):
                error_pct = row['relative_error_mean'] * 100
                product = row['actual_product_mean']
                f.write(f"{i:2d}. rec_rate={row['rec_rate']:.2e}, "
                       f"chromosome_length={row['chromosome_length']:.0e}, "
                       f"product={product:.2e}, trees={row['num_trees_mean']:.0f}, error={error_pct:.1f}%\n")

def main():
    """Main function to run the final corrected experiment."""
    print(" FINAL Corrected Tree Experiment: Finding ~10,000 Trees")
    print("Fixes all plotting and table issues")
    print("="*70)
    
    experiment = FinalCorrectedExperiment()
    
    # Run experiment
    results_df = experiment.run_experiment(replicates=2)
    
    # Create complete summary table
    summary_table = experiment.create_summary_table(results_df)
    
    # Visualize results
    experiment.visualize_results(results_df, summary_table)
    
    # Print complete results
    experiment.print_complete_results(summary_table)
    
    print(f"\n✅ Experiment complete! Check the 'final_corrected_experiment' directory for results.")

if __name__ == "__main__":
    main()