"""
Test script for the simulation functionality.
"""

import os
import sys
import unittest
import msprime
import numpy as np
from pathlib import Path

# Add parent directory to path to import twisstntern
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import twisstntern

class TestSimulation(unittest.TestCase):
    """Test cases for simulation functionality."""
    
    def setUp(self):
        """Set up test environment."""
        # Create test directory
        self.test_dir = Path("test_simulation_results")
        self.test_dir.mkdir(exist_ok=True)
        
        # Create test config file
        self.config_file = self.test_dir / "test_config.yaml"
        with open(self.config_file, "w") as f:
            f.write("""# Test configuration
# Population sizes
ne_p1 = 1000
ne_p2 = 1000
ne_p3 = 1000
ne_p12 = 1500
ne_p123 = 2000
ne_O = 2500

# Split times
t1 = 1000
t2 = 2000
t3 = 3000

# Simulation parameters
n_ind = 5
n_loci = 10
seq_length = 1000
rec_rate = 1e-8
mutation_rate = 1e-8
seed = 42
""")
    
    def tearDown(self):
        """Clean up after tests."""
        # Remove test files
        if self.config_file.exists():
            self.config_file.unlink()
        # Remove test directory
        if self.test_dir.exists():
            for file in self.test_dir.glob("*"):
                file.unlink()
            self.test_dir.rmdir()
    
    def test_locus_simulation(self):
        """Test locus simulation."""
        # Run simulation
        tree_sequences = twisstntern.simulate_locus(str(self.config_file))
        
        # Check results
        self.assertEqual(len(tree_sequences), 10)  # Should have 10 loci
        
        # Check each tree sequence
        for ts in tree_sequences:
            self.assertIsInstance(ts, msprime.TreeSequence)
            self.assertEqual(ts.num_samples, 20)  # 5 individuals * 4 populations
            self.assertEqual(ts.sequence_length, 1000)  # From config
    
    def test_chromosome_simulation(self):
        """Test chromosome simulation."""
        # Run simulation
        ts = twisstntern.simulate_chromosome(str(self.config_file))
        
        # Check results
        self.assertIsInstance(ts, msprime.TreeSequence)
        self.assertEqual(ts.num_samples, 20)  # 5 individuals * 4 populations
        self.assertEqual(ts.sequence_length, 1000)  # From config
    
    def test_save_tree_sequences(self):
        """Test saving tree sequences."""
        # Simulate and save loci
        tree_sequences = twisstntern.simulate_locus(str(self.config_file))
        output_dir = self.test_dir / "locus_sim"
        twisstntern.save_tree_sequences(tree_sequences, str(output_dir))
        
        # Check saved files
        self.assertTrue(output_dir.exists())
        self.assertEqual(len(list(output_dir.glob("*.ts"))), 10)
        
        # Simulate and save chromosome
        ts = twisstntern.simulate_chromosome(str(self.config_file))
        output_dir = self.test_dir / "chrom_sim"
        twisstntern.save_tree_sequences(ts, str(output_dir))
        
        # Check saved files
        self.assertTrue(output_dir.exists())
        self.assertEqual(len(list(output_dir.glob("*.ts"))), 1)
    
    def test_optional_parameters(self):
        """Test handling of optional parameters."""
        # Create config without mutation rate and seed
        with open(self.config_file, "w") as f:
            f.write("""# Test configuration without optional parameters
# Population sizes
ne_p1 = 1000
ne_p2 = 1000
ne_p3 = 1000
ne_p12 = 1500
ne_p123 = 2000
ne_O = 2500

# Split times
t1 = 1000
t2 = 2000
t3 = 3000

# Simulation parameters
n_ind = 5
n_loci = 10
seq_length = 1000
rec_rate = 1e-8
""")
        
        # Test locus simulation
        tree_sequences = twisstntern.simulate_locus(str(self.config_file))
        self.assertEqual(len(tree_sequences), 10)
        
        # Test chromosome simulation
        ts = twisstntern.simulate_chromosome(str(self.config_file))
        self.assertIsInstance(ts, msprime.TreeSequence)

if __name__ == "__main__":
    unittest.main() 