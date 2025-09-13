#!/usr/bin/env python
# coding: utf-8

import unittest
import sys
from pathlib import Path

# Add src to Python path
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))


class TestImports(unittest.TestCase):
    """Test that all modules can be imported correctly."""
    
    def test_import_main_module(self):
        """Test importing the main twisstntern module."""
        try:
            import twisstntern
            self.assertTrue(hasattr(twisstntern, '__version__'))
        except ImportError as e:
            self.fail(f"Failed to import twisstntern: {e}")
    
    def test_import_utils(self):
        """Test importing utils module."""
        try:
            from twisstntern import utils
            self.assertTrue(hasattr(utils, 'cartizian'))
            self.assertTrue(hasattr(utils, 'dump_data'))
        except ImportError as e:
            self.fail(f"Failed to import utils: {e}")
    
    def test_import_analysis(self):
        """Test importing analysis module."""
        try:
            from twisstntern import analysis
            self.assertTrue(hasattr(analysis, 'fundamental_asymmetry'))
            self.assertTrue(hasattr(analysis, 'triangles_analysis'))
        except ImportError as e:
            self.fail(f"Failed to import analysis: {e}")
    
    def test_import_visualization(self):
        """Test importing visualization module."""
        try:
            from twisstntern import visualization
            self.assertTrue(hasattr(visualization, 'plot'))
            self.assertTrue(hasattr(visualization, 'plot_fundamental_asymmetry'))
        except ImportError as e:
            self.fail(f"Failed to import visualization: {e}")
    
    def test_import_config(self):
        """Test importing config module.""" 
        try:
            from twisstntern import config
            self.assertTrue(hasattr(config, 'TwisstnternConfig'))
        except ImportError as e:
            self.fail(f"Failed to import config: {e}")
    
    def test_import_twisstntern_simulate(self):
        """Test importing the refactored twisstntern_simulate module."""
        try:
            from src.twisstntern_simulate import Config, TwisstnternSimulateConfig
            self.assertTrue(Config is not None)
            self.assertTrue(TwisstnternSimulateConfig is not None)
        except ImportError as e:
            self.fail(f"Failed to import twisstntern_simulate: {e}")
    
    def test_hydra_config_structure(self):
        """Test that Hydra config structure is valid."""
        try:
            from src.twisstntern_simulate.hydra_config import TwisstnternSimulateConfig
            from omegaconf import OmegaConf
            
            # Test that we can create a structured config
            cfg = OmegaConf.structured(TwisstnternSimulateConfig)
            self.assertTrue(cfg is not None)
            self.assertTrue(hasattr(cfg, 'simulation'))
            self.assertTrue(hasattr(cfg, 'analysis'))
            self.assertTrue(hasattr(cfg, 'visualization'))
        except Exception as e:
            self.fail(f"Failed to validate Hydra config structure: {e}")


if __name__ == '__main__':
    unittest.main()