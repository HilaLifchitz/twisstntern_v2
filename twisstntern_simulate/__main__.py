#!/usr/bin/env python
# coding: utf-8

import argparse
import sys
import os
from .main import main as twisstntern_simulate_main
from .utils import download_config_template


def main():
    """Main entry point for twisstntern_simulate command-line interface."""
    
    # Check for special commands first
    if len(sys.argv) > 1 and sys.argv[1] == '--get-config':
        print("🔧 TWISSTNTERN_SIMULATE Configuration Template")
        print("=" * 50)
        
        if len(sys.argv) > 2:
            # Custom destination path provided
            destination = sys.argv[2]
            downloaded_path = download_config_template(destination)
        else:
            # Default to current directory
            downloaded_path = download_config_template()
        
        if downloaded_path:
            print(f"📁 Template ready for use: {downloaded_path}")
            print("\n💡 Next steps:")
            print("   1. Edit the configuration file to match your simulation needs")
            print("   2. Run: python -m twisstntern_simulate -c config_template.yaml -o results")
        sys.exit(0)
    
    # Otherwise, run the normal main function
    twisstntern_simulate_main()


if __name__ == "__main__":
    main()
