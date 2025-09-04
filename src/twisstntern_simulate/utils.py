#!/usr/bin/env python
# coding: utf-8

"""
Package-specific utilities for twisstntern_simulate.

This module contains utilities that are specific to the twisstntern_simulate package
and are not shared with twisstntern. All shared utilities have been moved to the core module.
"""

def download_config_template(destination_path=None):
    """
    Download the config_template.yaml file from GitHub.
    
    This function downloads the latest config template directly from the
    GitHub repository, ensuring users always get the most up-to-date version.
    
    Args:
        destination_path (str, optional): Where to save the config file.
                                        Defaults to current directory.
    
    Returns:
        str: Path where the config file was saved
        
    Examples:
        >>> import twisstntern_simulate.utils as utils
        >>> config_path = utils.download_config_template()
        >>> print(f"Config template downloaded to: {config_path}")
    """
    import requests
    import os
    
    if destination_path is None:
        destination_path = os.path.join(os.getcwd(), 'config_template.yaml')
    
    # GitHub raw URL for the config template
    url = "https://raw.githubusercontent.com/HilaLifchitz/twisstntern_v2/main/config_template.yaml"
    
    try:
        print("📥 Downloading config template from GitHub...")
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        
        with open(destination_path, 'w', encoding='utf-8') as f:
            f.write(response.text)
        
        print(f"✅ Config template downloaded to: {destination_path}")
        return destination_path
        
    except requests.exceptions.RequestException as e:
        print(f"❌ Failed to download config template: {e}")
        print("🔗 Please download manually from:")
        print("   https://github.com/HilaLifchitz/twisstntern_v2/raw/main/config_template.yaml")
        return None
    except Exception as e:
        print(f"❌ Error saving config template: {e}")
        return None