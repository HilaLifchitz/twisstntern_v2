"""
Script to download twisst.py from GitHub.
"""

import os
import requests
from pathlib import Path

def download_twisst():
    """Download twisst.py from GitHub."""
    # GitHub raw content URL for twisst.py
    url = "https://raw.githubusercontent.com/simonhmartin/twisst/master/twisst.py"
    
    # Get the directory of this script
    script_dir = Path(__file__).parent
    
    # Download the file
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Save the file
        with open(script_dir / "twisst.py", "w") as f:
            f.write(response.text)
        
        print("Successfully downloaded twisst.py")
    except Exception as e:
        print(f"Error downloading twisst.py: {str(e)}")
        print("Please download twisst.py manually from:")
        print("https://github.com/simonhmartin/twisst/blob/master/twisst.py")
        print("and place it in the twisstntern/external directory.")

if __name__ == "__main__":
    download_twisst() 