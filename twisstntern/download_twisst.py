"""
Script to download twisst.py from GitHub and fix NumPy 2.0 compatibility.
"""

import os
import requests
from pathlib import Path


def fix_numpy_compatibility(content):
    """Fix NumPy 2.0 compatibility issues in twisst.py"""
    # Replace np.NaN with np.nan (NumPy 2.0 compatibility)
    content = content.replace("np.NaN", "np.nan")

    print("✓ Fixed NumPy 2.0 compatibility (np.NaN → np.nan)")
    return content


def download_twisst():
    """Download twisst.py from GitHub and apply fixes."""
    # GitHub raw content URL for twisst.py
    url = "https://raw.githubusercontent.com/simonhmartin/twisst/master/twisst.py"

    # Get the directory of this script and create external subdirectory
    script_dir = Path(__file__).parent
    external_dir = script_dir / "external"
    external_dir.mkdir(exist_ok=True)  # Create external directory if it doesn't exist
    twisst_path = external_dir / "twisst.py"

    # Download the file
    try:
        print("Downloading twisst.py from GitHub...")
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes

        # Apply fixes to the content
        print("Applying compatibility fixes...")
        content = fix_numpy_compatibility(response.text)

        # Save the fixed file
        with open(twisst_path, "w") as f:
            f.write(content)

        print(f"✓ Successfully downloaded and fixed twisst.py to {twisst_path}")

    except Exception as e:
        print(f"✗ Error downloading twisst.py: {str(e)}")
        print("Please download twisst.py manually from:")
        print("https://github.com/simonhmartin/twisst/blob/master/twisst.py")
        print("and place it in the twisstntern/external directory.")


def fix_existing_twisst():
    """Fix an existing twisst.py file for NumPy 2.0 compatibility."""
    script_dir = Path(__file__).parent
    external_dir = script_dir / "external"
    twisst_path = external_dir / "twisst.py"

    if not twisst_path.exists():
        print(f"✗ twisst.py not found at {twisst_path}")
        return False

    try:
        # Read the existing file
        with open(twisst_path, "r") as f:
            content = f.read()

        # Check if it needs fixing
        if "np.NaN" not in content:
            print("✓ twisst.py already compatible with NumPy 2.0")
            return True

        # Apply fixes
        print("Fixing existing twisst.py for NumPy 2.0 compatibility...")
        fixed_content = fix_numpy_compatibility(content)

        # Write back the fixed content
        with open(twisst_path, "w") as f:
            f.write(fixed_content)

        print("✓ Successfully fixed existing twisst.py")
        return True

    except Exception as e:
        print(f"✗ Error fixing twisst.py: {str(e)}")
        return False


def ensure_twisst_available():
    """Ensure twisst.py is available and working. Download if necessary."""
    script_dir = Path(__file__).parent
    external_dir = script_dir / "external"
    twisst_path = external_dir / "twisst.py"

    if twisst_path.exists():
        print("Found existing twisst.py")
        if not fix_existing_twisst():
            print("Trying to re-download...")
            download_twisst()
    else:
        print("twisst.py not found, downloading...")
        download_twisst()

    # Verify that twisst is now available
    if twisst_path.exists():
        print(f"✓ twisst.py is available at {twisst_path}")
        return True
    else:
        print("✗ Failed to make twisst.py available")
        return False


if __name__ == "__main__":
    ensure_twisst_available()
