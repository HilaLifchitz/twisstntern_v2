#!/usr/bin/env python3

import twisstntern
import matplotlib.pyplot as plt
import os
import glob
import sys

def list_csv_files():
    """List all CSV files in the current directory and subdirectories, excluding virtual environment, system directories, Results directory, and Jupyter checkpoints."""
    csv_files = []
    # Directories to exclude
    exclude_dirs = {
        '.venv', '__pycache__', '.git', 'venv', 'env', 'node_modules', 'Results',
        '.ipynb_checkpoints'  # Jupyter notebook checkpoints
    }
    
    # Get the absolute path of the current directory
    current_dir = os.path.abspath('.')
    
    for root, dirs, files in os.walk('.'):
        # Convert root to absolute path
        abs_root = os.path.abspath(root)
        
        # Skip if we're in a virtual environment or system directory
        if any(exclude_dir in abs_root.split(os.sep) for exclude_dir in exclude_dirs):
            continue
            
        # Skip excluded directories
        dirs[:] = [d for d in dirs if d not in exclude_dirs]
        
        for file in files:
            # Only include files that end with .csv (case insensitive)
            if file.lower().endswith('.csv'):
                # Get relative path from current directory
                rel_path = os.path.join(root, file)
                csv_files.append(rel_path)
                
    return sorted(csv_files)

def main():
    # Print instructions
    print("\n=== TWISSTNTERN Analysis Tool ===")
    print("\nAvailable CSV files:")
    csv_files = list_csv_files()
    
    if not csv_files:
        print("No CSV files found in the current directory or subdirectories.")
        return
        
    # Display numbered list of files
    for i, file in enumerate(csv_files, 1):
        print(f"{i}. {file}")
    
    # Get file selection
    while True:
        try:
            choice = input("\nEnter the number of the file you want to analyze (or 'q' to quit): ")
            if choice.lower() == 'q':
                return
            choice = int(choice)
            if 1 <= choice <= len(csv_files):
                input_file = csv_files[choice - 1]
                break
            print(f"Please enter a number between 1 and {len(csv_files)}")
        except ValueError:
            print("Please enter a valid number")
    
    # Get granularity
    print("\nSelect granularity:")
    print("1. Coarse (0.25)")
    print("2. Fine (0.1)")
    print("3. Superfine (0.05)")
    
    while True:
        try:
            choice = input("\nEnter your choice (1-3): ")
            choice = int(choice)
            if choice == 1:
                granularity_value = "coarse"
                break
            elif choice == 2:
                granularity_value = "fine"
                break
            elif choice == 3:
                granularity_value = "superfine"
                break
            print("Please enter 1, 2, or 3")
        except ValueError:
            print("Please enter a valid number")

    try:
        print(f"\nProcessing file: {input_file}")
        print("This may take a few moments...")
        
        # Process the data using twisstntern
        results = twisstntern.run_analysis(input_file, granularity_value)
        print(f"\nAnalysis complete! Results have been saved in the Results directory.")
        print("You can find:")
        print("1. CSV Results Table")
        print("2. Basic Asymmetry Plot")
        print("3. Ternary Coordinates Plot")
        print("4. Subtriangles Index Plot")
        print("5. Results Summary Plot")
        
    except Exception as e:
        print(f"\nError processing data: {str(e)}")
        print("Please make sure the file is a valid CSV with three columns (T1, T2, T3)")
        # Print more detailed error information
        import traceback
        print("\nDetailed error information:")
        traceback.print_exc()

if __name__ == "__main__":
    main() 