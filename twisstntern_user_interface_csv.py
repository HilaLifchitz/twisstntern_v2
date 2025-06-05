#!/usr/bin/env python
# coding: utf-8

import os
from pathlib import Path
import twisstntern
import glob

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
    print("\nChoose granularity level:")
    print("1. superfine (0.05)")
    print("2. fine (0.1)")
    print("3. coarse (0.25)")
    print("4. custom -- granuality alpha should be chosen such that 1/alpha is an even number")
    
    choice = input("\nEnter your choice (1-4): ")
    
    if choice == "1":
        granularity = "superfine"
    elif choice == "2":
        granularity = "fine"
    elif choice == "3":
        granularity = "coarse"
    elif choice == "4":
        try:
            value = float(input("Enter custom granularity value (must be even divisor of 1): "))
            if 1/value % 2 != 0:
                print("Error: 1/granularity must be even")
                return
            granularity = value
        except ValueError:
            print("Error: Invalid granularity value")
            return
    else:
        print("Error: Invalid choice")
        return

    try:
        print(f"\nProcessing file: {input_file}")
        print("This may take a few moments...")
        
        # Run analysis
        results, fundamental_results = twisstntern.run_analysis(input_file, granularity)
        
        # Print results
        print("\nAnalysis complete!")
        print("\nFundamental asymmetry results:")
        print(f"n_right: {fundamental_results[0]}")
        print(f"n_left: {fundamental_results[1]}")
        print(f"D_LR: {fundamental_results[2]:.4f}")
        print(f"G-test: {fundamental_results[3]:.4f}")
        print(f"p-value: {fundamental_results[4]:.4e}")
        
        print("\nResults have been saved to the Results directory")
        print("You can find:")
        print("1. CSV Results Table")
        print("2. Basic Asymmetry Plot")
        print("3. Ternary Coordinates Plot")
        print("4. Subtriangles Index Plot")
        print("5. Results Summary Plot")
        
    except Exception as e:
        print(f"\nError: {str(e)}")
        print("Please make sure the file is a valid CSV with three columns (T1, T2, T3)")
        # Print more detailed error information
        import traceback
        print("\nDetailed error information:")
        traceback.print_exc()

if __name__ == "__main__":
    main() 