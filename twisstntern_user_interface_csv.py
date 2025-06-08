#!/usr/bin/env python
# coding: utf-8

import os
from pathlib import Path
import twisstntern
import glob

def list_supported_files():
    """List all supported files (CSV and tree files) in the current directory and subdirectories, 
    excluding virtual environment, system directories, Results directory, and Jupyter checkpoints."""
    
    supported_files = []
    # Directories to exclude
    exclude_dirs = {
        '.venv', '__pycache__', '.git', 'venv', 'env', 'node_modules', 'Results',
        '.ipynb_checkpoints'  # Jupyter notebook checkpoints
    }
    
    # Supported file extensions
    csv_extensions = {'.csv'}
    tree_extensions = {'.trees', '.ts', '.newick', '.nwk', '.tree', '.nexus'}
    all_extensions = csv_extensions | tree_extensions
    
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
            # Check if file has supported extension (case insensitive)
            file_lower = file.lower()
            if any(file_lower.endswith(ext) for ext in all_extensions):
                # Get relative path from current directory
                rel_path = os.path.join(root, file)
                
                # Determine file type
                if any(file_lower.endswith(ext) for ext in csv_extensions):
                    file_type = "CSV"
                else:
                    file_type = "TREE"
                
                supported_files.append((rel_path, file_type))
                
    return sorted(supported_files, key=lambda x: (x[1], x[0]))  # Sort by type then name

def get_tree_parameters(file_path):
    """
    Interactive prompts to get taxon names and outgroup for tree files.
    Only prompts if the file is not a TreeSequence (.trees/.ts) file.
    
    Returns:
        tuple: (taxon_names, outgroup) or (None, None) for TreeSequence files
    """
    file_path = Path(file_path)
    
    # TreeSequence files don't need taxon names and outgroup
    if file_path.suffix.lower() in ['.trees', '.ts']:
        print("✅ TreeSequence file detected - population information will be extracted automatically")
        return None, None
    
    print(f"\n📄 Newick/Nexus file detected: {file_path.name}")
    print("For Newick/Nexus files, you need to specify taxon names and outgroup.")
    
    # Get taxon names
    print("\nEnter taxon names (space-separated):")
    print("Example: O P1 P2 P3")
    while True:
        taxon_input = input("Taxon names: ").strip()
        if taxon_input:
            taxon_names = taxon_input.split()
            if len(taxon_names) >= 3:
                break
            else:
                print("Error: Need at least 3 taxon names for topology analysis")
        else:
            print("Error: Please enter taxon names")
    
    # Get outgroup
    print(f"\nAvailable taxa: {', '.join(taxon_names)}")
    while True:
        outgroup = input("Enter outgroup taxon name: ").strip()
        if outgroup in taxon_names:
            break
        else:
            print(f"Error: '{outgroup}' not found in taxon names. Please choose from: {', '.join(taxon_names)}")
    
    print(f"✅ Configuration: Taxa={taxon_names}, Outgroup={outgroup}")
    return taxon_names, outgroup

def main():
    # Print instructions
    print("\n=== TWISSTNTERN Analysis Tool ===")
    print("\nSupported file formats:")
    print("📊 CSV files: Pre-computed topology weights")
    print("🌳 Tree files: TreeSequence (.trees, .ts), Newick (.newick, .nwk, .tree), Nexus (.nexus)")
    print("\nAvailable files:")
    
    supported_files = list_supported_files()
    
    if not supported_files:
        print("No supported files found in the current directory or subdirectories.")
        return

    # Display numbered list of files with types
    for i, (file_path, file_type) in enumerate(supported_files, 1):
        print(f"{i}. [{file_type}] {file_path}")
    
    # Get file selection
    while True:
        try:
            choice = input("\nEnter the number of the file you want to analyze (or 'q' to quit): ")
            if choice.lower() == 'q':
                return
            choice = int(choice)
            if 1 <= choice <= len(supported_files):
                input_file, file_type = supported_files[choice - 1]
                break
            print(f"Please enter a number between 1 and {len(supported_files)}")
        except ValueError:
            print("Please enter a valid number")
    
    # Get tree parameters if needed
    taxon_names = None
    outgroup = None
    if file_type == "TREE":
        taxon_names, outgroup = get_tree_parameters(input_file)
    
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
        if file_type == "TREE":
            print(f"File type: {file_type}")
            if taxon_names:
                print(f"Taxa: {taxon_names}")
                print(f"Outgroup: {outgroup}")
        print("This may take a few moments...")
        
        # Run analysis - FIXED: now unpacking 3 values
        results, fundamental_results, csv_file_used = twisstntern.run_analysis(
            file=input_file, 
            granularity=granularity,
            taxon_names=taxon_names,
            outgroup=outgroup
        )
        
        # Print results
        print("\nAnalysis complete!")
        print(f"Data file used: {csv_file_used}")
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
        if file_type == "CSV":
            print("Please make sure the file is a valid CSV with three columns (T1, T2, T3)")
        else:
            print("Please make sure the tree file is valid and parameters are correct")
        # Print more detailed error information
        import traceback
        print("\nDetailed error information:")
        traceback.print_exc()

if __name__ == "__main__":
    main() 