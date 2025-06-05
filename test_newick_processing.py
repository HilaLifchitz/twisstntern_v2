# Test script for newick processing with population labels
import sys

sys.path.append("twisstntern")
from twisstntern.tree_processing import (
    newick_to_twisst_weights,
    trees_to_twisst_weights_unified,
)

# Read a few trees
with open("tree files/ts_pop_plain.newick", "r") as f:
    trees = [line.strip() for line in f.readlines()[:5]]  # First 5 trees

print(f"Read {len(trees)} trees")
print(f"First tree sample: {trees[0][:100]}...")

# Test with population names (as they appear in the trees)
print("\n" + "=" * 50)
print("TESTING WITH POPULATION NAMES")
print("=" * 50)

taxon_names = ["O", "P1", "P2", "P3"]
outgroup = "O"

try:
    df = newick_to_twisst_weights(
        trees, taxon_names=taxon_names, outgroup=outgroup, verbose=True
    )
    print("\n✓ SUCCESS with population names!")
    print(f"DataFrame shape: {df.shape}")
    print("First few rows:")
    print(df.head())
except Exception as e:
    print(f"✗ Error with population names: {e}")
    import traceback

    traceback.print_exc()

# Test the unified function
print("\n" + "=" * 50)
print("TESTING UNIFIED FUNCTION")
print("=" * 50)

try:
    df2 = trees_to_twisst_weights_unified(
        "tree files/ts_pop_plain.newick",
        taxon_names=["O", "P1", "P2", "P3"],
        outgroup="O",
        output_file="test_newick_weights.csv",
        verbose=True,
    )
    print("\n✓ SUCCESS with unified function!")
    print(f"DataFrame shape: {df2.shape}")
    print("First few rows:")
    print(df2.head())
except Exception as e:
    print(f"✗ Error with unified function: {e}")
    import traceback

    traceback.print_exc()

print("\nTest completed!")
