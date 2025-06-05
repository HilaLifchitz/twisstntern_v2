# This script is used to generate different type of tree sequences in different formats
# using msprime- so that we can test where "tree_processing.py" works

import msprime
import tskit
import matplotlib.pyplot as plt
import matplotlib as ml
import demesdraw
import os
import ete3

# we first generate the demographic model for msprime
# we can change the parameters here to generate different tree sequences

####################################################################################
# PARAMETERS
####################################################################################


# KEEP THIS CONSTANT -- SHOULD BE KEPT SMALL FOR TESTING
nP1 = 10 # number of samples in population 1
nP2 = 10 # number of samples in population 2
nP3 = 10 # number of samples in population 3
n0 = 10 # number of samples in the outgroup

# Provide divergence times
# WE KEEP CONSTANT FOR THIS RUN
t1 = 100 # time of split between population 1 and population 2
t2 = 200 # time of split between population (1,2) and population 3
t3 = 300 # time of split between population (1,2,3) and the outgroup 0

# Provide migration rate
m = 0.0001 # migration rate between population 2 and population 3
mig_rate = m

Ne = 1000 # population size we keep constant for all populations for simplicity
# Provide population sizes
NeP1 = Ne
NeP2 = Ne
NeP3 = Ne
NeO = Ne
NeP12 = Ne
NeP123 = Ne
NeANC = Ne


# A demography where (1,2) first coalesce, with our given Ne

demography = msprime.Demography()

# initializing populations
demography.add_population(name="O", initial_size=NeO)
demography.add_population(name="P1", initial_size=NeP1)
demography.add_population(name="P2", initial_size=NeP2)
demography.add_population(name="P3", initial_size=NeP3)
demography.add_population(name="P13", initial_size=NeP12)
demography.add_population(name="P123", initial_size=NeP123)
demography.add_population(name="ANC", initial_size=NeANC)

# adding split times
demography.add_population_split(time=t1, derived=["P1", "P2"], ancestral="P13")
demography.add_population_split(time=t2, derived=["P13", "P3"], ancestral="P123")
demography.add_population_split(time=t3, derived=["P123", "O"], ancestral="ANC")

# setting up gene flow
demography.set_migration_rate("P2", "P3", mig_rate)

#Plot the demographic model
ml.rcParams['figure.figsize'] = (8.0, 5.0)
graph = msprime.Demography.to_demes(demography)
fig, ax = plt.subplots()  # use plt.rcParams["figure.figsize"]
demesdraw.tubes(graph, ax=ax, seed=1)
plt.show()


######################################################################################
# FUNCTIONS FOR POPULATION-LABELED NEWICK TREES
######################################################################################

def get_population_map(ts):
    """
    Create a mapping from sample ID to population name for a TreeSequence.
    
    Args:
        ts: msprime TreeSequence object
        
    Returns:
        dict: Mapping from sample_id (int) to population_name (str)
    """
    pop_map = {}
    
    # Create mapping from population ID to population name
    pop_id_to_name = {}
    for pop_id in range(ts.num_populations):
        pop = ts.population(pop_id)
        if pop.metadata and 'name' in pop.metadata:
            pop_name = pop.metadata['name']
        else:
            pop_name = str(pop_id)  # fallback to ID if no name
        pop_id_to_name[pop_id] = pop_name
    
    # Map each sample to its population name
    for sample_id in ts.samples():
        node = ts.node(sample_id)
        pop_name = pop_id_to_name[node.population]
        pop_map[sample_id] = pop_name
    
    return pop_map

def create_newick_with_sample_labels(tree, pop_map):
    """
    Create a Newick string from a tskit tree with proper sample labels.
    
    Args:
        tree: tskit Tree object
        pop_map: dict mapping sample_id to population_name
        
    Returns:
        str: Newick string with population labels for samples
    """
    def _get_newick_recursive(node):
        """Recursively build Newick string"""
        if tree.is_sample(node):
            # This is a sample - use population label
            pop_label = pop_map.get(node, f"Sample_{node}")
            return pop_label
        else:
            # This is an internal node - recurse on children
            children = tree.children(node)
            if len(children) == 0:
                return ""
            
            child_strings = []
            for child in children:
                child_str = _get_newick_recursive(child)
                if child_str:  # Only add if not empty
                    # Add branch length if available
                    branch_length = tree.branch_length(child)
                    if branch_length > 0:
                        child_str += f":{branch_length}"
                    child_strings.append(child_str)
            
            if len(child_strings) == 1:
                return child_strings[0]
            else:
                return f"({','.join(child_strings)})"
    
    # Start from root
    root = tree.root
    newick = _get_newick_recursive(root)
    
    # Add semicolon at the end
    if not newick.endswith(';'):
        newick += ';'
    
    return newick

def relabel_newick_with_populations(newick_string, pop_map):
    """
    Relabel a Newick tree string to use population names instead of sample IDs.
    This version handles both sample IDs and node names properly.
    
    Args:
        newick_string (str): Original Newick tree with sample IDs
        pop_map (dict): Mapping from sample_id to population_name
        
    Returns:
        str: Newick tree with population labels
    """
    try:
        # Parse the tree
        tree = ete3.Tree(newick_string)
        
        # Relabel only leaf nodes
        for leaf in tree.get_leaves():
            leaf_name = leaf.name
            
            # Try to extract node ID from 'nXX' format  
            if leaf_name.startswith('n') and leaf_name[1:].isdigit():
                node_id = int(leaf_name[1:])
                if node_id in pop_map:
                    leaf.name = pop_map[node_id]
                    continue
            
            # Try direct sample ID mapping
            try:
                sample_id = int(leaf_name)
                if sample_id in pop_map:
                    leaf.name = pop_map[sample_id]
                    continue
            except ValueError:
                pass
                
            # If we can't map it, leave it as is
        
        # Return the relabeled tree in Newick format
        return tree.write(format=1)  # format=1 includes branch lengths
        
    except Exception as e:
        print(f"Error relabeling tree: {e}")
        return newick_string  # Return original if there's an error

def get_taxon_names_and_outgroup(ts):
    """
    Extract taxon names and determine outgroup for twisst analysis.
    
    Args:
        ts: msprime TreeSequence object
        
    Returns:
        tuple: (taxon_names_list, outgroup_name)
    """
    pop_map = get_population_map(ts)
    
    # Get unique population names
    taxon_names = sorted(list(set(pop_map.values())))
    
    # Use the first population as outgroup (which should be "0" for population O)
    outgroup = "0"  # This corresponds to population "O"
    
    return taxon_names, outgroup

######################################################################################
# GENERATING THE TREE SEQUENCES
######################################################################################

# For the Chromosome usage -- ts is a tskit.TreeSequence object
ts = msprime.sim_ancestry(
    samples={"O": n0, "P1": nP1, "P2": nP2, "P3": nP3},
    demography=demography,
    sequence_length=100000000,
    recombination_rate=0.00000001,
    ploidy=1,
)
 
num_replicates=20
# For the locus usage its impossible to save the trees as a tskit.TreeSequence object!
# but we can save this as a Newick/nexus files
ts1 = msprime.sim_ancestry(
    samples={"O": n0, "P1": nP1, "P2": nP2, "P3": nP3},   
    demography=demography,
    num_replicates=num_replicates
)


######################################################################################
# SAVING THE TREE SEQUENCES
######################################################################################

# Create output directory if it doesn't exist
output_dir = "tree files"
os.makedirs(output_dir, exist_ok=True)


# 1. Save as tskit.TreeSequence object
ts.dump(os.path.join(output_dir, "ts.trees"))

# 2a. Saving marginal trees as Newick with population labels (one tree per line)
def save_ts_as_population_labeled_newick(ts, output_path):
    """
    Saves all marginal trees from a TreeSequence as Newick format with population labels.

    Parameters:
    - ts: msprime.TreeSequence object
    - output_path: full path to the output .newick file
    """
    pop_map = get_population_map(ts)
    
    with open(output_path, "w") as f:
        for i, tree in enumerate(ts.trees()):
            # Create Newick with proper sample labels
            labeled_newick = create_newick_with_sample_labels(tree, pop_map)
            f.write(f"[Tree {i} @ interval {tree.interval}]\n")
            f.write(labeled_newick + "\n\n")

save_ts_as_population_labeled_newick(ts, os.path.join(output_dir, "ts_pop_labeled.newick"))

# 2b. Save Newick trees without metadata headers but with population labels
def save_ts_as_plain_population_newick(ts, output_path):
    """
    Saves marginal trees as Newick format with population labels, no interval annotations.
    This format is ideal for twisst.
    """
    pop_map = get_population_map(ts)
    
    with open(output_path, "w") as f:
        for tree in ts.trees():
            labeled_newick = create_newick_with_sample_labels(tree, pop_map)
            f.write(labeled_newick + "\n")

save_ts_as_plain_population_newick(ts, os.path.join(output_dir, "ts_pop_plain.newick"))

# 2c. Save the OLD versions for comparison
def save_ts_as_newick(ts, output_path):
    """
    Saves all marginal trees from a TreeSequence as Newick format.

    Parameters:
    - ts: msprime.TreeSequence object
    - output_path: full path to the output .newick file
    """
    with open(output_path, "w") as f:
        for i, tree in enumerate(ts.trees()):
            f.write(f"[Tree {i} @ interval {tree.interval}]\n")
            f.write(tree.as_newick() + "\n\n")

save_ts_as_newick(ts, os.path.join(output_dir, "ts_old.newick"))

def save_ts_as_plain_newick(ts, output_path):
    """
    Saves marginal trees as Newick format without interval annotations.
    """
    with open(output_path, "w") as f:
        for tree in ts.trees():
            f.write(tree.as_newick() + "\n")

save_ts_as_plain_newick(ts, os.path.join(output_dir, "ts_old_plain.newick"))

# 3. Save as Nexus format with population labels
def save_ts_as_population_nexus(ts, output_path):
    """
    Save all trees from a TreeSequence in Nexus format with population labels.
    Each tree includes interval metadata in its label.
    """
    pop_map = get_population_map(ts)
    
    with open(output_path, "w") as f:
        f.write("#NEXUS\n\n")
        f.write("Begin trees;\n")
        for i, tree in enumerate(ts.trees()):
            interval = f"[{tree.interval[0]:.1f},{tree.interval[1]:.1f}]"
            labeled_newick = create_newick_with_sample_labels(tree, pop_map)
            f.write(f"  Tree TREE{i+1} {interval} = {labeled_newick};\n")
        f.write("End;\n")

save_ts_as_population_nexus(ts, os.path.join(output_dir, "ts_pop.nexus"))

# Old nexus for comparison
def save_ts_as_nexus(ts, output_path):
    """
    Save all trees from a TreeSequence in Nexus format.
    Each tree includes interval metadata in its label.
    """
    with open(output_path, "w") as f:
        f.write("#NEXUS\n\n")
        f.write("Begin trees;\n")
        for i, tree in enumerate(ts.trees()):
            interval = f"[{tree.interval[0]:.1f},{tree.interval[1]:.1f}]"
            f.write(f"  Tree TREE{i+1} {interval} = {tree.as_newick()};\n")
        f.write("End;\n")

save_ts_as_nexus(ts, os.path.join(output_dir, "ts_old.nexus"))


# 4. ts1 is a generator object-- we save each replicate separately

# Subdirectory to store each replicate's files
rep_dir = os.path.join(output_dir, "replicates")
os.makedirs(rep_dir, exist_ok=True)
ts1=list(ts1)
for i in range(len(ts1)):
    base_name = f"rep{i}"
    
    #4a. Save in .trees format
    ts_path = os.path.join(rep_dir, f"{base_name}.trees")
    ts_rep = ts1[i]
    ts_rep.dump(ts_path)

    #4b. Save in .newick format with population labels
    newick_path = os.path.join(rep_dir, f"{base_name}_pop.newick")
    save_ts_as_population_labeled_newick(ts_rep, newick_path)

    #4c. Save in .nexus format with population labels
    nexus_path = os.path.join(rep_dir, f"{base_name}_pop.nexus")
    save_ts_as_population_nexus(ts_rep, nexus_path)

# 5. Save all replicates in a single Newick file with population labels
combined_newick_path = os.path.join(output_dir, "all_replicates_ts1_pop.newick")

with open(combined_newick_path, "w") as f:
    for rep_idx, ts_rep in enumerate(ts1):
        pop_map = get_population_map(ts_rep)
        for tree_idx, tree in enumerate(ts_rep.trees()):
            labeled_newick = create_newick_with_sample_labels(tree, pop_map)
            f.write(f"[Rep {rep_idx} | Tree {tree_idx} @ interval {tree.interval}]\n")
            f.write(labeled_newick + "\n\n")

# 6. Save a simple version for twisst testing - just the trees, no metadata
simple_combined_path = os.path.join(output_dir, "simple_replicates_for_twisst.newick")
with open(simple_combined_path, "w") as f:
    for rep_idx, ts_rep in enumerate(ts1):
        pop_map = get_population_map(ts_rep)
        for tree in ts_rep.trees():
            labeled_newick = create_newick_with_sample_labels(tree, pop_map)
            f.write(labeled_newick + "\n")

######################################################################################
# PRINT INFORMATION FOR TWISST USAGE
######################################################################################

# Get the taxon names and outgroup for the first tree sequence
taxon_names, outgroup = get_taxon_names_and_outgroup(ts)

print("\n" + "="*60)
print("TWISST USAGE INFORMATION")
print("="*60)
print(f"Taxon names to use with twisst: {taxon_names}")
print(f"Suggested outgroup: '{outgroup}'")
print("\nGenerated files:")
print("  - ts_pop_plain.newick: Best for twisst (plain Newick with pop labels)")
print("  - simple_replicates_for_twisst.newick: Multiple replicates for twisst")
print(f"\nExample twisst usage:")
print(f"  taxon_names = {taxon_names}")
print(f"  outgroup = '{outgroup}'")
print("="*60)

# Also save this info to a file
info_path = os.path.join(output_dir, "twisst_info.txt")
with open(info_path, "w") as f:
    f.write("TWISST USAGE INFORMATION\n")
    f.write("========================\n\n")
    f.write(f"Taxon names: {taxon_names}\n")
    f.write(f"Outgroup: {outgroup}\n\n")
    f.write("Generated files for twisst:\n")
    f.write("  - ts_pop_plain.newick: Single TreeSequence with population labels\n")
    f.write("  - simple_replicates_for_twisst.newick: Multiple replicates\n\n")
    f.write("Example usage in tree_processing.py:\n")
    f.write(f"  newick_to_twisst_weights(trees, taxon_names={taxon_names}, outgroup='{outgroup}')\n")
