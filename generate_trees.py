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
# PARAMETERS + Setting and plotting the demographic model
####################################################################################


# KEEP THIS CONSTANT -- SHOULD BE KEPT SMALL FOR TESTING
nP1 = 10  # number of samples in population 1
nP2 = 10  # number of samples in population 2
nP3 = 10  # number of samples in population 3
n0 = 10  # number of samples in the outgroup

# Provide divergence times
# WE KEEP CONSTANT FOR THIS RUN
t1 = 100  # time of split between population 1 and population 2
t2 = 200  # time of split between population (1,2) and population 3
t3 = 300  # time of split between population (1,2,3) and the outgroup 0

# Provide migration rate
m = 0.00001  # migration rate between population 2 and population 3
mig_rate = m

Ne = 1000  # population size we keep constant for all populations for simplicity
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

# Plot the demographic model
ml.rcParams["figure.figsize"] = (8.0, 5.0)
graph = msprime.Demography.to_demes(demography)
fig, ax = plt.subplots()  # use plt.rcParams["figure.figsize"]
demesdraw.tubes(graph, ax=ax, seed=1)
plt.show()



######################################################################################
# GENERATING THE TREE SEQUENCES
######################################################################################

# For the LOCUS usage -- ts1 is a generator of tskit.TreeSequence objects

#######################################################
# Paramters
num_replicates = 1000
#######################################################

ts1 = msprime.sim_ancestry(
    samples={"O": n0, "P1": nP1, "P2": nP2, "P3": nP3},
    demography=demography,
    num_replicates=num_replicates,
    ploidy=1,  # Use haploid samples like DaSh-bash approach
)

# For the CHROMOSOME usage -- ts is a tskit.TreeSequence object
#######################################################
# Paramters
sequence_length = 100000000
recombination_rate = 0.00000001
#######################################################

ts = msprime.sim_ancestry(
    samples={"O": n0, "P1": nP1, "P2": nP2, "P3": nP3},
    demography=demography,
    sequence_length=sequence_length,
    recombination_rate=recombination_rate,
    ploidy=1,  # Use haploid samples like DaSh-bash approach
)



######################################################################################
# FUNCTIONS FOR POPULATION-LABELED NEWICK TREES
######################################################################################


def get_population_map(ts):
    """
    Create a mapping from sample ID to unique population sample name for a TreeSequence.

    Args:
        ts: msprime TreeSequence object

    Returns:
        dict: Mapping from sample_id (int) to population_sample_name (str) like "P1_1", "P1_2", etc.
    """
    pop_map = {}

    # Create mapping from population ID to population name
    pop_id_to_name = {}
    for pop_id in range(ts.num_populations):
        pop = ts.population(pop_id)
        if pop.metadata and "name" in pop.metadata:
            pop_name = pop.metadata["name"]
        else:
            pop_name = str(pop_id)  # fallback to ID if no name
        pop_id_to_name[pop_id] = pop_name

    # Count samples within each population to create unique identifiers
    pop_sample_counts = {}
    
    # Map each sample to its unique population sample name
    for sample_id in ts.samples():
        node = ts.node(sample_id)
        pop_name = pop_id_to_name[node.population]
        
        # Increment counter for this population
        if pop_name not in pop_sample_counts:
            pop_sample_counts[pop_name] = 0
        pop_sample_counts[pop_name] += 1
        
        # Create unique sample identifier like "P1_1", "P1_2", etc.
        unique_sample_name = f"{pop_name}_{pop_sample_counts[pop_name]}"
        pop_map[sample_id] = unique_sample_name

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
    if not newick.endswith(";"):
        newick += ";"

    return newick


######################################################################################
# DASHA'S APPROACH - EXACT MATCH TO WORKING VERSION
######################################################################################

def ts_newick_rename_dict_hun(): # replaces ,1:	--> ,P1_1:
    """
    Creates dictionary to replace ,1: patterns with ,P1_1: patterns in Newick strings.
    This matches the exact approach from DaSh-bash/LittorinaBrooding that works with twisst.
    """
    sample_names = dict()
    
    #Total number of empirical samples
    ntotal=nP1+nP2+nP3+n0

    #Creating list of newick names (patterns to find)
    ts_name=list()
    for i in range(1,ntotal+1):
        ts_name.append(","+str(i)+":")

    #Creating list of population samples (replacement patterns)
    n0_names=list()
    for i in range(1,n0+1):
        n0_names.append(",O_"+str(i)+":")
    
    
    nP1_names=list()
    for i in range(1,nP1+1):
        nP1_names.append(",P1_"+str(i)+":")

    nP2_names=list()
    for i in range(1,nP2+1):
        nP2_names.append(",P2_"+str(i)+":")

    nP3_names=list()
    for i in range(1,nP3+1):
        nP3_names.append(",P3_"+str(i)+":")



    pop_names=nP1_names+nP2_names+nP3_names+n0_names

    #Using dictionary comprehension to convert lists to dictionary
    sample_names = {ts_name[i]: pop_names[i] for i in range(len(ts_name))}
    
    return sample_names

def ts_newick_rename_dict_tens(): #  replaces (1: --> (P1_1:
    """
    Creates dictionary to replace (1: patterns with (P1_1: patterns in Newick strings.
    This matches the exact approach from DaSh-bash/LittorinaBrooding that works with twisst.
    """
    sample_names = dict()
    
    #Total number of empirical samples
    ntotal=nP1+nP2+nP3+n0

    #Creating list of newick names (patterns to find)
    ts_name=list()
    for i in range(1,ntotal+1):
        ts_name.append("("+str(i)+":")

    #Creating list of population samples (replacement patterns)   
    n0_names=list()
    for i in range(1,n0+1):
        n0_names.append("("+"O_"+str(i)+":")
    
    nP1_names=list()
    for i in range(1,nP1+1):
        nP1_names.append("("+"P1_"+str(i)+":")

    nP2_names=list()
    for i in range(1,nP2+1):
        nP2_names.append("("+"P2_"+str(i)+":")

    nP3_names=list()
    for i in range(1,nP3+1):
        nP3_names.append("("+"P3_"+str(i)+":")


    pop_names=nP1_names+nP2_names+nP3_names+n0_names

    #Using dictionary comprehension to convert lists to dictionary
    sample_names = {ts_name[i]: pop_names[i] for i in range(len(ts_name))}
    
    return sample_names


def save_ts_as_dasha_newick(genealogies, output_path):
    """
    Save TreeSequence genealogies as Newick format using DaSh-bash approach.
    This creates files that work properly with twisst.
    
    Args:
        genealogies: iterable of TreeSequence objects
        output_path: path to save the Newick file
    """
    with open(output_path, "w") as file:
        for replicate_index, ts in enumerate(genealogies):
            for t in ts.trees():
                newick = t.newick(precision=1)
                
                # Apply DaSh-bash renaming approach
                replace_strings_tens = ts_newick_rename_dict_tens()            
                replace_strings_hun = ts_newick_rename_dict_hun()
                
                # Apply replacements in correct order
                for word in replace_strings_tens.items():
                    newick = newick.replace(str(word[0]), str(word[1]))
                for word in replace_strings_hun.items():
                    newick = newick.replace(str(word[0]), str(word[1]))
                    
                file.write(newick + "\n")


######################################################################################
# SAVING THE TREE SEQUENCES
######################################################################################

# Create output directory if it doesn't exist
output_dir = "tree files"
os.makedirs(output_dir, exist_ok=True)

##################################################################################
# SAVING THE TREE SEQUENCES -- CHROMOSOME USAGE (ts)
##################################################################################


# 1. Save as tskit.TreeSequence object
ts.dump(os.path.join(output_dir, "CHROM.trees"))


# 2. Save Newick trees without metadata headers but with population labels

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


save_ts_as_plain_population_newick(ts, os.path.join(output_dir, "CHROM_pop_plain.newick"))

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


save_ts_as_population_nexus(ts, os.path.join(output_dir, "CHROM_pop.nexus"))


######################################################################################
# SAVING THE TREE SEQUENCES -- LOCUS USAGE (ts1)
######################################################################################

# Convert ts1 generator to a list so we can use it multiple times
ts1_list = list(ts1)

# 4. Save a simple version for twisst testing - just the trees, no metadata
simple_combined_path = os.path.join(output_dir, "simple_replicates_LOCUS.newick")
with open(simple_combined_path, "w") as f:
    for rep_idx, ts_rep in enumerate(ts1_list):
        pop_map = get_population_map(ts_rep)
        for tree in ts_rep.trees():
            labeled_newick = create_newick_with_sample_labels(tree, pop_map)
            f.write(labeled_newick + "\n")


# 5. Save using DaSh-bash approach for comparison and testing
dasha_combined_path = os.path.join(output_dir, "dasha_approach_LOCUS.newick")
save_ts_as_dasha_newick(ts1_list, dasha_combined_path)

######################################################################################
# PRINT INFORMATION FOR TWISST USAGE
######################################################################################

# Get the taxon names and outgroup for the first tree sequence
taxon_names=["O","P1","P2","P3"]
outgroup="O"

print("\n" + "=" * 60)
print("TWISST USAGE INFORMATION")
print("=" * 60)
print(f"Taxon names to use with twisst: {taxon_names}")
print(f"Suggested outgroup: '{outgroup}'")
print("\nGenerated files:")
print("  - ts_pop_plain.newick: Best for twisst (plain Newick with pop labels)")
print("  - simple_replicates_for_twisst.newick: Multiple replicates for twisst")
print("  - dasha_approach_for_twisst.newick: DaSh-bash approach (may work better)")
print(f"\nExample twisst usage:")
print(f"  taxon_names = {taxon_names}")
print(f"  outgroup = '{outgroup}'")
print("=" * 60)

# Also save this info to a file
info_path = os.path.join(output_dir, "twisst_info.txt")
with open(info_path, "w") as f:
    f.write("TWISST USAGE INFORMATION\n")
    f.write("========================\n\n")
    f.write(f"Taxon names: {taxon_names}\n")
    f.write(f"Outgroup: {outgroup}\n\n")
    f.write("Generated files for twisst:\n")
    f.write("  - ts_pop_plain.newick: Single TreeSequence with population labels\n")
    f.write("  - simple_replicates_for_twisst.newick: Multiple replicates\n")
    f.write("  - dasha_approach_for_twisst.newick: DaSh-bash approach (exact match)\n\n")
    f.write("Example usage in tree_processing.py:\n")
    f.write(
        f"  newick_to_twisst_weights(trees, taxon_names={taxon_names}, outgroup='{outgroup}')\n"
    )
