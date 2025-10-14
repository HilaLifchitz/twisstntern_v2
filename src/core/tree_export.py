"""Shared helpers for exporting tree data to Newick files."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable

import tskit


def get_population_map(ts: tskit.TreeSequence) -> Dict[int, str]:
    """Create a mapping from sample identifiers to population-aware labels."""

    pop_map: Dict[int, str] = {}

    pop_id_to_name: Dict[int, str] = {}
    for pop_id in range(ts.num_populations):
        pop = ts.population(pop_id)
        if pop.metadata and "name" in pop.metadata:
            pop_name = pop.metadata["name"]
        else:
            pop_name = str(pop_id)
        pop_id_to_name[pop_id] = pop_name

    pop_sample_counts: Dict[str, int] = {}

    for sample_id in ts.samples():
        node = ts.node(sample_id)
        pop_name = pop_id_to_name[node.population]
        count = pop_sample_counts.get(pop_name, 0) + 1
        pop_sample_counts[pop_name] = count
        pop_map[sample_id] = f"{pop_name}_{count}"

    return pop_map


def create_newick_with_sample_labels(tree: tskit.Tree, pop_map: Dict[int, str]) -> str:
    """Build a Newick string that uses population labels for samples."""

    def _get_newick_recursive(node: int) -> str:
        if tree.is_sample(node):
            return pop_map.get(node, f"Sample_{node}")

        children = list(tree.children(node))
        if not children:
            return ""

        child_strings = []
        for child in children:
            child_str = _get_newick_recursive(child)
            if not child_str:
                continue
            branch_length = tree.branch_length(child)
            if branch_length is not None and branch_length > 0:
                child_str += f":{branch_length}"
            child_strings.append(child_str)

        if len(child_strings) == 1:
            return child_strings[0]

        return f"({','.join(child_strings)})"

    root = tree.root
    newick = _get_newick_recursive(root)

    if not newick.endswith(";"):
        newick += ";"

    return newick


def save_ts_CHROM_as_newick(ts: tskit.TreeSequence, output_path: Path | str) -> str:
    """Write all marginal trees from a recombining chromosome to Newick."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    pop_map = get_population_map(ts)
    with output_path.open("w") as handle:
        for tree in ts.trees():
            handle.write(create_newick_with_sample_labels(tree, pop_map) + "\n")

    return str(output_path)


def save_ts_LOCUS_as_plain_newick(
    ts_list: Iterable[tskit.TreeSequence], output_path: Path | str
) -> str:
    """Write a collection of TreeSequences (locus mode) to a Newick file."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as handle:
        for ts in ts_list:
            pop_map = get_population_map(ts)
            for tree in ts.trees():
                handle.write(create_newick_with_sample_labels(tree, pop_map) + "\n")

    return str(output_path)


def save_newick_strings(newick_strings: Iterable[str], output_path: Path | str) -> str:
    """Persist an iterable of raw Newick strings to disk."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as handle:
        for tree_str in newick_strings:
            text = tree_str.strip()
            if not text:
                continue
            if not text.endswith(";"):
                text += ";"
            handle.write(text + "\n")

    return str(output_path)

