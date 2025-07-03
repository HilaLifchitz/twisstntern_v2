"""
Configuration module for twisstntern.
Handles loading and validation of configuration parameters.

The configuration system uses YAML format for its readability and support for comments.
Configuration files should be structured as follows:

STEP 1: Setting up demography

1. Split Times:
   - List of population split events with:
     * "time": Time of split in generations before present
     * "derived_pop1": Name of the first derived population for the split
     * "derived_pop2": Name of the second derived population for the split
     * "ancestral_pop": Name of the ancestral population for the split

2. Population Definitions:
   - List of populations with their properties:
     * "name": Population identifier
     * "Ne": Effective population size (in individuals of the specified ploidy)
     * "sample_size": Number of samples to take (default: 10) - only for extant populations
     * OPTIONAL: "growth_rate": Population growth rate (default: 0.0)

3. Migration Rates:
   * Dictionary of migration rates between populations
   * Format: "source>destination": rate
   * Rates are proportions (e.g., 0.001 means 0.1% of individuals replaced by migrants per generation)
   * The migration matrix excludes the most ancestral population (ANC) and includes only extant populations (p1, p2, p3, O) and ancestral populations (p12, p123).

STEP 2: Coalescent simulation
1. Simulation Mode:
   - "simulation_mode": Choose between "locus" or "chromosome"
   - "locus" mode: Simulates independent non-recombining loci
   - "chromosome" mode: Simulates a recombining chromosome

2. OPTIONAL Parameters:
   - "mutation_rate": Mutation rate per base per generation 
   - "seed": Random seed for reproducibility 

3. Locus Mode Parameters (required if mode is "locus"):
   - "n_loci": Number of loci / windows to simulate
   - "locus_length": Length of each locus in base pairs

4. Chromosome Mode Parameters (required if mode is "chromosome"):
   - "chromosome_length": Total length of chromosome in base pairs
   - "rec_rate": Recombination rate per base per generation

Example configuration (YAML):
    populations:
      - name: "p1"
        Ne: 10000
        growth_rate: 0.0
        sample_size: 10
    splits:
      - time: 1000
        derived_pop1: "p1"
        derived_pop2: "p2"
        ancestral_pop: "p12"
    migration:
      p1>p2: 0.001
      p2>p1: 0.001
    simulation_mode: "locus"
    n_loci: 10000
    locus_length: 10000
"""

import yaml
from dataclasses import dataclass
from typing import List, Optional, Dict, Any
import re


def is_ancestral_population(pop_name):
    """
    Returns True if the population name is ancestral (e.g., 'p12', 'p123', 'p23', 'ANC'), False otherwise.
    """
    if pop_name == "ANC":
        return True
    match = re.fullmatch(r"p\d{2,}", pop_name)
    return match is not None


@dataclass
class Population:
    """
    Represents a population in the demographic model.

    Attributes:
        name: Unique identifier for the population
        Ne: Effective population size (in individuals of the specified ploidy)
        growth_rate: Population growth rate (optional, default: 0.0)
        sample_size: Number of samples to take from this population (only for extant populations)
        label: Human-readable label for the population (optional)
        is_ancestral: Whether this is an ancestral population (cannot be sampled)
    """

    name: str
    Ne: float
    growth_rate: float = 0.0
    sample_size: Optional[int] = None
    label: Optional[str] = None
    is_ancestral: bool = False


@dataclass
class Split:
    """
    Represents a population split event in the demographic model.

    Attributes:
        time: Time of split in generations before present
        derived_pop1: Name of the first derived population for the split
        derived_pop2: Name of the second derived population for the split
        ancestral_pop: Name of the ancestral population for the split
    """

    time: float
    derived_pop1: str
    derived_pop2: str
    ancestral_pop: str


########################################################################################
# this is the main class that loads the config file and validates the parameters
# it is used in the simulation.py file to run the simulation
########################################################################################


class Config:
    """
    Configuration handler for twisstntern simulations.

    This class loads and validates configuration parameters from a YAML file.
    It supports two simulation modes:
    1. Locus mode: Simulates independent non-recombining loci
    2. Chromosome mode: Simulates a recombining chromosome

    The configuration file should be in YAML format and include all necessary
    parameters for the chosen simulation mode.
    """

    def __init__(self, config_file: str):
        """
        Initialize configuration from a YAML file.

        Args:
            config_file: Path to the configuration file

        Raises:
            FileNotFoundError: If the configuration file doesn't exist
            yaml.YAMLError: If the file is not valid YAML
            ValueError: If required parameters are missing or invalid
        """
        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)

        # Load population labels
        self.population_labels = self.config.get("population_labels", {})

        # Load ploidy
        self.ploidy = self.config.get(
            "ploidy", 1
        )  # Default to haploid if not specified

        # Load populations
        self.populations = []
        for pop in self.config.get("populations", []):
            # Get population label if available
            label = self.population_labels.get(pop["name"])

            # Use the new function to determine if this is an ancestral population
            is_ancestral = is_ancestral_population(pop["name"])

            # Only set sample_size for extant populations
            sample_size = None if is_ancestral else pop.get("sample_size", 10)

            self.populations.append(
                Population(
                    name=pop["name"],
                    Ne=pop["Ne"],
                    growth_rate=pop.get(
                        "growth_rate", 0.0
                    ),  # Optional, defaults to 0.0
                    sample_size=sample_size,
                    label=label,
                    is_ancestral=is_ancestral,
                )
            )

        # Load splits
        self.splits = []
        for split in self.config.get("splits", []):
            self.splits.append(
                Split(
                    time=split["time"],
                    derived_pop1=split["derived_pop1"],
                    derived_pop2=split["derived_pop2"],
                    ancestral_pop=split["ancestral_pop"],
                )
            )

        # Load migration rates
        self.migration = self.config.get("migration", {})

        # Build migration matrix
        self.migration_matrix = self._build_migration_matrix()

        # Common parameters
        self.mutation_rate = self.config.get("mutation_rate", 1e-8)
        self.seed = self.config.get("seed")

        # Initialize mode-specific parameters first
        self.n_loci = None
        self.locus_length = None
        self.chromosome_length = None
        self.rec_rate = None

        # Simulation parameters - set the initial mode and load parameters
        self._simulation_mode = self.config.get("simulation_mode", "locus")
        self._reload_mode_parameters()

        # Validate configuration
        self._validate()

    def _build_migration_matrix(self):
        """
        Build a migration matrix from the migration rates.

        Returns:
            A 2D array where each row represents outgoing migration rates from a population to all others.
            The matrix excludes the most ancestral population (ANC) and includes only extant populations (p1, p2, p3, O) and ancestral populations (p12, p123).
        """
        # Define the order of populations for the matrix
        matrix_order = ["p1", "p2", "p3", "O", "p12", "p123"]
        n = len(matrix_order)
        matrix = [[0.0 for _ in range(n)] for _ in range(n)]

        # print("\nMigration dictionary from YAML:")
        # print("-------------------------------")
        # for source_dest, rate in self.migration.items():
        #     print(f"{source_dest}: {rate}")

        # print("\nBuilding migration matrix...")
        # print("----------------------------")
        for source_dest, rate in self.migration.items():
            try:
                source, dest = source_dest.split(">")
                if source in matrix_order and dest in matrix_order:
                    source_idx = matrix_order.index(source)
                    dest_idx = matrix_order.index(dest)
                    matrix[source_idx][dest_idx] = rate
                #  print(f"Added migration rate {rate} from {source} to {dest} at position [{source_idx}][{dest_idx}]")
            except ValueError:
                raise ValueError(f"Invalid migration rate key format: {source_dest}")
        return matrix

    def get_population_label(self, pop_name: str) -> str:
        """
        Get the human-readable label for a population.

        Args:
            pop_name: Name of the population

        Returns:
            The population label if available, otherwise the population name
        """
        return self.population_labels.get(pop_name, pop_name)

    def _reload_mode_parameters(self):
        """
        Reload mode-specific parameters after simulation_mode has been changed.
        This allows for proper mode overriding after initialization.
        """
        # Reset all mode-specific parameters
        self.n_loci = None
        self.locus_length = None
        self.chromosome_length = None
        self.rec_rate = None

        # Load mode-specific parameters based on current simulation_mode
        if self.simulation_mode == "locus":
            if "n_loci" not in self.config:
                raise ValueError("n_loci is required for locus mode")
            if "locus_length" not in self.config:
                raise ValueError("locus_length is required for locus mode")
            self.n_loci = self.config["n_loci"]
            self.locus_length = self.config["locus_length"]

        elif self.simulation_mode == "chromosome":
            if "chromosome_length" not in self.config:
                raise ValueError("chromosome_length is required for chromosome mode")
            if "rec_rate" not in self.config:
                raise ValueError("rec_rate is required for chromosome mode")
            self.chromosome_length = float(self.config["chromosome_length"])
            self.rec_rate = float(self.config["rec_rate"])

        # Re-validate with new parameters
        self._validate()

    @property
    def simulation_mode(self):
        """Get the simulation mode."""
        return self._simulation_mode

    @simulation_mode.setter
    def simulation_mode(self, value):
        """Set the simulation mode and reload mode-specific parameters."""
        self._simulation_mode = value
        self._reload_mode_parameters()

    def _validate(self):
        """
        Validate the configuration parameters.

        This method checks that:
        1. The simulation mode is valid
        2. Required parameters are present and valid for the chosen mode
        3. Population definitions are valid
        4. Split events are valid
        5. Migration rates are valid

        Raises:
            ValueError: If any validation check fails
        """
        # Validate simulation mode
        if self.simulation_mode not in ["locus", "chromosome"]:
            raise ValueError("simulation_mode must be 'locus' or 'chromosome'")

        # Validate common parameters
        if self.mutation_rate < 0:
            raise ValueError("mutation_rate must be non-negative")

        # Validate locus mode parameters
        if self.simulation_mode == "locus":
            if self.locus_length is None:
                raise ValueError("locus_length is required for locus mode")
            if self.locus_length <= 0:
                raise ValueError("locus_length must be positive")

        # Validate chromosome mode parameters
        elif self.simulation_mode == "chromosome":
            if self.chromosome_length is None:
                raise ValueError("chromosome_length is required for chromosome mode")
            if self.rec_rate is None:
                raise ValueError("rec_rate is required for chromosome mode")
            if self.chromosome_length <= 0:
                raise ValueError("chromosome_length must be positive")
            if self.rec_rate < 0:
                raise ValueError("rec_rate must be non-negative")

        # Validate populations
        if not self.populations:
            raise ValueError("At least one population must be defined")

        # Validate population sizes and sample sizes
        for pop in self.populations:
            self._validate_population(pop.__dict__)

        # Validate splits
        if not self.splits:
            raise ValueError("At least one population split must be defined")

        for split in self.splits:
            if split.time < 0:
                raise ValueError("Split time must be non-negative")
            if split.derived_pop1 not in [p.name for p in self.populations]:
                raise ValueError(
                    f"First derived population {split.derived_pop1} not found"
                )
            if split.derived_pop2 not in [p.name for p in self.populations]:
                raise ValueError(
                    f"Second derived population {split.derived_pop2} not found"
                )
            if split.ancestral_pop not in [p.name for p in self.populations]:
                raise ValueError(
                    f"Ancestral population {split.ancestral_pop} not found"
                )

        # Validate migration rates
        for source_dest, rate in self.migration.items():
            if rate < 0:
                raise ValueError(f"Migration rate {rate} must be non-negative")
            try:
                source, dest = source_dest.split(">")
                if source not in [p.name for p in self.populations]:
                    raise ValueError(f"Source population {source} not found")
                if dest not in [p.name for p in self.populations]:
                    raise ValueError(f"Destination population {dest} not found")
            except ValueError:
                raise ValueError(f"Invalid migration rate key format: {source_dest}")

    def _validate_population(self, pop: Dict[str, Any]) -> None:
        """Validate a single population configuration."""
        required_fields = ["name", "Ne"]
        for field in required_fields:
            if field not in pop:
                raise ValueError(f"Population missing required field: {field}")

        # Validate Ne
        try:
            ne_val = float(pop["Ne"])
        except (ValueError, TypeError):
            raise ValueError(f"Population {pop['name']} Ne must be a number")
        if ne_val <= 0:
            raise ValueError(f"Population {pop['name']} Ne must be a positive number")
        pop["Ne"] = ne_val

        # Validate growth_rate if present
        if "growth_rate" in pop:
            try:
                growth_rate = float(pop["growth_rate"])
                if growth_rate < 0:
                    raise ValueError(
                        f"Population {pop['name']} growth_rate must be non-negative"
                    )
                pop["growth_rate"] = growth_rate  # Convert to float
            except (ValueError, TypeError):
                raise ValueError(
                    f"Population {pop['name']} growth_rate must be a number"
                )

        # Validate sample_size if present
        if "sample_size" in pop and pop["sample_size"] is not None:
            if not isinstance(pop["sample_size"], int) or pop["sample_size"] <= 0:
                raise ValueError(
                    f"Population {pop['name']} sample_size must be a positive integer"
                )
            if is_ancestral_population(pop["name"]):
                raise ValueError(
                    f"Sample size should not be specified for ancestral population {pop['name']}"
                )
            else:
                pop["sample_size"] = pop["sample_size"]  # Keep the original sample_size

        # Validate label
        if "label" in pop:
            if not isinstance(pop["label"], str):
                raise ValueError(f"Population {pop['name']} label must be a string")
        else:
            pop["label"] = pop[
                "name"
            ]  # Default label to population name if not provided

        # Validate is_ancestral
        if "is_ancestral" in pop:
            if not isinstance(pop["is_ancestral"], bool):
                raise ValueError(
                    f"Population {pop['name']} is_ancestral must be a boolean"
                )
        else:
            pop["is_ancestral"] = False  # Default to False if not provided

        # Validate growth_rate only for extant populations
        if (
            pop["name"] not in ["p12", "p123", "ANC"]
            and pop.get("growth_rate") is not None
        ):
            if pop["growth_rate"] < 0:
                raise ValueError(f"Growth rate for {pop['name']} must be non-negative")

    def _validate_taxon_names(self):
        """
        Validate the taxon names for the simulation.

        This method checks that:
        1. All populations with sample size > 0 are included in the taxonNames list.

        Raises:
            ValueError: If any validation check fails
        """
        taxonNames = [
            pop.name
            for pop in self.populations
            if not pop.is_ancestral
            and (pop.sample_size is not None and pop.sample_size > 0)
        ]
        if not taxonNames:
            raise ValueError("No populations with sample size > 0 found")
