"""
Logging Configuration Module for TWISSTNTERN

This module provides centralized logging configuration for the TWISSTNTERN package.
It sets up both console and file logging with appropriate formatting and levels.

Features:
- Console logging with colored output
- File logging to output directory
- Proper log level handling
- Thread-safe logging
- Automatic log file creation in output directory
"""

import logging
import sys
from pathlib import Path
from typing import Optional
import time


class ColoredFormatter(logging.Formatter):
    """Colored formatter for console output"""
    
    # Color codes for different log levels
    COLORS = {
        'DEBUG': '\033[36m',      # Cyan
        'INFO': '\033[32m',       # Green
        'WARNING': '\033[33m',    # Yellow
        'ERROR': '\033[31m',      # Red
        'CRITICAL': '\033[35m',   # Magenta
        'RESET': '\033[0m'        # Reset
    }
    
    def format(self, record):
        # Add color to the log level
        levelname = record.levelname
        if levelname in self.COLORS:
            record.levelname = f"{self.COLORS[levelname]}{levelname}{self.COLORS['RESET']}"
        
        return super().format(record)


def setup_logging(
    output_dir: Optional[str] = None,
    log_level: str = "INFO",
    console_output: bool = True,
    verbose: bool = False
) -> str:
    """
    Set up logging configuration for TWISSTNTERN.
    
    Args:
        output_dir: Directory where log file will be saved. If None, only console logging.
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        console_output: Whether to output logs to console
        verbose: If True, sets DEBUG level and detailed formatting
        
    Returns:
        Path to log file if created, None if only console logging
    """
    
    # Clear any existing handlers
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Set log level
    if verbose:
        log_level = "DEBUG"
    
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)
    root_logger.setLevel(numeric_level)
    
    # Create formatters
    if verbose:
        console_format = "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s"
        file_format = "%(asctime)s - %(name)s - %(levelname)s - %(module)s.%(funcName)s:%(lineno)d - %(message)s"
    else:
        console_format = "%(levelname)s - %(message)s"
        file_format = "%(asctime)s - %(levelname)s - %(message)s"
    
    # Console handler
    if console_output:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(numeric_level)
        
        if sys.stdout.isatty():  # Only use colors if output is a terminal
            console_formatter = ColoredFormatter(console_format)
        else:
            console_formatter = logging.Formatter(console_format)
            
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)
    
    # File handler (if output directory specified)
    log_file_path = None
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Create log file with timestamp
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        log_file_path = output_path / f"twisstntern_{timestamp}.log"
        
        file_handler = logging.FileHandler(log_file_path, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)  # Always save debug info to file
        
        file_formatter = logging.Formatter(file_format)
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)
    
    # Set specific loggers to appropriate levels
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)
    
    # Suppress verbose output from external libraries unless in debug mode
    if not verbose:
        logging.getLogger("msprime").setLevel(logging.WARNING)
        logging.getLogger("tskit").setLevel(logging.WARNING)
        logging.getLogger("ete3").setLevel(logging.WARNING)
    
    return str(log_file_path) if log_file_path else None


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger with the specified name.
    
    Args:
        name: Logger name (usually __name__)
        
    Returns:
        Configured logger instance
    """
    return logging.getLogger(name)


def log_system_info():
    """Log system and environment information"""
    logger = get_logger(__name__)
    
    import platform
    import sys
    
    logger.info("="*60)
    logger.info("TWISSTNTERN Analysis Session")
    logger.info("="*60)
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Platform: {platform.platform()}")
    logger.info(f"Working directory: {Path.cwd()}")
    
    # Log package versions if available
    try:
        import numpy as np
        logger.debug(f"NumPy version: {np.__version__}")
    except ImportError:
        pass
        
    try:
        import pandas as pd
        logger.debug(f"Pandas version: {pd.__version__}")
    except ImportError:
        pass
        
    try:
        import matplotlib
        logger.debug(f"Matplotlib version: {matplotlib.__version__}")
    except ImportError:
        pass


def log_analysis_start(input_file: str, output_dir: str, **kwargs):
    """Log the start of an analysis with parameters"""
    logger = get_logger(__name__)
    
    logger.info("Starting TWISSTNTERN analysis")
    logger.info(f"Input file: {input_file}")
    logger.info(f"Output directory: {output_dir}")
    
    # Log other parameters
    for key, value in kwargs.items():
        if value is not None:
            logger.info(f"{key}: {value}")


def log_analysis_complete(duration: float, output_files: list):
    """Log analysis completion with summary"""
    logger = get_logger(__name__)
    
    logger.info("="*60)
    logger.info("Analysis completed successfully!")
    logger.info(f"Total duration: {duration:.2f} seconds")
    logger.info("Generated files:")
    
    for file_path in output_files:
        if Path(file_path).exists():
            size = Path(file_path).stat().st_size
            logger.info(f"  ✓ {file_path} ({size:,} bytes)")
        else:
            logger.warning(f"  ⚠ {file_path} (file not found)")
    
    logger.info("="*60)


def log_simulation_config(config, overrides=None):
    """Log simulation configuration parameters"""
    logger = get_logger(__name__)
    
    logger.info("="*60)
    logger.info("SIMULATION CONFIGURATION")
    logger.info("="*60)
    
    # Basic simulation parameters
    logger.info(f"Simulation mode: {config.simulation_mode}")
    logger.info(f"Random seed: {config.seed}")
    
    # Population parameters
    logger.info(f"Number of populations: {len(config.populations)}")
    for pop_config in config.populations:
        # Use population label if available, otherwise just the name
        pop_display = pop_config.name
        if hasattr(pop_config, 'label') and pop_config.label and pop_config.label != pop_config.name:
            pop_display = f"{pop_config.name} ({pop_config.label})"
        
        if hasattr(pop_config, 'sample_size') and pop_config.sample_size is not None:
            logger.info(f"  {pop_display}: Ne={pop_config.Ne}, samples={pop_config.sample_size}")
        else:
            logger.info(f"  {pop_display}: Ne={pop_config.Ne} (ancestral)")
    
    # Migration parameters
    if hasattr(config, 'migration') and config.migration:
        logger.info("Migration rates:")
        for route, rate in config.migration.items():
            # Enhance migration logging with population labels if available
            if hasattr(config, 'get_population_label'):
                try:
                    source, dest = route.split('>')
                    source_label = config.get_population_label(source)
                    dest_label = config.get_population_label(dest)
                    if source_label != source or dest_label != dest:
                        # Show both short names and labels for clarity
                        route_display = f"{route} ({source_label} → {dest_label})"
                    else:
                        route_display = route
                except:
                    route_display = route
            else:
                route_display = route
            logger.info(f"  {route_display}: {rate}")
    
    # Simulation-specific parameters
    if config.simulation_mode == "chromosome":
        logger.info(f"Chromosome length: {config.chromosome_length}")
        logger.info(f"Recombination rate: {config.rec_rate}")
        if hasattr(config, 'mutation_rate'):
            logger.info(f"Mutation rate: {config.mutation_rate}")
    elif config.simulation_mode == "locus":
        logger.info(f"Number of loci: {config.n_loci}")
        logger.info(f"Locus length: {config.locus_length}")
    
    # Ploidy
    logger.info(f"Ploidy: {config.ploidy}")
    
    # Demographic events
    if hasattr(config, 'demographic_events') and config.demographic_events:
        logger.info("Demographic events:")
        for event in config.demographic_events:
            logger.info(f"  {event}")
    
    # Command-line overrides
    if overrides:
        logger.info("Command-line overrides applied:")
        for key, value in overrides.items():
            if value is not None:
                logger.info(f"  {key}: {value}")
    
    logger.info("="*60)


def log_error(error: Exception, context: str = ""):
    """Log an error with context"""
    logger = get_logger(__name__)
    
    if context:
        logger.error(f"Error in {context}: {str(error)}")
    else:
        logger.error(f"Error: {str(error)}")
    
    # Log traceback in debug mode
    logger.debug("Exception details:", exc_info=True)


def log_topologies(topos, simplified_topos, columns, logger, message_prefix="", population_labels=None):
    """
    Log topology information both as strings and ASCII trees to the log file ONLY.
    This avoids duplicating the console output that's already handled by print() statements.
    
    Args:
        topos: List of topology objects with get_ascii() method
        simplified_topos: List of simplified topology strings  
        columns: List of column names (T1, T2, T3, etc.)
        logger: Logger instance
        message_prefix: Optional prefix for log messages
        population_labels: Optional dict mapping population IDs to descriptive labels
    """
    # Get a file-only logger to avoid duplicating console output
    file_logger = logging.getLogger(f"{logger.name}.topologies")
    
    # Remove any console handlers from this specific logger
    for handler in file_logger.handlers[:]:
        if isinstance(handler, logging.StreamHandler) and handler.stream.name == '<stdout>':
            file_logger.removeHandler(handler)
    
    # Ensure it only logs to file by setting propagate=False and using parent's file handler
    file_logger.propagate = False
    
    # Add only the file handler from the parent logger
    root_logger = logging.getLogger()
    for handler in root_logger.handlers:
        if isinstance(handler, logging.FileHandler):
            file_logger.addHandler(handler)
            break
    
    file_logger.setLevel(logging.INFO)
    
    if message_prefix:
        file_logger.info(f"{message_prefix}")
    else:
        file_logger.info("Discovered topologies:")
    
    file_logger.info("="*50)
    
    # Always add population legend - extract all leaf names from topologies
    leaf_names = set()
    for topo in topos:
        leaf_names.update([leaf.name for leaf in topo.get_leaves()])
    
    # Sort leaf names for consistent output
    sorted_leaf_names = sorted(leaf_names)
    
    file_logger.info("Population Legend:")
    file_logger.info("-" * 20)
    
    if population_labels:
        # Use custom population labels if provided
        for leaf_name in sorted_leaf_names:
            if leaf_name in population_labels:
                file_logger.info(f"  {leaf_name} = {population_labels[leaf_name]}")
            else:
                file_logger.info(f"  {leaf_name} = Population {leaf_name}")
    else:
        # Create basic legend with numeric IDs
        for leaf_name in sorted_leaf_names:
            file_logger.info(f"  {leaf_name} = Population {leaf_name}")
    
    file_logger.info("")
    
    for i, (topo, simplified, col_name) in enumerate(zip(topos, simplified_topos, columns)):
        file_logger.info(f"{col_name}:")
        file_logger.info(f"  String: {simplified}")
        
        # Get ASCII representation
        try:
            ascii_tree = topo.get_ascii()
            # Log each line of the ASCII tree
            file_logger.info("  ASCII Tree:")
            for line in ascii_tree.split('\n'):
                if line.strip():  # Only log non-empty lines
                    file_logger.info(f"    {line}")
        except Exception as e:
            file_logger.warning(f"  Could not generate ASCII tree: {e}")
        
        file_logger.info("")  # Empty line for readability
    
    file_logger.info("="*50) 