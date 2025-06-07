#!/usr/bin/env python
# coding: utf-8

from setuptools import setup, find_packages
from pathlib import Path

# Read the contents of README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text() if (this_directory / "README.md").exists() else ""

setup(
    name="twisstntern",
    version="0.1.0",
    author="Hila Lifchitz",
    author_email="hila.lifchitz@ist.ac.at",
    description="A package for analyzing ternary data from topology weights",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/HilaLifchitz/twisstntern_v2",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "matplotlib>=3.4.0",
        "tskit>=0.4.0",
        "msprime>=1.0.0",
        "ete3>=3.1.0",
        "requests>=2.25.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "black",
            "flake8",
        ],
    },
    entry_points={
        "console_scripts": [
            "twisstntern=twisstntern.__main__:main",
        ],
    },
    include_package_data=True,
    package_data={
        "twisstntern": ["external/*"],
    },
)
