from setuptools import setup, find_packages

setup(
    name="twisstntern",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "scipy>=1.7.0",
        "sympy>=1.7.0",
    ],
    extras_require={
        "dev": [
            "black>=25.1.0",
            "pytest>=7.0.0",
            "flake8>=6.0.0",
            "mypy>=1.0.0",
        ],
    },
    python_requires=">=3.10",
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for ternary plot analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/twisstntern",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
    ],
) 