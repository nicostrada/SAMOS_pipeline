"""
SAMOS Pipeline Setup Script

Installation:
    pip install -e .                    # Development mode
    pip install -e .[dev]               # With development dependencies
    pip install -e .[pypeit]            # With PypeIt integration
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_file = Path(__file__).parent / "README.md"
if readme_file.exists():
    long_description = readme_file.read_text()
else:
    long_description = "SAMOS Data Reduction Pipeline"

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
if requirements_file.exists():
    with open(requirements_file) as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]
else:
    requirements = []

setup(
    name="samos-pipeline",
    version="2.0.0",
    author="SAMOS Team",
    author_email="",
    description="Data reduction pipeline for SAMOS multi-object spectrograph and imager",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/samos/pipeline",  # Update with actual URL
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    python_requires=">=3.10",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "sphinx>=5.0",
            "sphinx-rtd-theme>=1.0",
            "black>=23.0",
            "flake8>=6.0",
        ],
        "pypeit": [
            "pypeit>=1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "samos-reduce-spectroscopy=scripts.run_spectroscopy_pipeline:main",
        ],
    },
    include_package_data=True,
    package_data={
        "samos": [
            # Note: calibration_data archived - SAMOS line lists now in pypeit_wrapper.py
            # Add future SAMOS-specific data files here if needed
        ],
    },
)
