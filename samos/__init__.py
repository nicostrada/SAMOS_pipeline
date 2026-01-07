"""
SAMOS Data Reduction Pipeline

A modular pipeline for reducing multi-object spectroscopy and imaging data
from the SAMOS instrument.

Subpackages
-----------
core
    Core data processing modules (mosaic, cosmic rays, calibration)
spectroscopy
    Spectroscopic reduction tools (traces, wavelength calibration, extraction)
imaging
    Imaging reduction tools (photometry, astrometry)
utils
    Utility functions (I/O, display, configuration)
pipeline
    Automated pipeline orchestration

Examples
--------
>>> import samos
>>> from samos.core import mosaic, cosmic_rays
>>> from samos.utils import display, io

>>> # Read and display a SAMOS frame
>>> hdu = mosaic.read_sami_mosaic('science_001.fits')
>>> display.display_image(hdu.data, zmin=500, zmax=3000)

>>> # Remove cosmic rays
>>> cleaned = cosmic_rays.remove_cosmic_rays(hdu.data)
"""

__version__ = "2.0.0"
__author__ = "SAMOS Team"

# Import main subpackages for easier access
from . import core
from . import spectroscopy
from . import imaging
from . import utils
from . import pipeline

__all__ = ['core', 'spectroscopy', 'imaging', 'utils', 'pipeline']
