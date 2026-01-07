"""
SAMOS Core Modules

Core image processing functionality for SAMOS data reduction.

Modules
-------
mosaic
    Multi-CCD mosaic assembly
cosmic_rays
    Cosmic ray detection and removal
calibration
    Bias, dark, and flat field calibration (future)
"""

from . import mosaic
from . import cosmic_rays

__all__ = ['mosaic', 'cosmic_rays']
