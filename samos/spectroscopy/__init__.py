"""
SAMOS Spectroscopy Module

Provides functions for spectroscopic data reduction including trace detection,
extraction, and wavelength calibration.
"""

from . import trace_detection
from . import extraction
from . import pypeit_wrapper

__all__ = ['trace_detection', 'extraction', 'pypeit_wrapper']
