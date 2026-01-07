"""
SAMOS Utility Modules

General utility functions for I/O, display, and configuration.

Modules
-------
io
    FITS I/O and file management
display
    Image display and visualization
config
    Configuration file management (future)
"""

from . import io
from . import display

__all__ = ['io', 'display']
