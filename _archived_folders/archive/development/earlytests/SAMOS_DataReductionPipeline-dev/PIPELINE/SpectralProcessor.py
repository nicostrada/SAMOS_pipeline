from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .Spectroscopy.wavelength import WavelengthCalibration

from .SAMOS_mods import (classify_spectroscopic_data,
                    add_wcs_keys,
                    identify_targets,
                    trace_targets)
                    #record_trace_information,
                    #save_extracted) 
                    #extraction ---to be added to SAMOS_mods after I get it working

from .SAMOS_mods import (NoMatchFound,
                    ReferenceData)
                    #NoTargetException

import sys
import os
import argparse
import astropy.units as u
import logging
from ccdproc import CCDData

import matplotlib.pyplot as plt
import warnings

"""

 Extract spectra from slit cutouts (after all of the other initial calibration
steps), and perform the wavelength calibration.

"""

