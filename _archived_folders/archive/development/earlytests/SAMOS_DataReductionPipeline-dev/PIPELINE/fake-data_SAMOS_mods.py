import os
import sys

from pathlib import Path


import json
import glob
import argparse

from astropy.nddata import CCDData
from astropy.visualization import hist
import ccdproc as ccdp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

from PIL import Image as P


from matplotlib.colors import LogNorm
#from convenience_functions import show_image

from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
import matplotlib.gridspec as gridspec

import warnings

from itertools import groupby

warnings.filterwarnings("ignore")

import logging

from astropy.stats import mad_std



import pandas as pd



def stitch_SAMI_amp_output(multi_ext_SAMI_fits):

	

