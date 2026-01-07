import json
from pathlib import Path

import sys

import os

import re 

import ccdproc as ccdp
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.colors import LogNorm
#from convenience_functions import show_image


import matplotlib.gridspec as gridspec

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import importlib

from itertools import groupby


#general os
import zipfile
import urllib.request

#general plotting
from matplotlib.patches import Rectangle

params={'legend.fontsize':'18','axes.labelsize':'18',
		'axes.titlesize':'18','xtick.labelsize':'18',
		'ytick.labelsize':'18','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
plt.rcParams.update(params)
plt.rcParams.update({'figure.max_open_warning': 0})

#table/math handling
import pandas as pd
import numpy as np
np.seterr(all='ignore')  # hides irrelevant warnings about divide-by-zero, etc

#astropy
import astropy
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table, QTable
from astropy.io import fits,ascii
from astropy.nddata import StdDevUncertainty
from astropy.modeling import models, fitting
from astropy.visualization import quantity_support,astropy_mpl_style, simple_norm
from astropy import constants as const
from astropy.stats import mad_std
from astropy.nddata import CCDData
from astropy.visualization import hist

#specutils
import specutils
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import gaussian_smooth
from specutils.fitting import fit_generic_continuum
from specutils.fitting import find_lines_derivative
from specutils.fitting import find_lines_threshold
from specutils.fitting import fit_lines
from specutils.manipulation import noise_region_uncertainty
from specutils.analysis import centroid
from specutils.analysis import line_flux
from specutils.analysis import equivalent_width
from specutils.analysis import template_comparison


#######################

from ginga import GingaPlugin
#from ginga.misc import Widgets

from ginga.qtw.QtHelp import QtGui, QtCore
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
from ginga.mplw.ImageViewCanvasMpl import ImageViewCanvas
from ginga.mplw.FigureCanvasQt import FigureCanvas
from ginga.misc import log
from ginga.util.loader import load_data
from matplotlib.figure import Figure
	
from ginga import toolkit, colors
from ginga.gw import Viewers
from ginga.misc.log import get_logger 
from ginga.canvas.CanvasObject import get_canvas_types
from ginga import AstroImage, colors

########################






class SamosPlugin(GingaPlugin.LocalPlugin):

	def __init__(self, fv, fitsimage):


		super(SamosPlugin, self).__init__(fv, fitsimage)


		self.header_keys = ['naxis',
						 'date',
						 'slit',
						 'date-obs',
						 'obstype',
						 'object',
						 'exptime',
						 'obsra',
						 'obsdec',
						 'grating',
						 'cam_targ',
						 'grt_targ',
						 'filter',
						 'filter2',
						 'gain',
						 'rdnoise',
						 'ccdsum',
						 'wavmode']



		raw_data_path = "/Users/danakoeppe/PipelineProjects/Make_FITS_for_SAMOS/SAMI_reduction/GOODMAN_raw_data/"

		if not os.path.exists("plugin_test_products"):
			os.mkdir("plugin_test_products")

		self.raw_data_dir = raw_data_path
		self.procdir = "plugin_test_products"

		

	def build_gui(self, container):

		top = Widgets.VBox()
		top.set_border_width(4)

		# this is a little trick for making plugins that work either in
		# a vertical or horizontal orientation.  It returns a box container,
		# a scroll widget and an orientation ('vertical', 'horizontal')
		vbox, sw, orientation = Widgets.get_oriented_box(container)
		vbox.set_border_width(4)
		vbox.set_spacing(2)

		# Take a text widget to show some instructions
		self.msgFont = self.fv.getFont("sansFont", 12)
		tw = Widgets.TextArea(wrap=True, editable=False)
		tw.set_font(self.msgFont)
		self.tw = tw

		# Frame for instructions and add the text widget with another
		# blank widget to stretch as needed to fill emp
		fr = Widgets.Frame("Instructions")
		vbox2 = Widgets.VBox()
		vbox2.add_widget(tw)
		vbox2.add_widget(Widgets.Label(''), stretch=1)
		fr.set_widget(vbox2)
		vbox.add_widget(fr, stretch=0)

		# Add a spacer to stretch the rest of the way to the end of the
		# plugin space
		spacer = Widgets.Label('')
		vbox.add_widget(spacer, stretch=1)

		# scroll bars will allow lots of content to be accessed
		top.add_widget(sw, stretch=1)

		# A button box that is always visible at the bottom
		btns = Widgets.HBox()
		btns.set_spacing(3)

		# Add a close button for the convenience of the user
		btn = Widgets.Button("Close")
		btn.add_callback('activated', lambda w: self.close())
		btns.add_widget(btn, stretch=0)
		btns.add_widget(Widgets.Label(''), stretch=1)
		top.add_widget(btns, stretch=0)

		# Add our GUI to the container
		container.add_widget(top, stretch=1)
		# NOTE: if you are building a GUI using a specific widget toolkit
		# (e.g. Qt) GUI calls, you need to extract the widget or layout
		# from the non-toolkit specific container wrapper and call on that
		# to pack your widget, e.g.:
		#cw = container.get_widget()
		#cw.addWidget(widget, stretch=1)


		"""
		# Add dropdowns for choosing groups of organized raw data.
		wfilepick = QtGui.QComboBox()

		for name in self.drawcolors:
			wdrawcolor.addItem(name)
		index = self.drawcolors.index('lightblue')
		wdrawcolor.setCurrentIndex(index)
		wdrawcolor.activated.connect(self.set_drawparams)
		self.wdrawcolor = wdrawcolor
		"""

		hbox.addStretch(1)
		org_data_button = Widgets.Button("Organize Raw Data")
		wquit.add_callback('organize_data', self.organize_raw_data)
		hbox.add_widget(org_data_button, stretch=0)


	def close(self):
		"""
		Example close method.  You can use this method and attach it as a
		callback to a button that you place in your GUI to close the plugin
		as a convenience to the user.
		"""
		chname = self.fv.get_channel_name(self.fitsimage)
		self.fv.stop_local_plugin(chname, str(self))
		return True

	def start(self):
		"""
		This method is called just after ``build_gui()`` when the plugin
		is invoked.  This method may be called many times as the plugin is
		opened and closed for modal operations.  This method may be omitted
		in many cases.
		"""
		self.tw.set_text("""This plugin doesn't do anything interesting.""")
		self.resume()

	def pause(self):
		"""
		This method is called when the plugin loses focus.
		It should take any actions necessary to stop handling user
		interaction events that were initiated in ``start()`` or
		``resume()``.
		This method may be called many times as the plugin is focused
		or defocused.  It may be omitted if there is no user event handling
		to disable.
		"""
		pass

	def resume(self):
		"""
		This method is called when the plugin gets focus.
		It should take any actions necessary to start handling user
		interaction events for the operations that it does.
		This method may be called many times as the plugin is focused or
		defocused.  The method may be omitted if there is no user event
		handling to enable.
		"""
		pass

	def stop(self):
		"""
		This method is called when the plugin is stopped.
		It should perform any special clean up necessary to terminate
		the operation.  The GUI will be destroyed by the plugin manager
		so there is no need for the stop method to do that.
		This method may be called many  times as the plugin is opened and
		closed for modal operations, and may be omitted if there is no
		special cleanup required when stopping.
		"""
		pass

	def redo(self):
		"""
		This method is called when the plugin is active and a new
		image is loaded into the associated channel.  It can optionally
		redo the current operation on the new image.  This method may be
		called many times as new images are loaded while the plugin is
		active.  This method may be omitted.
		"""
		pass

		pass

	def organize_raw_data(self):
		

		SNight = SN.SAMOSNight(raw_data_dir=raw_data_path,
					   obsid="fake_samos_test",
					   proc_dir=self.procdir,
					   LOG_FILENAME=LOG_FILENAME,
					   ignore_bias=False,
					   ignore_flats=False)

		SNight()

		self.data_buckets = SNight.data_buckets



	def __str__(self):
		"""
		This method should be provided and should return the lower case
		name of the plugin.
		"""
		return 'SamosPlugin'
