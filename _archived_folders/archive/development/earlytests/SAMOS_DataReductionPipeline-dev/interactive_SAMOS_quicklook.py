import json
from pathlib import Path

import sys
sys.path.insert(0,'Users/danakoeppe/PipelineProjects/SAMOS_DRP/SAMOS_Draft/PIPELINE')
#print(sys.path)

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

quantity_support();  # auto-recognizes units on matplotlib plots

################
###ginga####
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


#########DATA REDUC MODS##############

print(os.getcwd())
from PIPELINE import SAMOS_NIGHT as SN



######################################


class FitsViewer(QtGui.QMainWindow):

	def __init__(self, logger):
		super(FitsViewer, self).__init__()


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



		self.logger = get_logger("my viewer", log_stderr=False, log_file='ginga.log', level="DEBUG")
		self.drawcolors = colors.get_colors()
		self.dc = get_canvas_types()
		

		fig = Figure()
		w = FigureCanvas(fig)

		self.fig = fig
		self.canvas = FigureCanvas(self.fig)

		vbox2 = QtGui.QWidget()
		layout = QtGui.QVBoxLayout()
		layout.addWidget(self.canvas, stretch=1)

		# create the ginga viewer and configure it
		fi = ImageViewCanvas(self.logger)
		fi.enable_autocuts('on')
		fi.set_autocut_params('zscale')
		fi.enable_autozoom('on')
		fi.enable_draw(False)
		fi.set_callback('drag-drop', self.drop_file_cb)
		fi.set_callback('cursor-changed', self.cursor_cb)
		fi.set_bg(0.2, 0.2, 0.2)
		fi.ui_set_active(True)
		self.fitsimage = fi
		fi.set_figure(fig)


		# enable some user interaction
		bd = fi.get_bindings()
		bd.enable_all(True)

		# canvas that we will draw on
		canvas = self.dc.DrawingCanvas()
		canvas.enable_draw(True)
		canvas.enable_edit(True)
		canvas.set_drawtype('rectangle', color='lightblue')
		canvas.set_surface(fi)
		self.canvas = canvas
		# add canvas to view
		fi.get_canvas().add(canvas)
		canvas.ui_set_active(True)
		canvas.register_for_cursor_drawing(fi)
		canvas.add_callback('draw-event', self.draw_cb)

		
		w.resize(512, 512)

		# add scrollbar interface around this viewer
		#si = ScrolledView(fi)

		vbox = QtGui.QVBoxLayout()
		vbox.setContentsMargins(QtCore.QMargins(2, 2, 2, 2))
		vbox.setSpacing(1)
		vbox.addWidget(w, stretch=0)

		self.readout = QtGui.QLabel("")
		vbox.addWidget(self.readout, stretch=0,
					   alignment=QtCore.Qt.AlignCenter)

		hbox = QtGui.QHBoxLayout()
		hbox.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))

		wdrawtype = QtGui.QComboBox()
		self.drawtypes = ['rectangle','square']#fi.get_drawtypes()
		for name in self.drawtypes:
			#if np.logical_or(name=='rectangle',name=='square'):
				
			wdrawtype.addItem(name)
		index = self.drawtypes.index('rectangle')
		wdrawtype.setCurrentIndex(index)
		wdrawtype.activated.connect(self.set_drawparams)
		self.wdrawtype = wdrawtype

		wdrawcolor = QtGui.QComboBox()
		for name in self.drawcolors:
			wdrawcolor.addItem(name)
		index = self.drawcolors.index('lightblue')
		wdrawcolor.setCurrentIndex(index)
		wdrawcolor.activated.connect(self.set_drawparams)
		self.wdrawcolor = wdrawcolor

		wfill = QtGui.QCheckBox("Fill")
		wfill.stateChanged.connect(self.set_drawparams)
		self.wfill = wfill

		walpha = QtGui.QDoubleSpinBox()
		walpha.setRange(0.0, 1.0)
		walpha.setSingleStep(0.1)
		walpha.setValue(1.0)
		walpha.valueChanged.connect(self.set_drawparams)
		self.walpha = walpha


		
		org_data_button = QtGui.QPushButton("Organize Raw Data")
		#org_data_button.add_callback('organize_data', self.organize_raw_data)
		org_data_button.clicked.connect(self.organize_raw_data)
		#hbox.add_widget(org_data_button, stretch=0)
		

		wclear = QtGui.QPushButton("Clear Canvas")
		wclear.clicked.connect(self.clear_canvas)


		wopen = QtGui.QPushButton("Open File")
		wopen.clicked.connect(self.open_file)


		wquit = QtGui.QPushButton("Quit")
		wquit.clicked.connect(self.quit)

		wgetimg = QtGui.QPushButton("Crop Data")
		wgetimg.clicked.connect(self.crop_image)

		woscan = QtGui.QPushButton("Subtract Overscan")
		woscan.clicked.connect(self.subtract_overscan)


		hbox.addStretch(1)
		for w in (wopen, wdrawtype, wdrawcolor, wfill,wgetimg,
				  QtGui.QLabel('Alpha:'), walpha, wclear, wquit,
				  org_data_button, woscan):
			hbox.addWidget(w, stretch=0)

		hw = QtGui.QWidget()
		hw.setLayout(hbox)
		vbox.addWidget(hw, stretch=0)

		
		mode = self.canvas.get_draw_mode()
		hbox = QtGui.QHBoxLayout()
		hbox.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))

		btn1 = QtGui.QRadioButton("Draw")
		btn1.setChecked(mode == 'draw')
		btn1.toggled.connect(lambda val: self.set_mode_cb('draw', val))
		btn1.setToolTip("Choose this to draw on the canvas")
		hbox.addWidget(btn1)

		btn2 = QtGui.QRadioButton("Edit")
		btn2.setChecked(mode == 'edit')
		btn2.toggled.connect(lambda val: self.set_mode_cb('edit', val))
		btn2.setToolTip("Choose this to edit things on the canvas")
		hbox.addWidget(btn2)


		btn3 = QtGui.QRadioButton("Pick")
		btn3.setChecked(mode == 'pick')
		btn3.toggled.connect(lambda val: self.set_mode_cb('pick', val))
		btn3.setToolTip("Choose this to pick things on the canvas")
		hbox.addWidget(btn3)


		hbox.addWidget(QtGui.QLabel(''), stretch=1)
		hw = QtGui.QWidget()
		hw.setLayout(hbox)
		vbox.addWidget(hw, stretch=0)


		vw = QtGui.QWidget()
		self.setCentralWidget(vw)
		vw.setLayout(vbox)

		



	def set_drawparams(self, kind):
		index = self.wdrawtype.currentIndex()
		kind = self.drawtypes[index]
		index = self.wdrawcolor.currentIndex()
		fill = (self.wfill.checkState() != 0)
		alpha = self.walpha.value()

		params = {'color': self.drawcolors[index],
				  'alpha': alpha,
				  }
				  #'x1'	 : self.canvas.D
		#if kind in ('circle', 'rectangle', 'polygon', 'triangle',
		#			'righttriangle', 'ellipse', 'square', 'box'):
		if kind in ('rectangle', 'square'):
			params['fill'] = fill
			params['fillalpha'] = alpha

		self.canvas.set_drawtype(kind, **params)

	def get_drawparams(self):
		
		print(self.doodle_dict)
		


	def clear_canvas(self):
		self.canvas.delete_all_objects()

	def load_file(self, filepath):
		image = load_data(filepath, logger=self.logger)
		self.fitsimage.set_image(image)
		self.setWindowTitle(filepath)
		rows = image.shape[0]
		cols = image.shape[1]
		imgtxt = "%s,%s"%(rows,cols)
		self.readout.setText(imgtxt)
		self.ccd = CCDData.read(filepath)
		print(rows,cols)

	def open_file(self):
		res = QtGui.QFileDialog.getOpenFileName(self, "Open FITS file",
												".", "FITS files (*.fits)")
		if isinstance(res, tuple):
			fileName = res[0]
		else:
			fileName = str(res)
		if len(fileName) != 0:
			self.load_file(fileName)

	def drop_file_cb(self, viewer, paths):
		fileName = paths[0]
		self.load_file(fileName)

	def cursor_cb(self, viewer, button, data_x, data_y):
		"""This gets called when the data position relative to the cursor
		changes.
		"""
		# Get the value under the data coordinates
		try:
			# We report the value across the pixel, even though the coords
			# change halfway across the pixel
			value = viewer.get_data(int(data_x + viewer.data_off),
									int(data_y + viewer.data_off))

		except Exception:
			value = None

		fits_x, fits_y = data_x + 1, data_y + 1

		# Calculate WCS RA
		try:
			# NOTE: image function operates on DATA space coords
			image = viewer.get_image()
			if image is None:
				# No image loaded
				return
			ra_txt, dec_txt = image.pixtoradec(fits_x, fits_y,
											   format='str', coords='fits')
		except Exception as e:
			#self.logger.warning("Bad coordinate conversion: %s" % (
			#	str(e)))
			ra_txt = 'BAD WCS'
			dec_txt = 'BAD WCS'

		text = "RA: %s  DEC: %s  X: %.2f  Y: %.2f  Value: %s" % (
			ra_txt, dec_txt, fits_x, fits_y, value)
		self.readout.setText(text)

	def quit(self, *args):
		self.logger.info("Attempting to shut down the application...")
		self.deleteLater()


	def set_mode_cb(self, mode, tf):
		self.logger.info("canvas mode changed (%s) %s" % (mode, tf))
		if not (tf is False):
			self.canvas.set_draw_mode(mode)
		return True


	def draw_cb(self, canvas, tag):
		obj = canvas.get_object_by_tag(tag)
		index = self.wdrawtype.currentIndex()
		kind = self.drawtypes[index]

		obj.add_callback('pick-down', self.pick_cb, 'down')
		obj.add_callback('pick-up', self.pick_cb, 'up')
		obj.add_callback('pick-move', self.pick_cb, 'move')
		obj.add_callback('pick-hover', self.pick_cb, 'hover')
		obj.add_callback('pick-enter', self.pick_cb, 'enter')
		obj.add_callback('pick-leave', self.pick_cb, 'leave')
		obj.add_callback('pick-key', self.pick_cb, 'key')
		obj.pickable = True
		obj.add_callback('edited', self.edit_cb)
		
		doodle_points = np.asarray(obj.get_points())
		doodle_txt = "X{:n},Y{:n} : ({:.3f}, {:.3f})"
		#print(doodle_points.shape)
		

		for i in range(doodle_points.shape[0]):
			print(doodle_txt.format(i,i,
									doodle_points[i][0],
									doodle_points[i][1]))

		self.drawing = obj




	def pick_cb(self, obj, canvas, event, pt, ptype):
		self.logger.info("pick event '%s' with obj %s at (%.2f, %.2f)" % (
			ptype, obj.kind, pt[0], pt[1]))
		return True

	def edit_cb(self, obj):
		self.logger.info("object %s has been edited" % (obj.kind))
		return True



	def crop_image(self):

		image = self.fitsimage.get_image()
		#print(type(image))
		image_data = image.get_data()

		draw_pts = self.drawing.get_points()
		print(self.drawing.get_points())

		xs = sorted(list(set([int(i[0]) for i in draw_pts])))
		ys = sorted(list(set([int(i[1]) for i in draw_pts])))

		print(xs,ys)

		cropped_image = image_data[ys[0]:ys[1],xs[0]:xs[1]]
		self.fitsimage.set_data(cropped_image)



	def subtract_overscan(self):

		image = self.fitsimage.get_image()

		raw_ccd = self.ccd 

		#print(raw_ccd.header['trimsec'])

		trim_section = raw_ccd.header['trimsec']
		x,y = trim_section.strip("[").strip("]").split(",")
		x1,x2 = x.split(":")
		os_x_start = x2
		os_x_end = raw_ccd.data.shape[1]
		osec = "[%s:%s,%s]"%(os_x_start,os_x_end,y)
		os_ccd = ccdp.subtract_overscan(raw_ccd,fits_section=osec,overscan_axis=1)
		tos_ccd = ccdp.trim_image(os_ccd,fits_section=trim_section)

		self.fitsimage.set_data(tos_ccd.data)
		self.ccd = tos_ccd

	def organize_raw_data(self):
		

		SNight = SN.SAMOSNight(raw_data_dir=self.raw_data_dir,
					   obsid="fake_samos_test",
					   proc_dir=self.procdir,
					   LOG_FILENAME='ginga_quicklook.log',
					   ignore_bias=False,
					   ignore_flats=False)

		SNight()

		self.data_buckets = SNight.data_buckets




def main(options, args):

	app = QtGui.QApplication(sys.argv)

	# ginga needs a logger.
	# If you don't want to log anything you can create a null logger by
	# using null=True in this call instead of log_stderr=True
	logger = log.get_logger("test", log_stderr=True, level=40)

	w = FitsViewer(logger)
	w.resize(524, 524)
	w.show()
	app.setActiveWindow(w)
	w.raise_()
	w.activateWindow()

	if len(args) > 0:
		w.load_file(args[0])

	app.exec_()


if __name__ == '__main__':
	print(sys.argv[1:])
	main(None, sys.argv[1:])

