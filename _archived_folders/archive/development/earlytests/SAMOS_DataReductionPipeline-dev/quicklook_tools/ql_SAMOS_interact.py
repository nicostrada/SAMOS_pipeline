#! /usr/bin/env python
#
# example3_mpl.py -- Copy attributes from a Ginga Qt widget into a Matplotlib
#                       figure.
#
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
#
"""
   $ ./example3_mpl.py [fits file]
example3 displays a native ginga widget beside a matplotlib figure as two
panes.  A fits file can be dropped into the left pane and manipulated using
the standard Ginga interactive controls
see (http://ginga.readthedocs.io/en/latest/quickref.html).
Drop down boxes allow the color map to be changed.
The right pane has two buttons under it: pressing each button sets up a
different kind of plot in the mpl pane based on the current state of the
ginga pane.
You need Qt4 with python bindings (or pyside) installed to run this example.
"""
import sys
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

from ginga.qtw.ImageViewCanvasQt import ImageViewCanvas
from ginga.qtw.QtHelp import QtGui, QtCore
from ginga import cmap, imap
from ginga.misc import log
from ginga.util.loader import load_data
from ginga import AstroImage, colors
from ginga.canvas.CanvasObject import get_canvas_types


from matplotlib.patches import Rectangle
from matplotlib.widgets import Slider

from astropy.visualization import quantity_support,astropy_mpl_style, simple_norm
import ccdproc as ccdp
from astropy import constants as const
from astropy.stats import mad_std
from astropy.nddata import CCDData
from astropy.visualization import hist
from astropy.modeling import models, fitting

import numpy as np

STD_FORMAT = '%(asctime)s | %(levelname)1.1s | %(filename)s:%(lineno)d (%(funcName)s) | %(message)s'


class FitsViewer(QtGui.QMainWindow):

	def __init__(self, logger):
		super(FitsViewer, self).__init__()
		self.logger = logger

		

		# Add matplotlib color maps to our built in ones
		cmap.add_matplotlib_cmaps()
		self.cmaps = cmap.get_names()
		self.imaps = imap.get_names()

		self.drawcolors = colors.get_colors()
		self.dc = get_canvas_types()

		wd, ht = 500, 500

		###########begin vbox1###############


		# Create a Ginga widget
		fi = ImageViewCanvas(logger, render='widget')

		fi.enable_autocuts('on')
		fi.set_autocut_params('zscale')
		fi.enable_autozoom('on')
		fi.enable_draw(False)
		fi.set_callback('drag-drop', self.drop_file_cb)
		fi.set_bg(0.2, 0.2, 0.2)
		fi.ui_set_active(True)
		self.fitsimage = fi

		fi.show_color_bar(True)

		# enable various key and mouse controlled actions
		bd = fi.get_bindings()
		bd.enable_all(True)

		# canvas that we will draw on
		canvas1 = self.dc.DrawingCanvas()
		canvas1.enable_draw(True)
		canvas1.enable_edit(True)
		canvas1.set_drawtype('rectangle', color='lightblue')
		canvas1.set_surface(fi)
		self.canvas1 = canvas1
		# add canvas to view
		fi.get_canvas().add(canvas1)
		canvas1.ui_set_active(True)
		canvas1.register_for_cursor_drawing(fi)
		canvas1.add_callback('draw-event', self.draw_cb)

		self.cp_tag = 'compass'

		# pack widget into layout
		gingaw = fi.get_widget()
		gingaw.resize(wd, ht)



		vbox1 = QtGui.QWidget()
		layout = QtGui.QVBoxLayout()
		layout.addWidget(gingaw, stretch=1)


		#wopen = QtGui.QPushButton("Open File")
		#wopen.clicked.connect(self.open_file)

		# add buttons to layout
		hbox = QtGui.QHBoxLayout()
		hbox.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))
		hbox.addStretch(1)
		
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

		wopen = QtGui.QPushButton("Open File")
		wopen.clicked.connect(self.open_file)


		wquit = QtGui.QPushButton("Quit")
		wquit.clicked.connect(self.close)

		wsliceplot = QtGui.QPushButton("Plot Coadd")
		wsliceplot.setVisible(False)
		wsliceplot.clicked.connect(self.plot_slice)
		self.wsliceplot = wsliceplot

		wcropimg = QtGui.QPushButton("Crop Data")
		wcropimg.clicked.connect(self.crop_image)


		woscan = QtGui.QPushButton("Subtract Overscan")
		woscan.clicked.connect(self.subtract_overscan)
		self.woscan = woscan

		undo_woscan = QtGui.QPushButton("Undo Subtract Overscan")
		undo_woscan.setVisible(False)
		undo_woscan.clicked.connect(self.undo_subtract_overscan)
		self.undo_woscan = undo_woscan

		wclear = QtGui.QPushButton("Clear Canvas")
		wclear.clicked.connect(self.clear_canvas)


		hbox.addStretch(1)
		for w in (wopen, wdrawtype, wdrawcolor, wfill, wclear,
					walpha, QtGui.QLabel('Alpha:'),
					wquit):
		
			hbox.addWidget(w, stretch=0)

		

		hw = QtGui.QWidget()
		hw.setLayout(hbox)
		layout.addWidget(hw, stretch=0)

		hbox = QtGui.QHBoxLayout()
		hbox.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))
		hbox.addStretch(1)
		for w in (undo_woscan, woscan, wsliceplot, wcropimg):
		
			hbox.addWidget(w, stretch=0)

		hw = QtGui.QWidget()
		hw.setLayout(hbox)
		layout.addWidget(hw, stretch=0)

		vbox1.setLayout(layout)
		vbox1.setFixedWidth(700)

		#####################################


		###########begin vbox2###############

		# Create a matplotlib Figure
		#self.fig = matplotlib.figure.Figure(figsize=(wd, ht))
		self.fig1 = plt.figure()
		canvas2 = FigureCanvas(self.fig1)
		self.canvas2 = canvas2		

		vbox2 = QtGui.QWidget()
		layout2 = QtGui.QVBoxLayout()
		layout2.addWidget(self.canvas2, stretch=2)


		self.fig2 = plt.figure() #
		self.canvas3 = FigureCanvas(self.fig2) #




		vbox3 = QtGui.QWidget()
		#layout3 = QtGui.QVBoxLayout()
		layout2.addWidget(self.canvas3,stretch=2)
	
		hbox = QtGui.QHBoxLayout()
		hbox.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))


		wgetimg = QtGui.QPushButton("Fit Cross Disp.")
		wgetimg.setVisible(False)
		wgetimg.clicked.connect(self.fit_cross_disp_slice)
		self.wgetimg = wgetimg
		#wquit = QtGui.QPushButton("Quit")
		#wquit.clicked.connect(self.close)

		hbox.addStretch(1)
		#for w in (wgetimg):
		hbox.addWidget(wgetimg, stretch=0)

		hw = QtGui.QWidget()
		hw.setLayout(hbox)
		layout2.addWidget(hw, stretch=0)

		
		
		vbox2.setLayout(layout2)

		vbox = QtGui.QVBoxLayout()
		vbox.setContentsMargins(QtCore.QMargins(2, 2, 2, 2))
		vbox.setSpacing(1)

		#####################################



	
  

	
		#################end vbox2,3###################


		w = QtGui.QWidget()
		layout = QtGui.QHBoxLayout()
		layout.addWidget(vbox1, stretch=1.0)
		layout.addWidget(vbox2, stretch=1.0)
		#layout.addWidget(vbox3, stretch=1.0)
		w.setLayout(layout)

		vbox.addWidget(w, stretch=1)

		mode = self.canvas1.get_draw_mode()
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
		vw.setLayout(vbox)
		self.setCentralWidget(vw)


	def clear_canvas(self):
		self.fitsimage.delete_all_objects()

	def load_file(self, filepath):
		image = load_data(filepath, logger=self.logger)
		self.fitsimage.set_image(image)
		self.setWindowTitle(filepath)
		self.ccd = CCDData.read(filepath)

		"""
		# create compass
		try:
			try:
				self.fitsimage.delete_object_by_tag(self.cp_tag)
			except KeyError:
				pass

			width, height = image.get_size()
			x, y = width / 2.0, height / 2.0
			# radius we want the arms to be (approx 1/4 the largest dimension)
			radius = float(max(width, height)) / 4.0

			Compass = self.fitsimage.get_draw_class('compass')
			self.fitsimage.add(Compass(x, y, radius, color='skyblue',
									   fontsize=14), tag=self.cp_tag)
		except Exception as e:
			self.logger.warning("Can't calculate compass: %s" % (
				str(e)))
		"""


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
		filename = paths[0]
		self.load_file(filename)

	def closeEvent(self, ce):
		self.close()



	def set_drawparams(self, kind):
		index = self.wdrawtype.currentIndex()
		kind = self.drawtypes[index]
		index = self.wdrawcolor.currentIndex()
		fill = (self.wfill.checkState() != 0)
		alpha = self.walpha.value()

		params = {'color': self.drawcolors[index],
				  'alpha': alpha,
				  'fill' : fill,
				  }
				  #'x1'  : self.canvas.D
		#if kind in ('circle', 'rectangle', 'polygon', 'triangle',
		#           'righttriangle', 'ellipse', 'square', 'box'):
		if kind in ('rectangle', 'square'):
			params['fill'] = fill
			params['fillalpha'] = alpha

		self.canvas1.set_drawtype(kind, **params)


	def get_image(self):
		fi = self.fitsimage
		# clear previous image
		self.fig1.clf()

		ax = self.fig1.add_subplot(111)
		ax.autoscale(True, tight=True)

		x0, y0, x1, y1 = tuple(map(int, fi.get_datarect()))
		#extent = (x0, x1, y0, y1)

		image = fi.get_image()
		arr = image.cutout_data(x0, y0, x1, y1)

		extent = self.get_wcs_extent(image, x0, y0, x1, y1)

		# get cut levels
		loval, hival = fi.get_cut_levels()

		# make the equivalent color map for matplotlib
		cm = self.make_mpl_colormap(fi)

		# add the image to the figure
		interp = 'nearest'
		img = ax.imshow(arr, interpolation=interp, origin="lower",
						vmin=loval, vmax=hival, cmap=cm,
						aspect="equal", extent=extent)

		# add a colorbar
		self.fig1.colorbar(img, orientation='vertical')

		# force an update of the figure
		self.fig1.canvas.draw()

	def set_mode_cb(self, mode, tf):
		self.logger.info("canvas mode changed (%s) %s" % (mode, tf))
		if not (tf is False):
			self.canvas1.set_draw_mode(mode)
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

	def clear_canvas(self):
		self.canvas1.delete_all_objects()
		self.fig1.clf()
		self.fig2.clf()
		self.fig1.canvas.draw()
		self.fig2.canvas.draw()

	
	def updateLabel(self, value):

		self.label.setText(str(value))

	def popup_button(self):

		self.wpb = QtGui.QWidget()
		self.wpb.setGeometry(QtCore.QRect(100, 100, 400, 200))

		psf_gtmp = "     Model: {:s}\n amplitude: {:5.3f} \n      mean: {:5.3f} \n    stddev: {:5.3f}"
		

		psf_mtmp = "     Model: {:s}\n amplitude: {:5.3f} \n       x_0: {:5.3f} \n     gamma: {:5.3f} \n     alpha: {:5.3f}"

		gauss_txt = psf_gtmp.format("Gauss",self.psf_gauss_template.amplitude.value, 
								   self.psf_gauss_template.mean.value, 
								   self.psf_gauss_template.stddev.value)

		moff_txt = psf_mtmp.format("Moffat",self.psf_moff_template.amplitude.value, 
								   self.psf_moff_template.x_0.value, 
								   self.psf_moff_template.gamma.value,
								   self.psf_moff_template.alpha.value)


		self.wgauss = QtGui.QLabel(gauss_txt)
		self.wmoff = QtGui.QLabel(moff_txt)


		for lbl in (self.wgauss, self.wmoff):
				lbl.setAlignment(QtCore.Qt.AlignCenter)


		gauss_button = QtGui.QPushButton("Fit with Gaussian")
		gauss_button.clicked.connect(self.trace_with_gauss_model)
		

		moff_button = QtGui.QPushButton("Fit with Moffat")
		moff_button.clicked.connect(self.trace_with_moffat_model)

		

		gvbox = QtGui.QGridLayout()
		gvbox.addWidget(self.wgauss, 0, 0)
		gvbox.addWidget(gauss_button, 1, 0)
		gvbox.addWidget(self.wmoff, 0, 1)
		gvbox.addWidget(moff_button, 1, 1)

		

		self.wpb.setLayout(gvbox)
		self.wpb.show()

	def trace_with_gauss_model(self):

		background_poly = models.Polynomial1D(2)
		self.fitter = fitting.LevMarLSQFitter()
		self.extraction_kernel = self.psf_gauss_template + background_poly
		self.fit_extraction_kernel = self.fitter(self.extraction_kernel, self.slice_xpix, self.slice_1)
		self.fit_line = self.fit_extraction_kernel(self.slice_xpix)

		self.fig2.clf()
		gax = self.fig2.gca()

		psf_g = gax.plot(self.slice_xpix, self.fit_extraction_kernel[0](self.slice_xpix), label="Gaussian PSF")
		poly_g = gax.plot(self.slice_xpix, self.fit_extraction_kernel[1](self.slice_xpix), label="Background")
		sum_g = gax.plot(self.slice_xpix, self.fit_line, label="Composite Kernel")
		lgda_g = gax.legend()

		self.fig2.canvas.draw()
		self.wpb.close()

	def trace_with_moffat_model(self):

		background_poly = models.Polynomial1D(2)
		self.fitter = fitting.LevMarLSQFitter()
		self.extraction_kernel = self.psf_moff_template + background_poly
		self.fit_extraction_kernel = fitter(self.extraction_kernel, self.slice_xpix, self.slice_1)
		self.fit_line = self.fit_extraction_kernel(self.slice_xpix)

		self.fig2.clf()
		max = self.fig2.gca()

		psf_g = max.plot(self.slice_xpix, self.fit_extraction_kernel[0](self.slice_xpix), label="Gaussian PSF")
		poly_g = max.plot(self.slice_xpix, self.fit_extraction_kernel[1](self.slice_xpix), label="Background")
		sum_g = max.plot(self.slice_xpix, self.fit_line, label="Composite Kernel")
		lgda_g = max.legend()

		self.fig2.canvas.draw()
		self.wpb.close()

	def trace_full_spectrum_test(self):

		trace_center_model = models.Polynomial1D(0) #we use a constant because the spectrum has already been rectified
		trace_center_model.c0 = self.fit_extraction_kernel.mean_0.value # use the parameter for center of the PSF profile
		print(trace_center_model)

		spectrum = np.zeros(self.er_width, dtype=float) #initialize our spectrum with zeros
		column_pixels = np.arange(self.er_width)
		trace_centers = trace_center_model(column_pixels) # calculate our trace centers array

		# Loop over columns
		for x in column_pixels:
		    # create the kernel for this column, using the fit trace centers
		    kernel_column = self.fitter(self.extraction_kernel,
		    							self.slice_xpix, 
		    							self.slice_1) #fit_extraction_kernel.copy()
		    
		    kernel_column.mean_0 = trace_centers[x]
		    # kernel_column.stddev_0 = fwhm_fit(x) # if accounting for a varying FWHM, uncomment this line.
		    kernel_values = kernel_column(self.slice_xpix)
		    
		    # isolate the relevant column in the spectrum and variance images
		    #variance_column = variance_image[:, x] # remember that numpy arrays are row, column
		    image_pixels = self.extractregdata[:, x]
		    
		    # calculate the kernal normalization
		    #g_x = np.ma.sum(kernel_values**2 / variance_column)
		    #if np.ma.is_masked(g_x): #this column isn't valid, so we'll skip it
		    #    continue
		    
		    # and now sum the weighted column
		    weighted_column = np.ma.divide(image_pixels * kernel_values, variance_column)
		    spectrum[x] = np.ma.sum(weighted_column) / g_x




	def trace_full_spectrum(self):

		n_bin = 100
		bin_width = self.er_width // n_bin
		bin_centers = np.arange(0, self.er_width, bin_width+1, dtype=float) + bin_width // 2
		binned_spectrum = np.hstack([self.extractregdata[:, i:i+bin_width+1].sum(axis=1)[:, None] 
		                                 for i in range(0, self.er_width, bin_width+1)])
		bin_fwhms = np.zeros_like(bin_centers, dtype=float)

		for y in range(bin_centers.size):
		    bin_fit = fitter(self.fit_extraction_kernel, self.slice_xpix, binned_spectrum[:, y])
		    bin_fwhms[y] = bin_fit.stddev_0.value
		    
		bin_ny, bin_nx = binned_spectrum.shape
		bin_ar = bin_nx / (3 * bin_ny)

		
		trace_center_model = models.Polynomial1D(0) #we use a constant because the spectrum has already been rectified
		trace_center_model.c0 = self.fit_extraction_kernel.mean_0.value 
		self.spectrum = np.zeros(self.er_width, dtype=float) #initialize our spectrum with zeros
		column_pixels = np.arange(self.er_width)
		trace_centers = trace_center_model(column_pixels) # calculate our trace centers array

		# Loop over columns
		for x in column_pixels:
		    # create the kernel for this column, using the fit trace centers
		    kernel_column = self.fit_extraction_kernel.copy()
		    kernel_column.mean_0 = trace_centers[x]
		    # kernel_column.stddev_0 = fwhm_fit(x) # if accounting for a varying FWHM, uncomment this line.
		    kernel_values = kernel_column(self.slice_xpix)
		    
		    # isolate the relevant column in the spectrum and variance images
		    #variance_column = variance_image[:, x] # remember that numpy arrays are row, column
		    image_pixels = self.extractregdata[:, x]
		    #self.spectrum[x] = 

	def plot_slice(self):


		if self.extractregdata is None:
			fi = self.fitsimage
			image = fi.get_image()
			image_data = image.get_data()
		else:
			image_data = self.extractregdata

		image_shape = image_data.shape
		image_h,image_w = image_shape

		self.er_width = image_w
		self.er_height = image_h

		
		# clear previous image
		self.fig2.clf()

		ax = self.fig2.add_subplot(111)
		ax.autoscale(True, tight=True)

		
		axt = self.fig1.gca()
		print(axt)
		axx1 = self.fig1.add_axes([0.25, 0.10, 0.65, 0.03], facecolor='r')
		axx2 = self.fig1.add_axes([0.25, 0.05, 0.65, 0.03], facecolor='r')

		self.scol = Slider(axx1, 'Column idx', 0, image_w, 
							valinit=image_w//2, valstep=1)
		self.swidth = Slider(axx2, 'Width', 0, image_w, 
							valinit=30, valstep=1)

		col_idx = self.scol.val
		width = self.swidth.val
		self.col_idx = col_idx 
		self.width = width

		sx, sy, sw, sh = self.make_slice(col_idx,width)

		self.slice_rectangle = Rectangle((sx, sy), 
										sw, sh, 
						facecolor='none', edgecolor='cyan', linestyle='--')

		
		axt.add_patch(self.slice_rectangle)

		x,y = self.slice_rectangle.xy
		w = self.slice_rectangle.get_width() 
		h = self.slice_rectangle.get_height() 

		
		slice_er_y, slice_er_x = np.mgrid[y:y+h, x:x+w]
		extract_slice = self.extractregdata[slice_er_y, slice_er_x]
		extract_slice_shape = extract_slice.shape
		slice_xpix = np.arange(self.er_height)#np.arange(extract_slice_shape[0])
		
		slice_1 = self.slice_coadd()#extract_slice.sum(axis=1)
		self.slice_1 = slice_1
		self.slice_xpix = slice_xpix


		lin_raw, = ax.plot(slice_xpix, slice_1, 'k-')
		ax.set_xlabel('Cross-dispersion pixel')
		ax.axes.set_ylabel('Coadded signal')

		self.scol.on_changed(self.update_slider)
		self.swidth.on_changed(self.update_slider)
		

		self.canvas2.draw()

		
		self.fig2.canvas.draw()
		self.fig1.canvas.draw()

		self.wgetimg.setVisible(True)

	def make_slice(self,col_idx,width):

		sy, sh, sw = 0, self.er_height, width
		#print(self.col_idx)
		sx = col_idx - (width // 2)
		return sx, sy, sw, sh

	def slice_coadd(self):

		half_width = self.width // 2
		to_coadd = np.arange(max(0, self.col_idx - half_width), 
						 min(self.er_width-1, self.col_idx + half_width))
		#print(to_coadd.shape)
		return self.extractregdata[:, to_coadd].sum(axis=1) / self.width

	def update_slider(self,val):

		col_idx = self.scol.val
		width = self.swidth.val
		new_slice_box = self.make_slice(col_idx,width)
		self.slice_rectangle.set_bounds(*new_slice_box)

		x,y = self.slice_rectangle.xy
		w = self.slice_rectangle.get_width() 
		h = self.slice_rectangle.get_height() 

		
		slice_er_y, slice_er_x = np.mgrid[y:y+h, x:x+w]
		extract_slice = self.extractregdata[slice_er_y, slice_er_x]
		extract_slice_shape = extract_slice.shape
		slice_xpix = np.arange(self.er_height)#np.arange(extract_slice_shape[0])
		
		self.col_idx = col_idx 
		self.width = width
		#update line plot
		self.slice_1 = self.slice_coadd()

		#lin4.set_ydata(self.slice_1)
		#self.slice_xpix = slice_xpix
		self.fig2.clf()
		axb = self.fig2.gca()
		lin_coadd, = axb.plot(self.slice_xpix, self.slice_1, 'k-')
		axb.set_xlabel('Cross-dispersion pixel')
		axb.axes.set_ylabel('Coadded signal')
		#l.set_ydata(amp*np.sin(2*np.pi*freq*t))
		self.fig1.canvas.draw()
		self.fig2.canvas.draw()


	def fit_cross_disp_slice(self):

		max_pixel = np.argmax(self.slice_1)
		fwhm = 5.

		self.moffat_profile = models.Moffat1D(amplitude=1, gamma=fwhm, x_0=max_pixel, alpha=1)
		self.gauss_profile = models.Gaussian1D(amplitude=1, mean=max_pixel, stddev=fwhm)

		self.fig2.clf()
		axb = self.fig2.gca()
		kern5 = axb.plot(self.slice_xpix, self.slice_1 / self.slice_1[max_pixel], label='Kernel Slice')
		moff5 = axb.plot(self.slice_xpix, self.moffat_profile(self.slice_xpix), label='Moffat Profile')
		gaus5 = axb.plot(self.slice_xpix, self.gauss_profile(self.slice_xpix), label='Gaussian Profile')

		axb.legend()
		self.fig2.canvas.draw()
		self.psf_gauss_template = self.gauss_profile  
		self.psf_gauss_template.amplitude = self.slice_1[max_pixel]
		self.psf_moff_template = self.moffat_profile
		self.psf_moff_template.amplitude = self.slice_1[max_pixel]
		print(self.psf_gauss_template)
		self.popup_button()
		self.wgetimg.setDisabled(False)



	def crop_image(self):

		image = self.fitsimage.get_image()
		#print(type(image))
		image_data = image.get_data()
		image_shape = image_data.shape

		draw_pts = self.drawing.get_points()
		print(self.drawing.get_points())

		xs = sorted(list(set([int(i[0]) for i in draw_pts])))
		ys = sorted(list(set([int(i[1]) for i in draw_pts])))

		print(xs,ys)

		#extraction_region = image_data[ys[0]:ys[1],xs[0]:xs[1]]
		#self.fitsimage.set_data(cropped_image)
		
		xs0 = np.max((0,xs[0]))
		xs1 = np.min((image_shape[1],xs[1]))
		ys0 = np.max((0,ys[0]))
		ys1 = np.min((image_shape[0],ys[1]))

		er_y, er_x = np.mgrid[ys0:ys1,xs0:xs1]
		extraction_region = image_data[er_y, er_x]
		er_ny, er_nx = extraction_region.shape

		aspect_ratio = er_nx / (3. * er_ny)

		er_norm = simple_norm(extraction_region, stretch='log')

		ax = self.fig1.add_subplot(111)
		ax.autoscale(True, tight=True)

		img = ax.imshow(extraction_region, cmap='gray', aspect=aspect_ratio, 
						  norm=er_norm, interpolation='none')
		
		
		#cax = self.fig1.add_axes([0.1, 0.9, 0.8, 0.05])
		cb = self.fig1.colorbar(img,orientation='vertical')
		

		#locator = mpl.ticker.MaxNLocator(nbins=3)
		#cb.locator = locator
		#cb.update_ticks()
		#self.fitsimage.set_data(cropped_image)
		#self.clear_canvas()
		self.fig1.canvas.draw()

		self.extractregdata = extraction_region
		self.wsliceplot.setVisible(True)



	def subtract_overscan(self):

		image = self.fitsimage.get_image()
		raw_ccd = self.ccd

		self.raw_ccd = raw_ccd

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

		self.woscan.setDisabled(True)
		self.undo_woscan.setVisible(True)



	def undo_subtract_overscan(self):

		self.fitsimage.set_data(self.raw_ccd.data)
		self.ccd = self.raw_ccd

		self.undo_woscan.setVisible(False)
		self.woscan.setDisabled(False)


	#def show_headers(self):



	def find_matching_obj_arc_slit_pair(self):


		draw_pts = self.drawing.get_points()
		print(self.drawing.get_points())

		xs = sorted(list(set([int(i[0]) for i in draw_pts])))
		ys = sorted(list(set([int(i[1]) for i in draw_pts])))

		print(xs,ys)

		#extraction_region = image_data[ys[0]:ys[1],xs[0]:xs[1]]
		#self.fitsimage.set_data(cropped_image)
		
		xs0 = np.max((0,xs[0]))
		xs1 = np.min((image_shape[1],xs[1]))
		ys0 = np.max((0,ys[0]))
		ys1 = np.min((image_shape[0],ys[1]))

		er_y, er_x = np.mgrid[ys0:ys1,xs0:xs1]
		extraction_region = image_data[er_y, er_x]
		er_ny, er_nx = extraction_region.shape



def main(options, args):

	app = QtGui.QApplication(args)

	logger = log.get_logger(name="example3", options=options)
	w = FitsViewer(logger)
	w.resize(1024, 540)
	w.show()
	app.setActiveWindow(w)
	w.raise_()
	w.activateWindow()

	if len(args) > 0:
		w.load_file(args[0])

	app.exec_()


if __name__ == "__main__":

	# Parse command line options
	from argparse import ArgumentParser

	argprs = ArgumentParser()

	argprs.add_argument("--debug", dest="debug", default=False,
						action="store_true",
						help="Enter the pdb debugger on main()")
	argprs.add_argument("--profile", dest="profile", action="store_true",
						default=False,
						help="Run the profiler on main()")
	log.addlogopts(argprs)

	(options, args) = argprs.parse_known_args(sys.argv[1:])

	# Are we debugging this?
	if options.debug:
		import pdb

		pdb.run('main(options, args)')

	# Are we profiling this?
	elif options.profile:
		import profile

		print(("%s profile:" % sys.argv[0]))
		profile.run('main(options, args)')

	else:
		main(options, args)

# END