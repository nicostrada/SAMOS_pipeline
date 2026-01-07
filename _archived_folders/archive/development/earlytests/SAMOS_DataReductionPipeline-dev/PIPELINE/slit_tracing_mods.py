import sys
import os
import json
import numpy as np
import re
import pickle
from astropy.io import fits
from astropy import units
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.coordinates import EarthLocation
from astropy.modeling import (models, fitting, Model)
import datetime
from astroscrappy import detect_cosmics
from astropy.stats import sigma_clip
import subprocess
import time
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.wcs import WCS

import glob
import pandas as pd
import ccdproc
from ccdproc import CCDData,ImageFileCollection
import logging
import argparse
from .DataBucket import DataBucket
from .SAMOS_NIGHT import SAMOSNight,GenerateDcrParFile
from .SAMOSHelpers import is_file_saturated
import warnings
from threading import Timer
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Qt5Agg')
import scipy
from scipy import signal
warnings.filterwarnings('ignore')



### Pixel-World header ###


pixscale = 0.14814815 
pixscale_deg = pixscale/3600.
decpangle = 0.
cd11 = pixscale_deg*np.cos(decpangle)
cd22 = cd11
cd12 = pixscale_deg*np.sin(decpangle)
cd21 = -pixscale_deg*np.sin(decpangle)

fake_SAMI_header_img = { 'NAXIS' : 2,
						  'NAXIS1' : 1215,
						  'NAXIS2' : 1215,
						  'PIXSCALE' : pixscale,
						  'CD1_1' : cd11,
						  'CD1_2' : cd12,
						  'CD2_1' : cd21,
						  'CD2_2' : cd22,
						  'CTYPE1' : 'RA---TAN',
						  'CTYPE2' : 'DEC--TAN',
						  'CRVAL1' : 92.0215197, #reference RA (deg)
						  'CRVAL2' : 24.3394499, #reference dec (deg)
						  'CRPIX1' : 654, #reference X val (pixel)
						  'CRPIX2' : 433  #reference Y val (pixel)
						}


### DMD-World header ###


mirror_scale = 1./6. #arcsec per mirror
mscale_deg = mirror_scale/3600.
#pretend decpangle (=PAngle) is the same
m_cd11 = mscale_deg*np.cos(decpangle)
m_cd22 = m_cd11
m_cd12 = mscale_deg*np.sin(decpangle)
m_cd21 = -mscale_deg*np.sin(decpangle)



wrld_to_DMD_header = {'NAXIS' : 2,
					  'NAXIS1' : 1080,
					  'NAXIS2' : 1080,
					  'PIXSCALE' : mirror_scale, 
					  'CD1_1' : m_cd11,
					  'CD1_2' : m_cd12,
					  'CD2_1' : m_cd21,
					  'CD2_2' : m_cd22,
					  'CTYPE1' : 'RA---TAN',
					  'CTYPE2' : 'DEC--TAN',
					  'CRVAL1' : 92.0215197, #reference RA (deg)
					  'CRVAL2' : 24.3394499, #reference dec (deg)
					  'CRPIX1' : 300, #reference X val (pixel)
					  'CRPIX2' : 200  #reference Y val (pixel)
						}


def sky_coords_to_DMD(skycoords,dmd_header=wrld_to_DMD_header):


	wcs = WCS(dmd_header)

	return wcs.world_to_array_index(skycoords)

def DMD_to_sky(dmd_col,dmd_row,dmd_header=wrld_to_DMD_header):

	wcs = WCS(dmd_header)

	return wcs.array_index_to_world(dmd_row,dmd_col)


def sky_coords_to_SAMI(skycoords,sami_header=fake_SAMI_header_img):

	wcs = WCS(sami_header)

	return wcs.world_to_array_index(skycoords)


def DMD_to_SAMI(dmd_col,dmd_row,dmd_header=wrld_to_DMD_header,
				sami_header=fake_SAMI_header_img):


	sky_from_dmd = DMD_to_sky(dmd_col,dmd_row,dmd_header)

	sami_from_sky = sky_coords_to_SAMI(sky_from_dmd,sami_header)

	print("DMD (row,col) = %s,%s"%(dmd_row,dmd_col)) 
	print("RA,DEC = %s,%s"%(sky_from_dmd.ra.deg,sky_from_dmd.dec.deg))
	print("SAMI (row,col) = %s,%s"%(sami_from_sky[0],sami_from_sky[1]))

	return sami_from_sky


def read_initial_slit_edges_from_header(img):

	slit_df_data =[]
	
	for slit_def in list(img.header["slit_**"].values()):

		
		specsec_xcols = int(img.header['specsec'].split(",")[0].split(":")[0].strip("'["))
		specsec_yrows = int(img.header['specsec'].split(",")[1].split(":")[0])

		cutout_size = 10
		binsize = 30



		x_left_start = int(slit_def.split(',')[0].split(':')[0].strip("'").strip("['").strip("]'"))-specsec_xcols
		x_right_start = int(slit_def.split(',')[0].split(':')[1].strip("'").strip("['").strip("]'"))-specsec_xcols
		y_mid_start = int(slit_def.split(',')[1].strip("'").strip("['").strip("]'"))-specsec_yrows


		slit_df_data.append([x_left_start,x_right_start,y_mid_start])

	header_slit_edges = pd.DataFrame(slit_df_data,
		columns=["x_edge_left","x_edge_right","y_mid"])

	return header_slit_edges


def identify_slits(master_flat,x_init,y_init,binsize=30,cutout_size=10):

	previous_y_coord = y_init
	previous_x_coord = x_init

	x_edges = []
	y_edges = []
	
	
	N_pixel_x = master_flat.shape[1]
	N_pixel_y = master_flat.shape[0]
	while previous_y_coord < N_pixel_y:

		end_pixel = np.asarray([previous_y_coord + binsize-1, N_pixel_y-1]).min()

		cutout_bin = master_flat[previous_y_coord:end_pixel,
						 np.int(previous_x_coord-cutout_size/2.):
							np.int(previous_x_coord+cutout_size/2.)]

		#print(np.median(cutout_bin, axis=0).shape)
		#break
		if cutout_bin.shape[0]>0:

			profile = np.median(cutout_bin, axis=0)
			derivative = np.abs(np.roll(profile,1)-profile)
			derivative[0] = 0
			derivative[-1] = 0
			peak,peak_location = derivative.max(),np.argmax(derivative)
			#print(derivative)

			peak_location += (previous_x_coord - (cutout_size/2.))
			if peak_location < 0:
				peak_location = 0
			if peak_location > N_pixel_x-1:
				peak_location = N_pixel_x-1

			#print("peak location updated: %s"%(peak_location))

			y_edges.append(int(np.round(0.5*(previous_y_coord + end_pixel),0)))
			x_edges.append(int(peak_location))

			previous_x_coord = peak_location

		previous_y_coord += binsize

	previous_y_coord = y_init
	previous_x_coord = x_init


	while previous_y_coord > 0:

		end_pixel = np.asarray([previous_y_coord - binsize-1, 15]).max()

		cutout_bin = master_flat[end_pixel:previous_y_coord,
						 np.int(previous_x_coord-cutout_size/2.):
							np.int(previous_x_coord+cutout_size/2.)]

		#break
		if cutout_bin.shape[0]>0:

			profile = np.median(cutout_bin, axis=0)
			derivative = np.abs(np.roll(profile,1)-profile)
			#print(cutout_bin.shape)
			derivative[0] = 0
			derivative[-1] = 0
			peak,peak_location = derivative.max(),np.argmax(derivative)
			#print(peak)

			peak_location += (previous_x_coord - (cutout_size/2.))
			if peak_location < 0:
				peak_location = 0
			if peak_location > N_pixel_x-1:
				peak_location = N_pixel_x-1

			#print("peak location updated: %s"%(peak_location))

			y_edges.insert(0,int(np.round(0.5*(previous_y_coord + end_pixel),0)))
			x_edges.insert(0,int(peak_location))
			#print(previous_ycoord)

			previous_x_coord = peak_location
			#print(peak_location,int(np.round(0.5*(previous_y_coord + end_pixel),0)))

		previous_y_coord -= binsize
		
	return x_edges,y_edges


def cutout_slit(ccd,nslit,x_edges_left,x_edges_right,
				 y_edges_left,y_edges_right):

	first_row,last_row = y_edges_left[nslit][0],y_edges_left[nslit][-1]
	rows = last_row-first_row

	first_col,last_col = x_edges_left[nslit][0],x_edges_right[nslit][0]
	cols = last_col-first_col


	last_row_first_col,last_row_last_col = x_edges_left[nslit][-1],x_edges_right[nslit][-1]
	
	last_row_cols = last_row_last_col-last_row_first_col


	blank_slit = None


	for row_chunk in range(len(y_edges_left[nslit])-1):
		
		slit_chunk = ccd[y_edges_left[nslit][row_chunk]:y_edges_left[nslit][row_chunk+1],
						x_edges_left[nslit][row_chunk]:x_edges_right[nslit][row_chunk]]
		#print(slit_chunk.shape)
		
		
		if blank_slit is None:
			blank_slit = np.array(slit_chunk)
			cols_in_chunks = slit_chunk.shape[1]
			
		if slit_chunk.shape[1]>cols_in_chunks:
			slit_chunk = slit_chunk[:,:cols_in_chunks]
			#print("new chunk shape",slit_chunk.shape)

			
		elif slit_chunk.shape[1]<cols_in_chunks:
			
			num_zero_cols = cols_in_chunks-slit_chunk.shape[1]
			slit_chunk = np.hstack((slit_chunk,np.zeros((slit_chunk.shape[0],num_zero_cols))))

			#print("new chunk shape",slit_chunk.shape)

		blank_slit = np.vstack((blank_slit,slit_chunk))

	fits_slit_trim_sections = '[{:d}:{:d},{:d}:{:d}]'.format(\
							int(first_col),int(last_col),\
							int(first_row),int(last_row))
							#top of slit y and bottom of slit y

	

	#log.debug("N Slits Trimmed: {:s}".format(str(len(fits_slit_trim_sections))))

	return blank_slit,fits_slit_trim_sections



def cutout_slit_old(ccd,x_edges,y_edges):

	#x_edges, y_edges = get_edges(input,'LMask2_ycoords_c1.txt','LMask2')

	margin = 5 #extra space to include in cutout_bin
	slit_height = 50#found by inspection)
	sz = ccd.shape
	N_pixel_x = sz[1]
	N_pixel_y = sz[0]
	print("number of cols: %f"%N_pixel_x)
	print("number of rows: %f"%N_pixel_y)
	x_axis = np.linspace(0,N_pixel_x-1,num=N_pixel_x,dtype=int)
	y_axis = np.linspace(0,N_pixel_y-1,num=N_pixel_y,dtype=int)
	pixel_x = np.zeros((N_pixel_y,N_pixel_x),dtype=int)
	for row in range(N_pixel_y):
		for col in range(N_pixel_x):
			pixel_x[row,col] = col

	pixel_y = np.zeros((N_pixel_y,N_pixel_x),dtype=int)
	for row in range(N_pixel_y):
		for col in range(N_pixel_x):
			pixel_y[row,col] = row


	cutout_slit_arrays_full = []
	new_y_edges = np.empty_like(y_edges)

	new_top_ys = np.linspace(slit_height,ccd.shape[0],num=y_edges.shape[0],dtype=int)

	for row in range(y_edges.shape[0]):
		new_y_edges[row].fill(new_top_ys[row])
	#print(new_y_edges)

	#make a copy of the full ccd for each slit and
	#and populate new ccd arrays with data where the slit
	#is located, and nans everywhere else.
	ccd_slice = np.asarray(ccd[:,:])
	for slit in range(y_edges.shape[0]):
		im = np.copy(ccd_slice)#np.asarray(ccd_slice.copy())
		im[:,:] = np.nan
		top_y = y_edges[slit]
		top_y = np.asarray(top_y)
		bottom_y = top_y-slit_height
		bottom_y = np.asarray(bottom_y).astype(int)
		#in the new image array, I want to shift the columns of included data
		#so they all align in the same range of rows. That way the slit cutouts
		#won't be all slanted.

		new_top_y = np.asarray(new_y_edges[slit])

		new_bot_y = new_top_y-slit_height
		new_bot_y = np.asarray(new_bot_y).astype(int)
		for col in range(im.shape[1]):
			col = int(col)

			im[int(new_bot_y[col]):int(new_top_y[col]),col] = \
				   ccd_slice[int(bottom_y[col]):int(top_y[col]),col]

		slit_cutout_im = im[~np.isnan(im).all(axis=1)]
		cutout_slit_arrays_full.append(slit_cutout_im)
		#final slit cutout array is trimmed of nans

	fits_slit_trim_sections = ['[1:{:d},{:d}:{:d}]'.format(int(x_edges[i][-1]),
							int(y_edges[i][0]),int(y_edges[i][0])+slit_height) for i in \
							range(len(x_edges))]
													#top of slit y and bottom of slit y
	log.debug("N Slits Trimmed: {:s}".format(str(len(fits_slit_trim_sections))))

	return cutout_slit_arrays_full,fits_slit_trim_sections
