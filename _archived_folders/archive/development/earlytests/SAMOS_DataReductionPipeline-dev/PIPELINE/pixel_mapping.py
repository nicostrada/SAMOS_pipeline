import numpy as np
#import pyfits

from astropy.io import fits as pyfits

from astroquery.simbad import Simbad
import astropy.units as u
import astropy.coordinates as coord

from astroquery.vizier import Vizier

import os, sys
import glob, math
from math import sin, cos, tan, asin, sqrt, pi
import time, datetime
import functools
import string
import urllib


def create_points_array(nx, ny):

    """

    Input shape of original array (nx,ny). 
    Outputs matrix with row vectors of points in the format (row, col, index)


    """

    points = []

    ind = -1
    for row in range(nx):
        for col in range(ny):
            ind += 1
            point = (row,col,ind)
            points.append(point)


    A = np.mat(np.array([p for p in points]))

    return A


class DMD:
        def __init__(self, dmd_arr):
            self.nx_mirrors = self.ny_mirrors = int(1080) # DMD field of view

            self.points_mat = create_points_array(self.nx_mirrors,self.ny_mirrors)

            self.dmd_scale = 1./6. # arcsec per mirror

            self.scale_mat = get_scale_matrix(self.dmd_scale,self.dmd_scale)

            self.dmd_arr = dmd_arr

    
class ImCam:

    def __init__(self,im_arr):

        self.nx_pix = self.ny_pix = 1024 

        self.ic_scale = 180./1024. # arcsec/pixel (3 arcmin FOV)

        self.points_mat = create_points_array(self.nx_pix,self.ny_pix)

        self.scale_mat = get_scale_matrix(self.ic_scale,self.ic_scale)

        self.im_arr = im_arr




def do_transform(t_mat, input_arr, new_arr_shape):

    new_nrows, new_ncols = new_arr_shape
    arr_transformed = np.empty((new_nrows, new_ncols), dtype=np.uint8)

    for i, row in enumerate(input_arr):
        for j, col in enumerate(row):
            arr_data = input_arr[i, j]
            input_coords = np.array([i, j, 1])
            i_out, j_out = t_mat @ input_coords
            arr_transformed[i_out, j_out] = arr_data
    
    return arr_transformed



def get_rotation_matrix(ang_deg):

    ang_rad = ang_deg * math.pi/180.

    sin_ang = sin(ang_rad)
    cos_ang = cos(ang_rad)

    return np.mat(np.array([[cos_ang, sin_ang, 0], 
                            [-sin_ang, cos_ang, 0],
                            [0, 0, 1]]))


def get_scale_matrix(scale_x,scale_y):

    return np.mat(np.array([[scale_x, 0, 0], 
                            [0, scale_y, 0], [0, 0, 1]],dtype=float))






