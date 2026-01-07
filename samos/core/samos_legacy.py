#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 07:02:45 2024

@author: robberto
"""

from astropy.io import fits
#import matplotlib.pyplot as plt

import os
#import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style

import lacosmic

#import pandas as pd
#import copy
#import heapq

class SAMOS:
    def __init__(self, dir_name):
        self.dir_name = dir_name
        print("working on the directory:", dir_name)
        #self.file = Path.cwd() / file_name
        #self.db_dict = self._load_db()

    """ A Routine to read and rearrange correctly the 4 SAMI CCDs """
    def read_SAMI_mosaic(self,file):
        hdul = fits.open(file)
        hdr = hdul[0].header
        data1 = hdul[1].data ; data1=data1[:,54:2101]  #top left, flipped y
    #    plt.imshow(data1, vmin=500, vmax=3000)
    #    plt.show()
        data2 = hdul[2].data ; data2=data2[:,64:2111]  #top right, flipped y 
    #    plt.imshow(data2, vmin=500, vmax=3000)
    #    plt.show()
        data3 = hdul[3].data ; data3=data3[:,54:2101]  #bottom left, flipped y
    #    plt.imshow(data3, vmin=500, vmax=3000)
    #    plt.show()
        data4 = hdul[4].data ; data4=data4[:,64:2111]  #bottom right, flipped y
    #    plt.imshow(data4, vmin=500, vmax=3000)
    #    plt.show()
    #    print(data1.shape)
        # Create a blank mosaic image
        mosaic_data = np.zeros((data1.shape[0]*2, data1.shape[1]+data2.shape[1]))
        # Place the images into the mosaic
        mosaic_data[data1.shape[0]:, 0:data1.shape[1]] += np.flip(data1,axis=0)  #  top left                  
        mosaic_data[data1.shape[0]:, data1.shape[1]:] += np.flip(data2,axis=0)   #  top right
        mosaic_data[0:data1.shape[0], 0:data3.shape[1]] = np.flip(data3,axis=0)  #  bottom left
        mosaic_data[0:data1.shape[0]:, data3.shape[1]:] += np.flip(data4,axis=0) #  bottom right                    
        mosaic = np.flip(mosaic_data,axis=0)
    #    plt.imshow(mosaic, vmin=500, vmax=3000)
    #    plt.show()
    
        """fix line 2046"""
        mosaic[:,2046] = (mosaic[:,2045]+mosaic[:,2047])/2
        
    
        hdu = fits.PrimaryHDU(data=mosaic, header=hdr)
        return hdu

        
    def display_image(self, image,zmin,zmax):
        # Set the visualization style
        plt.figure(figsize=(14,8))
        plt.style.use(astropy_mpl_style)
        # Display the image
        plt.imshow(image, origin='lower', cmap='gray',vmin=zmin, vmax=zmax)
        #the next two lines scale the color bar to the height of the figure
        im_ratio = image.shape[0]/image.shape[1]
        plt.colorbar(fraction=0.046*im_ratio, pad=0.04)
        plt.show()
        

# Load your image data
    def CR_correct(self,image):
        
        # Remove cosmic rays
        cleaned_data, mask = lacosmic.lacosmic(image, contrast=5, cr_threshold=20, neighbor_threshold=5, readnoise=3, effective_gain=2)
        
        # Save the cleaned image
        #hdu[0].data = cleaned_data
        #hdu.writeto("cleaned_image.fits", overwrite=True)
        return cleaned_data