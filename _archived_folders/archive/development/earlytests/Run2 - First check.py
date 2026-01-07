#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 07:02:45 2024

@author: robberto
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from skimage.measure import label#, regionprops

import os
import glob
import numpy as np
import pandas as pd
import copy
import heapq


""" A Routine to read and rearrange correctly the 4 SAMI CCDs """
def read_SAMI_mosaic(file):
    hdul = fits.open(file)
    hdr = hdul[0].header
    data1 = hdul[1].data ; data1=data1[:,54:2101]
    data2 = hdul[2].data ; data2=data2[:,54:2101]
    data3 = hdul[3].data ; data3=data3[:,54:2101]
    data4 = hdul[4].data ; data4=data4[:,54:2101]
    # Create a blank mosaic image
    mosaic_data = np.zeros((data1.shape[0]*2, data1.shape[1]*2))
    # Place the images into the mosaic
    mosaic_data[0:data1.shape[0], 0:data1.shape[1]] = np.flip(data3,axis=0)  #top left
    mosaic_data[data1.shape[0]:, 0:data1.shape[1]] += np.flip(data1,axis=0)   #bottom left                  
    mosaic_data[0:data1.shape[0]:, data1.shape[1]:] += np.flip(data4,axis=0)                     
    mosaic_data[data1.shape[0]:, data1.shape[1]:] += np.flip(data2,axis=0)
    hdu = fits.PrimaryHDU(data=mosaic_data, header=hdr)
    return hdu




""" DATA ANALYSIS R136A

=> Night 20241016

"""

#CREATE A WORKING DIRECTORY FIRST IF DOES NOT EXIST
working_directory = "/Users/robberto/Library/CloudStorage/Box-Box/My Documents - Massimo Robberto/@Massimo/_Science/2. Projects_HW/SAMOS/SAMOS_DATA_ANALYSIS"
os.chdir(working_directory) 
print(working_directory)

data_directory = "/Users/robberto/Library/CloudStorage/Box-Box/My Documents - Massimo Robberto/@Massimo/_Science/2. Projects_HW/SAMOS/SAMOS_DATA/RUN2/SAMI/20241016"


#FIRST STEP: CREAT A MASK FOR MAP T00
File_NR = '124'
files = []
files.append(os.path.join(data_directory,'target.'+File_NR+'.fits'))
hdu = read_SAMI_mosaic(files[0])
mask = hdu.data
mask[83:1275,2046:2056]=65535 #CLEANUP THE SPACE BETWEEN SAMI CCDS  
hdu.data = mask
hdul = fits.HDUList([hdu]) ; hdul.writeto('test.fits',overwrite=True)
plt.imshow(mask, vmin=1500, vmax=18000)

binary_mask = mask*0  # WE WANT MASKS of 0's and 1's

#Set the imaging area, no edges that are confusing. By eye!
y0 = 85
y1 = 1270
dy = y1-y0   # 1185
x0 = 1490
x1 = 2662
dx = x1-x0   # 1172

for i in range(dx):
    """ 
    plt.figure(figsize=(12, 4))
    #TO FIND THE EDGES WE SCROLL OVER ALL VERTICAL PIXELS
    row_min = np.min(mask[:,1490+i:1491+i+1], axis=1) ;
    row_min_p1 = np.roll(row_min,+1)
    diff = row_min - row_min_p1
    max_index = np.argmax(diff)
    min_index = np.argmin(diff)
    print(max_index, min_index)
    plt.plot(row_min-row_min_p1,c='black')
    plt.show()
    #AND WE FIND THAT THE EDGES AT
    #  ROWS 79 AND 1273 FOR COLUMN 1490
    #  ROWS 81 AND 1272 FOR COLUMN 1490+1172  = 2662
    """
    
    #"""
    #TO FIND THE EDGES WE GO LET TO RIGHT (ON THE X AXIS) CHECKING EACH COLUMN...
    #plt.figure(figsize=(12, 4))
    col_min = np.min(mask[y0:y1,x0+i:x0+1+i], axis=1) ;
    col_min_p1 = np.roll(col_min,+1)    
    diff = col_min - col_min_p1
    max_index = np.argmax(diff)
    min_index = np.argmin(diff)
    print(x0+i, min_index, max_index)
    #plt.plot(col_min-col_min_p1,c='black')
    #plt.show()
    
    
    if (max_index - min_index) >5: # We found a mask?
        print(x0+i,min_index, max_index)
        for j in range(len(col_min)):
            if diff[j] != 0:
                binary_mask[j+y0-1:j+y0+2,x0+i-1:x0+i+2] += 1
                binary_mask[j+y0,x0+i:x0+i+1] += 1
        print("row done")    
        
        
   #"""
""" LABEL THE SLITS (CLUSTER) 
create a copy of the mask where all the pixels in a clusters/slit share the same interger
"""
all_masks = label(binary_mask > 0)
hdu.data = all_masks
hdul = fits.HDUList([hdu]) ; hdul.writeto('/Users/robberto/Desktop/arr.fits',overwrite=True)

""" CREATION OF THE MASKS"""
#hdul = fits.open('/Users/robberto/Desktop/arr.fits')
#all_masks =   hdul[0].data
Nr_of_masks = all_masks.max()
for i_mask in range(Nr_of_masks):
    mask= copy.deepcopy(all_masks)
    mask[mask != (i_mask+1)]=0
    pixel_mask = np.argwhere(mask != 0)
#check
    for index in pixel_mask:
        row, col = index
        #print(f"Value at ({row}, {col}): {mask[row, col]}")
    x0,x1 = min(pixel_mask[:,1]),max(pixel_mask[:,1])
    xc = int(np.round((x0+x1)/2))
    y0,y1 = min(pixel_mask[:,0]),max(pixel_mask[:,0])
    #creation of the trace
    trace = copy.deepcopy(mask)
    trace[y0:y1,:]=1
    hdu.data = trace
    hdul = fits.HDUList([hdu]) ; hdul.writeto('trace_'+str(i_mask+1)+'.fits',overwrite=True)


    
print("check")
            
#EXTRACT SPECTRA

#


File_NR = '088'
files = []
files.append(os.path.join(data_directory,'target.'+File_NR+'.fits'))
hdu = read_SAMI_mosaic(files[0])
all_data = hdu.data
plt.imshow(all_data,vmin=500,vmax=700)
plt.show()


#BIAS IMAGE
File_NR = '074'
files = []
files.append(os.path.join(data_directory,'bias.'+File_NR+'.fits'))
hdu = read_SAMI_mosaic(files[0])
bias = hdu.data
plt.imshow(bias, vmin=500, vmax=800)

#FLATFIELD IMAGE
File_NR = '117'
files = []
files.append(os.path.join(data_directory,'target.'+File_NR+'.fits'))
hdu = read_SAMI_mosaic(files[0])
flat = hdu.data
plt.imshow(flat, vmin=500, vmax=800)

#ARC LAMP IMAGE
File_NR = '116'
files = []
files.append(os.path.join(data_directory,'target.'+File_NR+'.fits'))
hdu = read_SAMI_mosaic(files[0])
arc  = hdu.data
plt.imshow(arc, vmin=500, vmax=800)


#SUBTRACT BIAS FROM IMAGE
all_data_biassub = all_data-bias
plt.imshow(all_data_biassub, vmin=0, vmax=300)

#SUBTRCT BIAS FROM FLATS
flat_biassub=flat-bias
plt.imshow(flat_biassub, vmin=0, vmax=30000)
plt.show()

#SUBTRCT ARCFROM FLATS
arc_biassub=arc-bias
plt.imshow(arc_biassub, vmin=0, vmax=30000)
plt.show()


"""
PREPARE FOR FOR THE FINAL PLOT
"""
#fix, axes = plt.subplots(nrows=10, ncols=1, figsize = (12,20))


#EXTRACT FRAMES WITH SINGLE SPECTRA
for i_mask in range(Nr_of_masks):
#    plt.imshow(all_data_biassub,vmin=0,vmax=200)
#    plt.show()
    hdul = fits.open('trace_'+str(i_mask+1)+'.fits')
    trace =   hdul[0].data
    data_masked = all_data_biassub*trace
#    z0 = 0#data_masked.min()
#    z1 = 200#data_masked.max()
#    plt.imshow(data_masked, vmin=z0, vmax=z1)
#    plt.show()
    hdu.data = data_masked
    hdul = fits.HDUList([hdu]) ; hdul.writeto('data_masked_'+str(i_mask+1)+'.fits',overwrite=True)
    print('done mask',i_mask)
    
    flat_masked = flat_biassub * trace
#    z0 = 0#data_masked.min()
#    z1 = 20000#data_masked.max()
#    plt.imshow(flat_masked, vmin=z0, vmax=z1)
#    plt.show()
    hdu.data = flat_masked
    hdul = fits.HDUList([hdu]) ; hdul.writeto('flat_masked_'+str(i_mask+1)+'.fits',overwrite=True)
    print('done flat mask',i_mask)
    
    pixel_mask = np.argwhere(trace != 0)
    row=[]
    col=[]
    for index in pixel_mask:
        row.append(index[0])
        col.append(index[1])
        #print(f"Value at ({row}, {col}): {mask[row, col]}")
    spec2d = data_masked[min(row):max(row),min(col):max(col)]   
    flat2d = flat_masked[min(row):max(row),min(col):max(col)]   
    flat2d_mean = np.mean(flat2d)
    flat2d_norm = flat2d/np.mean(flat2d)
    print('flat normalized')
    
    spec2d_ff = spec2d / flat2d_norm
    hdu.data = spec2d_ff
    hdul = fits.HDUList([hdu]) ; hdul.writeto('spec2d_ff_'+str(i_mask+1)+'.fits',overwrite=True)
    print('spectrum flatfielded', i_mask)
    #plt.imshow(spec2d_ff, vmin=0, vmax=30000)
    #plt.show()
    spec1d_ff = np.sum(spec2d_ff, axis=0)
    plt.ylim(0, np.median(spec1d_ff)*8)  # Sets the y-axis range from 0 to 1
#    plt.plot(spec1d_ff)
#    plt.show()
    
#   ARC TRACES
    arc_masked  = arc_biassub * trace
    plt.ylim(0, 100)
    plt.plot(arc_masked)
    plt.show()
    hdu.data = arc_masked
    hdul = fits.HDUList([hdu]) ; hdul.writeto('arc_masked_'+str(i_mask+1)+'.fits',overwrite=True)
    print('done arc masked',i_mask)
    arc2d = flat_masked[min(row):max(row),min(col):max(col)]   
    arc2d_ff = arc2d / flat2d_norm
    hdu.data = arc2d_ff
    hdul = fits.HDUList([hdu]) ; hdul.writeto('arc2d_ff_'+str(i_mask+1)+'.fits',overwrite=True)
    print('spectrum flatfielded', i_mask)


"""
    print('done spec2d_ff')
    #plt.figure(figsize=(12, 3))
    #all plots in pages
    counter=i_mask % 10
    plot_nr = int(i_mask/10)
    print(counter,plot_nr)
    plt.axes[counter].plot(spec1d_ff)
    plt.axes[counter].set_title(f'Slit {i_mask}' )
    ymax = min((heapq.nlargest(50, spec1d_ff)))*2
    plt. axes[counter].axis(ymin=0, ymax=ymax)#np.median(spec1d_ff)*3)  # Sets the y-axis range from 0 to 1
    #plt.subplot(10,1,counter+1)
    #plt.plot(spec1d_ff)
    if counter == 9:
        plt.tight_layout()
        plt.savefig('spectra_ff_'+str(plot_nr)+'.')        
        plt.show()
        fix, axes = plt.subplots(nrows=10, ncols=1, figsize = (12,20))

    

# Calculate the average of each row
row_averages = np.mean(flat_biassub[:,1000:1200], axis=1)

plt.figure(figsize=(12, 4))
plt.plot(row_averages[200:300])
"""

