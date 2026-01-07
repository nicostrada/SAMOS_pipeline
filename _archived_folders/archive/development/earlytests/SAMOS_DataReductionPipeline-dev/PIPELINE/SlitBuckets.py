from __future__ import print_function
import sys
import os
import json
from astropy.io import fits
from astropy import units
import glob
import pandas as pd
import numpy as np
import ccdproc
from ccdproc import CCDData
import logging
import argparse
from astropy.io.fits.verify import VerifyError
from ccdproc import ImageFileCollection
from .SAMOS_mods import cutout_slit,read_fits,write_fits
import PIPELINE.slit_tracing_mods as STM 
from PIPELINE.SAMOSHelpers import MakeThumbnail


log = logging.getLogger(__name__)

class SlitBuckets:
    """
    This class is to organize the slit images from a
    comp/targ bucket and pair each slit with its
    respective comparison lamp and object image
    __init__ arg takes an image bucket that has
    already been processed by ImageProcessor
    """

    def __init__(self,imgproc_night_bucket):

        
        self.raw_data_dir = imgproc_night_bucket.raw_data_dir
        self.processing_dir = imgproc_night_bucket.processing_dir
        self.work_dir = imgproc_night_bucket.work_dir
        self.mask_proc_dir = None
        self.gain = imgproc_night_bucket.gain
        self.rdnoise = imgproc_night_bucket.rdnoise
        self.ccdsum = imgproc_night_bucket.ccdsum
        self.empty_bucket = True
        self.slit_num = None
        self.slit_targs = None
        self.slit_comps = None
        self.slit_mask = None
        self.slit_verts = imgproc_night_bucket.mf_slit_verts
        self.slit_polys = imgproc_night_bucket.mf_slit_polys
        self.slit_edges_x = imgproc_night_bucket.slit_edges_x
        self.x_edges_all_left = imgproc_night_bucket.x_edges_all_left
        self.x_edges_all_right = imgproc_night_bucket.x_edges_all_right
        self.y_edges_all_left = imgproc_night_bucket.y_edges_all_left
        self.y_edges_all_right = imgproc_night_bucket.y_edges_all_right

        self.jpeg_dir = imgproc_night_bucket.jpeg_dir
        self.slit_edges_y = imgproc_night_bucket.slit_edges_y
        self.spec_buckets = imgproc_night_bucket.spec_buckets
        self.in_prefix = imgproc_night_bucket.out_prefix
        self.combined_bucket = imgproc_night_bucket.combined_bucket
        self.out_prefix_list = None
        #each new slit ccd file will have sNNN_prefix_filename.fits, where N
        #is the slit number
    def __call__(self):


        log.debug("Working on slit cutouts.")

        if self.combined_bucket is not None:
            #combined_bucket is a single dataframe, not a list of dataframes
            image_list = self.combined_bucket['file'].tolist()
            sample_img = os.path.join(self.processing_dir,
                                    np.random.choice(image_list))
        else:
            image_list = self.spec_buckets[0]['file'].tolist()

            sample_img = os.path.join(self.raw_data_dir,
                                    np.random.choice(image_list))

        sample_ccd = read_fits(sample_img)
        self.slit_mask =  sample_ccd.header['SLIT']
        self.mask_proc_dir = os.path.join(self.processing_dir,"{:s}".format(self.slit_mask))
        if not os.path.exists(self.mask_proc_dir):
            log.debug("making slitmask directory {:s}".format(self.mask_proc_dir))
            os.mkdir(self.mask_proc_dir)
        elif len(os.listdir(self.mask_proc_dir))>1:
            log.debug("cleaning slitmask directory of previous files {:s}".format(self.mask_proc_dir))
            [os.remove(os.path.join(self.mask_proc_dir,i)) for i in os.listdir(self.mask_proc_dir)]

        for fname in image_list:
            if self.combined_bucket is not None:
                ccdpath = os.path.join(self.processing_dir,fname)

            else:
                ccdpath = os.path.join(self.processing_dir,"{:s}{:s}".format(self.in_prefix,fname))

            ccd = read_fits(ccdpath)

            
            self.slit_num = 1

            num_slits = len(self.x_edges_all_left)

            for nslit in range(num_slits):

                slit,trimsecs = STM.cutout_slit(ccd,nslit,self.x_edges_all_left,
                                                     self.x_edges_all_right,
                                                     self.y_edges_all_left,
                                                     self.y_edges_all_right)

                if self.combined_bucket is not None:
                    out_prefix = "s{:s}_".format(str(self.slit_num))
                else:
                    out_prefix = "s{:s}_".format(str(self.slit_num))+self.in_prefix
                
                slit_data = np.rot90(slit,k=-1)
                slit_ccd = CCDData(slit_data,unit='adu')
                new_slit_fname = os.path.join(self.mask_proc_dir,"{:s}{:s}".format(out_prefix,fname))

                slit_ccd.header = ccd.header.copy()
                slit_ccd.header.set('SLITNUM',
                                    value="{:s}".format(str(self.slit_num).zfill(2)),
                                    comment='')

                slit_ccd.header.set("GRP_ORDR",value=3)
                slit_ccd.header.set("GRP_FUNC",value='Chebyshev1D')
                #pipeline will have to add these keywords for wcs
                #according to user input
                #using default for now.


                slit_ccd.header.add_comment('Data rotated 90 degrees to match how Goodman treats dispersion axis.')
                slit_ccd.write(new_slit_fname,overwrite=True)

                MakeThumbnail(new_slit_fname,self.jpeg_dir)

                self.slit_num+=1

                if self.out_prefix_list is None:
                    self.out_prefix_list = [out_prefix]
                else:
                    self.out_prefix_list.append(out_prefix)

                if slit_ccd.header['OBSTYPE']=='OBJECT':
                    log.debug("appending slit target")
                    self.add_slit_targ(new_slit_fname)

                elif slit_ccd.header['OBSTYPE']=='COMP':
                    log.debug("appending slit comparison lamp")
                    self.add_slit_comp(new_slit_fname)



    def add_slit_comp(self,slit_comp):

        if self.slit_comps is None:
            self.slit_comps = [slit_comp]
        else:
            self.slit_comps.append(slit_comp)

    def add_slit_targ(self,slit_targ):

        if self.slit_targs is None:
            self.slit_targs = [slit_targ]
        else:
            self.slit_targs.append(slit_targ)
