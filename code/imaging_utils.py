# This is file contains the program that imaging point source in ALMACAL project
#
# Author: Jianhang Chen
# Email: cjhastro@gmail.com
# LastUpdate: 2020-11-05

import os
import re
import json
import numpy as np
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt 
import analysisUtils as au

def read_spw(vis):
    """read the spectral windows
    """
    tb = tbtool()
    tb.open(vis + '/SPECTRAL_WINDOW')
    col_names = tb.getvarcol('NAME')
    col_freq = tb.getvarcol('CHAN_FREQ')
    tb.close()

    spw_specrange = {}

    for key in col_names.keys():
        freq_max = np.max(col_freq[key]) / 1e9
        freq_min = np.min(col_freq[key]) / 1e9
        spw_specrange[key] = [freq_min, freq_max]

    return spw_specrange.values()

def efficient_imsize(imsize):
    """This function try to optimize the imsize that can be divided by 2,3,5,7 only
    Which is much faster for casa
    """
    while True:
        if ((imsize%3 == 0) | (imsize%5 == 0) | (imsize%7 == 0)) & (imsize%2 == 0):
            return imsize
        else:
            print(imsize)
            imsize += 1
            continue
        

def make_cont_img(vis=None, basename=None, dirty_image=False, clean_image=False, myimagename=None, mycell=None, myimsize=None, outdir='./', dry_run=False, **kwargs):
    """This function is used to make the continuum image
    
    dirty_image: generate the dirty image

    clean_image: make the clean image with tclean
        You need to specify: weighting, niter, interactive
        The imagename, imsize, cell and restfreq will be set automatically if no value is provided by kwargs
    """
    spw_specrange = read_spw(vis)
    freq_mean = np.mean(spw_specrange) # in GHz
    
    # read the antenna_diameter
    tb = tbtool()
    tb.open(vis + '/ANTENNA')
    antenna_diameter = np.median(tb.getcol('DISH_DIAMETER'))
    tb.close()


    wavelength = (const.c / (freq_mean * u.GHz)).to(u.um) # in um
    fov = 2.0 * 1.22 * wavelength / (antenna_diameter*u.m).to(u.um) * 206265

    # get the baselines
    # baselines = au.getBaselineLengths(vis)

    # another way is to use getBaselineStats directly
    baseline_stat = au.getBaselineStats(vis)
    baseline_max = baseline_stat[2] * u.m
    baseline_90 = baseline_stat[-1] * u.m

    baseline_typical = baseline_90

    # if baseline_stat[0] < 300:
        # baseline_typical = baseline_max
    # else:
        # baseline_typical = baseline_90
    # calcuate the cell size
    if mycell is None:
        cellsize = 206265 / (baseline_typical / wavelength).decompose() / 6
        mycell = "{:.4f}arcsec".format(cellsize)
    
    if myimsize is None:
        myimsize = efficient_imsize(int(fov / cellsize))
        myrestfreq = str(freq_mean)+'GHz'


    if basename is None:
        basename = os.path.basename(vis)
    if myimagename is None:
        myimagename = os.path.join(outdir, basename + '.cont.auto')
    rmtables(tablenames=myimagename + '.*')

    if dry_run:
        print("Mean frequecy:", freq_mean)
        print("Maximum baseline:", baseline_max)
        print("My cell size:", mycell)
        print("My image size:", myimsize)
        return

    if dirty_image:
        tclean(vis=vis, spw="",
               datacolumn="data",
               imagename=myimagename,
               imsize=myimsize, cell=mycell, 
               restfreq=myrestfreq, phasecenter="", 
               specmode="mfs", outframe="LSRK",
               pblimit=0.1,
               weighting='natural',
               # weighting="briggs", robust=1.5,
               niter=0,
               # usemask='user',
               # usemask="auto-multithresh",
               interactive=False,
               savemodel="none", **kwargs)

    if clean_image:
        tclean(vis=vis, spw="",
               datacolumn="data",
               imagename=myimagename,
               imsize=myimsize, cell=mycell, 
               restfreq=myrestfreq, phasecenter="", 
               specmode="mfs", outframe="LSRK",
               pblimit=0.1,
               **kwargs)

