# This is file contains the program that imaging interferometer data automatically
#
# Author: Jianhang Chen
# Email: cjhastro@gmail.com

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

    spw_specrange = {}

    for key in col_names.keys():
        freq_max = np.max(col_freq[key]) / 1e9
        freq_min = np.min(col_freq[key]) / 1e9
        spw_specrange[key] = [freq_min, freq_max]

    return spw_specrange.values()

def effective_imsize(imsize):
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
        



def make_cont_img(vis=None, basename=None, antenna_diameter=12, dirty_image=False, debug=False,  **kwargs):
    """This function is used to make the continuum image

    """
    spw_specrange = read_spw(vis)
    freq_mean = np.mean(spw_specrange) # in GHz
    
    wavelength = (const.c / (freq_mean * u.GHz)).to(u.um) # in um
    fov = 2.0 * 1.22 * wavelength / (antenna_diameter*u.m).to(u.um) * 206265

    # get the baselines
    # baselines = au.getBaselineLengths(vis)
    baseline_stat = au.getBaselineStats(vis)
    baseline_max = baseline_stat[2] * u.m
    baseline_90 = baseline_stat[-1] * u.m

    # calcuate the cell size
    cellsize = 206265 / (baseline_90 / wavelength).decompose() / 6
    mycell = "{:.4f}arcsec".format(cellsize)
    
    myimsize = effective_imsize(int(fov / cellsize))
    myrestfreq = str(freq_mean)+'GHz'


    if basename is None:
        basename = os.path.basename(vis)
    myimagename = basename + '.cont.auto'
    rmtables(tablenames=myimagename + '.*')


    if debug:
        print("Mean frequecy:", freq_mean)
        print("Maximum baseline:", baseline_max)
        print("My cell size:", mycell)
        print("My image size:", myimsize)

    if dirty_image:

        print('vis', vis)
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

