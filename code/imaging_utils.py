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

from readms import read_spw

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
        

def make_cont_img(vis=None, basename=None, dirty_image=False, clean_image=False, myimagename=None, mycell=None, myimsize=None, outdir='./', dry_run=False, 
                  baseline_percent=90, **kwargs):
    """This function is used to make the continuum image
    
    dirty_image: generate the dirty image

    clean_image: make the clean image with tclean
        You need to specify: weighting, niter, interactive
        The imagename, imsize, cell and restfreq will be set automatically if no value is provided by kwargs
    """
    spw_specrange = read_spw(vis)
    freq_mean = np.mean(spw_specrange) # in GHz
    
    # get the baselines
    tb = tbtool()
    if isinstance(vis, list):
        baselines_list = []
        for v in vis:
            obs_baselines = au.getBaselineLengths(v, returnLengthsOnly=True)
            baselines_list.append(obs_baselines)
        # unfold the list
        baselines_list = [item for sublist in baselines_list for item in sublist]
    if isinstance(vis, str):
        baselines_list = au.getBaselineLengths(vis, returnLengthsOnly=True)
        
        # # another way is to use getBaselineStats directly
        # baseline_stat = au.getBaselineStats(vis)
        # baseline_max = baseline_stat[2] * u.m
        # baseline_90 = baseline_stat[-1] * u.m
        # baseline_typical = baseline_90
    baseline_typical = np.percentile(baselines_list, baseline_percent) * u.m

    # read the antenna_diameter
    if isinstance(vis, list):
        antenna_diameter_list = []
        for v in vis:
            tb.open(v + '/ANTENNA')
            antenna_diameter_list.append(tb.getcol('DISH_DIAMETER'))
            tb.close()
        antenna_diameter_list = [item for sublist in antenna_diameter_list for item in sublist]
    if isinstance(vis, str):
        tb.open(vis + '/ANTENNA')
        antenna_diameter_list = tb.getcol('DISH_DIAMETER')
        tb.close()
    antenna_diameter = np.median(antenna_diameter_list) * u.m

    wavelength = const.c / (freq_mean * u.GHz) # in um
    fov = (2.0 * 1.22 * wavelength / antenna_diameter * 206265).decompose()

    # calcuate the cell size
    if mycell is None:
        cellsize = 206265 / (baseline_typical / wavelength).decompose() / 6
        mycell = "{:.4f}arcsec".format(cellsize)
    
    if myimsize is None:
        myimsize = efficient_imsize(int(fov / cellsize))
        myrestfreq = str(freq_mean)+'GHz'


    if basename is None:
        if isinstance(vis, list):
            basename = 'concat'
        else:
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

def image_selfcal(vis=None, ncycle=3, ):
    # tclean and self-calibration
    os.system('rm -rf tclean/cycle0.*')
    tclean(vis=calmsfile,
           imagename='./tclean/cycle0',
           field='',
           spw='',
           datacolumn='data',
           specmode='mfs',
           deconvolver='hogbom',
           nterms=1,
           gridder='standard',
           imsize=320,
           cell='0.1arcsec',
           pblimit=0.1,
           weighting='natural',
           threshold='0mJy',
           niter=50000,
           interactive=True,
           savemodel='modelcolumn')


    os.system("rm -rf cycle1_phase.cal")
    gaincal(vis=calmsfile,
            caltable="cycle1_phase.cal",
            field="",
            solint="2min",
            calmode="p",
            refant=myrefant,
            gaintype="G")

    applycal(vis=calmsfile,
             field="",
             gaintable=["cycle1_phase.cal"],
             interp="linear")

    os.system("rm -rf {0} {0}.flagversions".format(calmsfile+'.selfcal'))
    split(vis=calmsfile,
          outputvis=calmsfile+'.selfcal',
          datacolumn="corrected")

    # make the cycle 1 image
    os.system('rm -rf tclean/cycle1.*')
    tclean(vis=calmsfile+'.selfcal',
           imagename='./tclean/cycle1',
           field='',
           spw='',
           specmode='mfs',
           deconvolver='hogbom',
           datacolumn='data',
           nterms=1,
           gridder='standard',
           imsize=320,
           cell='0.1arcsec',
           pblimit=0.1,
           weighting='natural',
           threshold='0mJy',
           niter=50000,
           interactive=True,
           savemodel='modelcolumn')


