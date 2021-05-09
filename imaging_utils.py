# This is file contains the program that related to imaging
# It is designed for ALMACAL project
#
# Author: Jianhang Chen
# Email: cjhastro@gmail.com
# LastUpdate: 2020-11-05
#
# History:
# 2020-11-05: First release

import os
import re
import json
import numpy as np
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt 
from cleanhelper import cleanhelper

# # tbtool
# try:
    # # for casa6
    # from casatool import table as tbtool
    # from casatasks import exportfits
    # from casatasks import rmtables
    # # from casatasks import tclean
# except:
    # from taskinit import tbtool
    # from exportfits_cli import exportfits_cli as exportfits # used by makeSimulatedImage()
    # from rmtables_cli import rmtables_cli as rmtables # used by addGaussianToFITSImage()
    # # from tclean_cli import tclean
try: 
    import analysisUtils as au
    has_analysisUtils = True
except:
    has_analysisUtils = False


# from ms_utils import read_spw

def calculate_sensitivity(vis, full_pwv=False, debug=True):
    """calculate the sensitivity of ALMA data, wrapper of analysisUtils.sensitivity

    """
    # Octile PWV(mm)
    if not isinstance(vis, str):
        raise ValueError("Only single visibility is supported.")
    spw_list = read_spw(vis)
    central_freq = "{:.2f}GHz".format(np.mean(spw_list))
    band_width = "{:.2f}GHz".format(np.sum(np.diff(spw_list)))
    
    antennalist = au.buildConfigurationFile(vis)
    
    #remove the cfg file
    if full_pwv: 
        pwv_list = [0.472, 0.658, 0.913, 1.262, 1.796, 2.748, 5.186]
        sensitivity = []
        for pwv in pwv_list:
            sensitivity.append(au.sensitivity(central_freq, band_width, '60s', pwv=pwv, 
                               antennalist=vis+'.cfg'))
    else:
        pwv = 1.262
        sensitivity = au.sensitivity(central_freq, band_width, '60s', pwv=pwv, antennalist=vis+'.cfg')
    os.system('rm -f {}.cfg'.format(vis))
    return sensitivity

def get_baselines(vis=None,):
    """calculate the baselines
    """
    tb = tbtool()
    tb.open(os.path.join(vis, 'ANTENNA'))
    positions = np.transpose(tb.getcol('POSITION'))
    #stations = tb.getcol('STATION')
    n_ant = len(positions)
    lengths = []
    for i in range(n_ant):
        for j in range(i+1, n_ant):
            x = positions[i][0]-positions[j][0]
            y = positions[i][1]-positions[j][1]
            z = positions[i][2]-positions[j][2]
            print([i,j], x,y,z)
            lengths.append((x**2 + y**2 + z**2)**0.5)
    return lengths

def make_cont_img(vis=None, basename=None, clean=False, myimagename=None, baseline_percent=90,
                  mycell=None, myimsize=None, outdir='./', fov_scale=2.0, imgsize_scale=1, 
                  myuvtaper=None, uvtaper=[], pbcor=True,
                  cellsize_scale=1, datacolumn="corrected", specmode='mfs', outframe="LSRK", 
                  weighting='natural', niter=0, interactive=False, usemask='auto-multithresh', 
                  threshold=None, auto_threshold=5.0, only_fits=False, suffix='', 
                  uvtaper_scale=None, debug=False, **kwargs):

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
    antenna_diameter = np.max(antenna_diameter_list) * u.m

    wavelength = const.c / (freq_mean * u.GHz) # in um
    fov = (fov_scale * 1.22 * wavelength / antenna_diameter * 206265).decompose()

    # calcuate the cell size
    if mycell is None:
        if myuvtaper:
            if isinstance(myuvtaper, list):
                myuvtaper_value = u.Unit(myuvtaper[0]).to(u.arcsec)
            else:
                myuvtaper_value = u.Unit(myuvtaper).to(u.arcsec)
            uvtaper = str(myuvtaper_value)+'arcsec'
            cellsize = np.sqrt((206265 / (baseline_typical / wavelength).decompose())**2 
                               + myuvtaper_value**2)/ 7.0
        else:
            cellsize = 206265 / (baseline_typical / wavelength).decompose()/ 7.0
        mycell = "{:.4f}arcsec".format(cellsize)
    # calculate image size 
    if myimsize is None:
        myimsize = cleanhelper.getOptimumSize(int(imgsize_scale * fov / cellsize))
        print(">>>", myimsize)
    # calculate frequecy
    myrestfreq = str(freq_mean)+'GHz'
    # calcuate threshold
    if not threshold:
        if has_analysisUtils:
            threshold = "{}mJy".format(1000.0*auto_threshold*calculate_sensitivity(vis))
        else:
            print("Warning: no analysisUtils found, set threshold to 0.0!")
            threshold = 0.0

    if basename is None:
        if isinstance(vis, list):
            basename = 'concat'
        else:
            basename = os.path.basename(vis)
    if myimagename is None:
        myimagename = os.path.join(outdir, basename + suffix)
    # else:
        # myimagename = o.path.join(outdir, myimagename)
    if debug:
        print("mean frequecy:", freq_mean)
        print("maximum baseline:", baseline_typical)
        print("cell size:", mycell)
        print("image size:", myimsize)

    if isinstance(vis, list):
        if len(vis) > 4:
            vis_combined = '/tmp/vis_combined.ms'
            concat(vis=vis, concatvis=vis_combined)
            vis = vis_combined
    if clean:
        rmtables('{}.*'.format(myimagename))
        os.system('rm -rf {}.fits'.format(myimagename))
        tclean(vis=vis,
               imagename=myimagename,
               imsize=myimsize, 
               cell=mycell, 
               restfreq=myrestfreq, 
               datacolumn=datacolumn, 
               specmode=specmode, 
               outframe=outframe, 
               weighting=weighting,
               uvtaper=uvtaper,
               niter=niter, 
               threshold=threshold,
               interactive=interactive, 
               usemask=usemask,
               **kwargs)
        if pbcor:
            impbcor(imagename=myimagename+'.image',
                    pbimage=myimagename+'.pb',
                    outfile=myimagename+'.pbcor.image')

        if uvtaper_scale:
            if isinstance(uvtaper_scale, (int, float)):
                uvtaper_scale = [uvtaper_scale,]
            for uvt_param in uvtaper_scale:
                uvt_scale = np.sqrt(1.0*uvt_param**2 - 1.0)
                img_header = imhead(myimagename+'.image')
                #print(img_header)
                restoringbeam = img_header['restoringbeam']
                #print(restoringbeam)
                beam_major = restoringbeam['major']
                beam_minor = restoringbeam['minor']
                beam_pa = restoringbeam['positionangle']
                # set-up the uvtaper
                bmaj = "{:.6f}{}".format(beam_major['value']*uvt_scale, beam_major['unit'])
                bmin = "{:.6f}{}".format(beam_minor['value']*uvt_scale, beam_minor['unit'])
                bpa = "{:.6f}{}".format(beam_pa['value'], beam_pa['unit'])
                # set up the resolution and imsize
                mycell = "{:.6f}arcsec".format(cellsize * uvt_param)
                myimsize = cleanhelper.getOptimumSize(int(imgsize_scale * fov / (cellsize*uvt_param)))
                # myimsize = efficient_imsize(int(myimsize / np.sqrt(uvt_scale**2+1.)))
                if debug:
                    print("uvtaper parameters: \n")
                    print("uvtaper:", [bmaj, bmin, bpa])
                    print("uvtaper cell size:", mycell)
                    print("uvtaper image size:", myimsize)

                os.system('rm -rf {}.uvscale{}.fits'.format(myimagename, uvt_param))
                myimagename_uvtaper = myimagename+'.uvscale{}'.format(uvt_param)
                tclean(vis=vis,
                      imagename=myimagename_uvtaper,
                      imsize=myimsize, 
                      cell=mycell, 
                      restfreq=myrestfreq, 
                      datacolumn=datacolumn, 
                      specmode=specmode, 
                      outframe=outframe, 
                      weighting=weighting,
                      niter=niter, 
                      threshold=threshold,
                      interactive=interactive, 
                      usemask=usemask,
                      uvtaper=[bmaj, bmin, bpa],
                      **kwargs)
                if pbcor:
                    impbcor(imagename=myimagename_uvtaper+'.image',
                            pbimage=myimagename_uvtaper+'.pb',
                            outfile=myimagename_uvtaper+'.pbcor.image')

    if isinstance(vis, list):
        os.system('rm -rf /tmp/vis_combined.ms')
    
    if only_fits:
        for image in glob.glob(myimagename+'*.image'):
            exportfits(imagename=image, fitsimage=image+'.fits')
        rmtables(tablenames=myimagename+'.*')

def image_selfcal(vis=None, ncycle=3, ):
    # tclean and self-calibration

    pass

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


