# The main class file for ALMACAL project

import os
import sys
import re
import numpy as np

from astropy.coordinates import SkyCoord

# import casatools and casatasks
# try:
    # import casalith
    # casaVersion = casalith.version_string()
# except:
    # import casadef
    # casaVersion = casadef.casa_version

# if casaVersion > '5.9.9':
    # from casatools import table as tbtool
    # from casatasks import exportfits
    # from casatasks import importfits
    # from casatasks import rmtables
# else:
    # from taskinit import tbtool
    # from exportfits_cli import exportfits_cli as exportfits # used by makeSimulatedImage()
    # from importfits_cli import importfits_cli as importfits # used by addGaussianToFITSImage()
    # from rmtables_cli import rmtables_cli as rmtables # used by addGaussianToFITSImage()

# from imaging_utils import make_cont_img
# from ms_utils import read_spw


class ALMACAL(object):
    """The Project class
    """
    def __init__(self, basedir=None):
        self.basedir = os.path.abspath(basedir)
        self.objlist = self.get_objlist()

    def get_objlist(self):
        p_obj = re.compile('J\d+[+-]\d')
        objlist = []
        for item in os.listdir(self.basedir):
            if p_obj.match(item):
                objlist.append(os.path.join(self.basedir, item))
        return objlist
        

class Calibrator(object):
    """The Calibrator class
    """

    def __init__(self, basedir=None, name=None):
        p_obj = re.compile('J\d+[+-]\d')
        self.basedir = os.path.abspath(basedir)
        if not name:
            if p_obj.match(os.path.basename(basedir)):
                self.name = os.path.basename(basedir)
            else:
                self.name = 'Undefined'
        self.obslist = self.get_obslist()
        
    def get_obslist(self, obslist=None, fullpath=True):
        p_obs = re.compile('uid___')
        if obslist:
            if isinstance(obslist, str):
                obslist_fromfile = []
                with open(obslist) as f:
                    obslist_lines = f.readlines()
                for item in obslist_lines:
                    line = item.strip()
                    if p_obs.match(line):
                        obslist_fromfile.append(line)
                obslist = obslist_fromfile
        else:
            obslist = []
            obslist_fromdir = os.listdir(self.basedir)
            for item in obslist_fromdir:
                if p_obs.match(item):
                    obslist.append(item)

        if fullpath:
            obslist_fullpath = []
            for item in obslist:
                obslist_fullpath.append(os.path.join(self.basedir, item))
            return obslist_fullpath
        return obslist

    def get_obs(self, idx):
        return SingleObs(self.obslist[idx])

class SingleObs(object):
    """The single observation
    """
    def __init__(self, vis=None):
        self.vis = os.path.abspath(vis)
        self.dirname = os.path.dirname(self.vis)
        self.name = os.path.basename(self.vis)
        p_obs = re.compile('(?P<obsname>uid___\w*\.ms(\.split\.cal)?\.(?P<objname>[\s\w+-]+)_(?P<band>B\d+))')
        if p_obs.match(self.name):
            obs_matched = p_obs.search(self.name).groupdict()
            self.objname = obs_matched['objname']
            self.band = obs_matched['band']

    def copy(self, outdir='./', outfilename=None, overwrite=False):
        """copy the original visibility
        """
        if outfilename is None:
            outfilename = self.name
        outvis = os.path.join(outdir, outfilename)
        if overwrite:
            os.system('cp -rf {} {}'.format(self.vis, outvis))
        else:
            os.system('cp -rn {} {}'.format(self.vis, outvis))
        return outvis
    


    def gen_image(self, outdir='./', imagename=None, niter=0, exclude_aca=False, 
                  debug=False, **kwargs):
        """make images for one calibrator on all or specific band

        Params:
        objfolder: the folder that contains all the observation of one calibrator
        band: the band want to imaged
        outdir: the directory where all the fits image will be placed
        **kwargs: the additional parameters supported by make_cont_img

        """
        if not imagename:
            imagename = os.path.join(outdir, self.name + '.cont.auto')

        if exclude_aca:
            try:
                tb.open(self.vis + '/ANTENNA')
                antenna_diameter = np.mean(tb.getcol('DISH_DIAMETER'))
                tb.close()
            except:
                return False
            if antenna_diameter < 12.0:
                if debug:
                    print("Excuding data from {}".format(antenna_diameter))
                return False
        try:
            make_cont_img(vis=self.vis, clean=True, myimagename=imagename, outdir=outdir, 
                          niter=niter, **kwargs)
        except:
            print("Error in imaging {}".format(self.vis))
        exportfits(imagename=imagename+'.image', fitsimage=imagename+'.fits')
        rmtables(tablenames=imagename+'.*')

    def remove_calibrator(self, copy=True, outdir='./'):
        if copy:
            vis = self.copy(outdir=outdir)
        else:
            vis = self.vis
        outfile = os.path.join(outdir, self.name+".point.cl")
        uvmodelfit(vis=vis, niter=5, comptype="P", sourcepar=[1.0, 0.0, 0.0], 
                   varypar=[True,False,False], outfile=outfile)
        # uvmodelfit(vis=vis, niter=5, comptype="G", outfile=outfile,
                   # sourcepar=[1.0, 0.0, 0.0, 1.0, 0.9, 0], 
                   # varypar=[True, False, False, True, True, True])
        ft(vis=vis, complist=outfile)
        uvsub(vis=vis,reverse=False)
        # outputvis = os.path.join(outdir, self.name)
        # os.system('mkdir -p {}'.format(os.path.dirname(outputvis)))
        split(vis=vis, outputvis=outputvis)
       
