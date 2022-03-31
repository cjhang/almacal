# This file include all the tasks for ALMACAL number counts projects

# The task list include:
# 1. gen_obstime: generate the observational time for the whole dataset
# 2. gen_dirty_image: generate the dirty image for all the observations of one calibrator
# 3. show_image: interactive way to inspect the images manually 
import os
import glob
import re
import pickle
import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib import patches
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit
from scipy import interpolate
from astropy.stats import sigma_clipped_stats, sigma_clip

import multiprocessing as mp
import analysisUtils as au
# from analysisUtils import tbtool
# from imaging_utils import make_cont_img

tb = au.tbtool()

def flatten(t):
    flat_list = [item for sublist in t for item in sublist]
    return flat_list

def gen_filenames(listfile=None, dirname=None, list_array=None, basedir=None, debug=False, 
        exclude_aca=True, suffix='', remove_prefix=True):
    """generate all the valid files names

    """
    file_list = []
    if dirname:
        obj_match = re.compile('^J\d*[+-]\d*$')
        obs_match = re.compile('(?P<obsname>uid___\w*\.ms(\.split\.cal)?\.(?P<objname>[\s\w+-]+)_(?P<band>B\d+))')
        if not os.path.isdir(dirname):
            raise ValueError('Invalid directory!')
        for item in os.listdir(dirname):
            if debug:
                print(item)
            if obj_match.match(item):
                obj_path = os.path.join(dirname, item)
                for obs in os.listdir(obj_path):
                    if obs_match.match(obs):
                        if exclude_aca:
                            try:
                                tb.open(vis + '/ANTENNA')
                                antenna_diameter = np.mean(tb.getcol('DISH_DIAMETER'))
                                tb.close()
                            except:
                                continue
                            if antenna_diameter < 12.0:
                                if debug:
                                    print("Excuding data from {}".format(antenna_diameter))
                                continue
                        file_list.append(os.path.join(obj_path, obs))
            elif obs_match.match(item):
                file_list.append(os.path.join(dirname, item))
    elif listfile:
        with open(listfile) as f:
            listfile_lines = f.readlines()
        for l in listfile_lines:
            if '/' in l:
                line = l.split('/')[-1].strip()
            else:
                line = l.strip()
            if basedir:
                line = os.path.join(basedir, line)
            file_list.append(line+suffix)
    elif list_array:
        for item in list_array:
            if basedir:
                item_fullpath = os.path.join(basedir, item)
            file_list.append(item_fullpath+suffix)
    return file_list

def savelist(l, filename=None, outdir='./', file_mode='a+'):
    with open(os.path.join(outdir, filename), file_mode) as f:
        for item in l:
            f.write(item+'\n')

def search_obs(basedir, config='main', band='B6', debug=False, config_select=True,
        time_select=True, start_time='2010-01-01T00:00:00', end_time='2050-01-01T00:00:00',
        freq_select=False, star_freq=100, end_freq=101):
    band_match = re.compile('_(?P<band>B\d{1,2})$')
    obs_match = re.compile('^uid___')
    if config == 'main':
        dish_diameter = 12.0
    elif config == 'aca':
        dish_diameter = 7.0
    filelist = []
    for obs in os.listdir(basedir):
        if debug:
            print(obs)
        obs_fullpath = os.path.join(basedir, obs)
        if obs_match.match(obs):
            if band_match.search(obs):
                obs_band = band_match.search(obs).groupdict()['band']
                if obs_band != band:
                    continue
            if time_select:
                start_time = Time(start_time)
                end_time = Time(end_time)
                try:
                    tb.open(obs_fullpath)
                    obs_time = Time(tb.getcol('TIME').max()/24/3600, format='mjd')
                    tb.close()
                except:
                    if debug:
                        print("Error in opening the visibility!")
                    continue
                if debug:
                    print('> obs_time', obs_time.iso)
                    print('> start_time', start_time.iso)
                    print('> end_time', end_time.iso)
                if obs_time < start_time or obs_time > end_time:
                    if debug:
                        print(">> Skip by wrong observation time: {}".format(obs))
                    continue

            if config_select:
                try:
                    tb.open(obs_fullpath + '/ANTENNA')
                    antenna_diameter = np.mean(tb.getcol('DISH_DIAMETER'))
                    tb.close()
                except:
                    continue

                if np.abs(antenna_diameter - dish_diameter) > 1e-4:
                    continue

            if freq_select:
                is_freq_covered = False
                spw_list = spw_stat(obs_fullpath)
                for freq in spw_list[band][0]:
                    if freq[0] <= star_freq and freq[-1] >= end_freq:
                        is_freq_covered = True
                        break
                if not is_freq_covered:
                    continue
            filelist.append(obs)
    return filelist

def gen_image(vis=None, band=None, outdir='./', niter=0, suffix='.cont.auto',
              check_validity=False, overwrite=False,
              check_aca=True, check_sensitivity=True, check_gaussian_noise=True,
              check_flux=True, minimal_fluxval=0.001,
              debug=False, **kwargs):
    """make images for visibility and do the first step checking for the valid

    Params:
    objfolder: the folder that contains all the observation of one calibrator
    band: the band want to imaged
    outdir: the directory where all the fits image will be placed
    **kwargs: the additional parameters supported by make_cont_img

    """
    p_obs = re.compile('uid___')
    
    if not p_obs.search(vis):
        print("Invalid name!")
        return False
    basename = os.path.basename(vis)
    myimagename = os.path.join(outdir, basename + suffix)
    fitsimage = myimagename+'.image.fits'

    is_valid = True
    if check_validity:
        if check_aca:
            try:
                tb.open(vis + '/ANTENNA')
                antenna_diameter = np.mean(tb.getcol('DISH_DIAMETER'))
                tb.close()
            except:
                is_valid = False
            if antenna_diameter < 12.0:
                if debug:
                    print("Excuding data from {}".format(antenna_diameter))
                is_valid = False

        if check_flux:
            # check the flux value of the subtracted point source
            p_uidimg = re.compile('(?P<uidname>uid___\w+.ms.split.cal.J\d+[+-]\d+_B\d+).cont.auto.fits')
            if p_uidimg.search(fitsimage):
                uidname = p_uidimg.search(fitsimage).groupdict()['uidname']
                fluxval = search_flux(uidname)
                if debug:
                    print(uidname)
                    print("fluxval: {}".format(fluxval))
                if fluxval > 0:
                    if fluxval/minimal_fluxval < 1:
                        is_valid = False
        # if check_sensitivity:
            # is_valid = check_vis2image(vis, imagefile=fitsimage)
        # if check_resolution:
            # with fits.open(fitsimage) as hdu:
                # header = hdu[0].header
    # if not is_valid:
        # print("Not passing validity check!")
    if os.path.isfile(fitsimage) and not overwrite:
        if not is_valid:
            print("Removing invalid image: {}".format(fitsimage))
            os.system('rm {}'.format(fitsimage))
        else:
            print('Using existing image {}'.format(fitsimage))
            return fitsimage
    elif is_valid:
        # check the intent
        intent_list = check_intent(vis)
        intent_list_full = []
        for i in intent_list:
            intent_list_full.append('CALIBRATE_{}#ON_SOURCE'.format(i))
        intent = ','.join(intent_list_full)
        if debug:
            print('intent', intent)
        try:
            make_cont_img(vis=vis, clean=True, myimagename=myimagename, outdir=outdir, niter=niter, intent=intent, **kwargs)
            exportfits(imagename=myimagename+'.image', fitsimage=fitsimage)
            rmtables(tablenames=myimagename+'.*')

        except:
            print("Error in imaging {}".format(vis))
 
    return fitsimage

def gen_images(vis=None, dirname=None, outdir='./', bands=None, exclude_aca=True, 
                  debug=False, **kwargs):
    """generate the images of all calibrators

    Args:
        allcal_dir: the root directory contaCODEins all the measurements of calibrators
        vis: the single visibility or list
        outdir: the output directory
        bands: the bands to be imaged
        **kwargs: the additional parameters of make_cont_img

    Examples:
        1. to get all the images for the whole project
            
            gen_all_image('/science_ALMACAL/data', '/tmp/all_images')

        2. get all the images just for one folder
    """
    os.system('mkdir -p {}'.format(outdir))
    if vis:
        if isinstance(vis, str):
            filelist = [vis,]
        elif isinstance(vis, list):
            filelist = vis
    if dirname:
        filelist = []
        obs_match = re.compile('^uid___')
        for obs in os.listdir(dirname):
            if debug:
                print(obs)
            if obs_match.match(obs):
                filelist.append(os.path.join(dirname, obs))
    if debug:
        print('filelist:', filelist)
    if bands is not None:
        band_match = re.compile('_(?P<band>B\d{1,2})$')
        for infile in filelist:
            if band_match.search(vis):
                obs_band = band_match.search(vis).groupdict()['band']
                if obs_band not in bands:
                    if debug:
                        print("Skip {}, mis-match band".format(infile))
                    continue
                    os.system('mkdir -p {}'.format(outdir_band))
                    gen_image(infile, band=band, outdir=outdir_band, **kwargs)
    else:
        for infile in filelist:
            gen_image(infile, outdir=outdir, **kwargs)

def show_images(fileglob=None, filelist=None, basedir=None, mode='auto', nrow=3, ncol=3, savefile=None, debug=False):
    """show images in an interactive ways, and record the input from inspector
    
    Parameters:
        fileglob: the match patter of the image filenames, like: './*.image'
    """
    
    if fileglob:
        if debug:
            print(fileglob)
        all_files = glob.glob(fileglob)
    if filelist:
        if debug:
            print(filelist)
        if isinstance(filelist, str):
            if not os.path.isfile(filelist):
                raise ValueError("{} is not found".format(filelist))
            try:
                flist = []
                with open(filelist) as f:
                    filelist_lines = f.readlines()
                for line in filelist_lines:
                    flist.append(line.strip())
            except:
                raise ValueError("file {} cannot be open".format(filelist))
        elif isinstance(filelist, list):
            flist = filelist
        else:
            raise ValueError("Wrong file type of filelist!")
        if basedir:
            all_files = []
            for item in flist:
                all_files.append(os.path.join(basedir, item))
        else:
            all_files = flist

    total_num = len(all_files)
    if total_num == 1:
        ncol = nrow = 1
    select_num = 0
    print("Find {} files".format(total_num))
    all_select = []
    if mode == 'random':
        all_files = np.random.choice(all_files, nrow * ncol)
    for i in range(0, len(all_files), ncol*nrow):
        fig = plt.figure(figsize=(4*nrow,4*ncol))
        for j in range(0, nrow*ncol):
            if (i+j)>=len(all_files):
                continue
            try:
                with fits.open(all_files[i+j]) as fitsimage:
                    imageheader = fitsimage[0].header
                    imagedata = fitsimage[0].data
                    wcs = WCS(imageheader)
                    # wcs2 = wcs.dropaxis(3)
                    # wcs2 = wcs2.dropaxis(2)

                    #with wcs projection
                    scale = np.abs(imageheader['CDELT1'])*3600
                    ny, nx = imagedata.shape[-2:]
                    x_index = (np.arange(0, nx) - nx/2.0) * scale
                    y_index = (np.arange(0, ny) - ny/2.0) * scale
                    x_map, y_map = np.meshgrid(x_index, y_index)

            except:
                print("Error in reading: {}".format(all_files[i+j]))
                continue
            ax = fig.add_subplot(nrow, ncol, j+1)#, projection=wcs2, slices=(0, 0, 'x', 'y'))
            # ax.text(10, 10, str(j), fontsize=20)
            ax.set_title(str(j+1))
            #ax.imshow(imagedata[0,0,:,:], origin='lower')#, cmap='viridis')
            imagedata = np.ma.masked_invalid(imagedata)
            #imagedata = imagedata.filled(0)
            ax.pcolormesh(x_map, y_map, imagedata[0,0,:,:])
            ax.text(0, 0, '+', color='r', fontsize=24, fontweight=100, horizontalalignment='center',
                    verticalalignment='center')
            ax.set_xlabel('RA [arcsec]')
            ax.set_ylabel('Dec [arcsec]')
            ax.set_title(os.path.basename(all_files[i+j][:15]))
        # show the image and record the 
        plt.show()
        print('Input the index of images (1-9), seperate with comma:')
        try:
            find_zero = False
            idx_input = input()
            if idx_input == 0:
                print("Currently at {}/{}".format(i, len(all_files)))
                plt.close('all')
                break
            if isinstance(idx_input, int):
                idx_input = [idx_input]
            select_num += len(idx_input)
            for ind in idx_input:
                if ind == 0:
                    print("Currently at {}/{}".format(i, len(all_files)))
                    plt.close('all')
                    find_zero = True
                    break
                all_select.append(all_files[i+ind-1])
            if find_zero:
                break
        except:
            plt.close('all')
            continue
        # plt.clf()
        plt.close('all')
    
    print("Totally {:.2}% of data have been selected.".format(100.*select_num/(total_num+1e-6)))
    #print(all_select)
    if savefile:
        obsname_match = re.compile('(?P<obsname>uid___\w*\.ms\.split\.cal\.J\d*[+-]+\d*_B\d+)')
        with open(savefile, 'w+') as f:
            for item in all_select:
                try:
                    obsname = obsname_match.search(item).groupdict()['obsname']
                    f.write(obsname+'\n')
                except:
                    print("Error in matching the obs name for filname: {}".format(item))
                    continue
    else:
        for item in all_select:
            print(item)

def make_combine_obs(obj, basedir=None, band=None, badfiles=None):
    badfiles_list = []
    if badfiles is not None:
        with open(badfiles) as f:
            badfiles_readlines = f.readlines()
        for item in badfiles_readlines:
            # remove the '\n' at the end of item
            badfiles_list.append(item.strip())

    all_files = os.listdir(obj)
    valid_files = []
    for obs in all_files:
        if obs in badfiles_list:
            continue
        band_match = re.compile('_(?P<band>B\d{1,2})$')
        try:
            obs_band = band_match.search(obs).groupdict()['band']
        except:
            continue
        if obs_band == band:
            if basedir is not None:
                valid_files.append(os.path.join(basedir, obs))
            else:
                valid_files.append(obs)

    # print(valid_files)
    return valid_files
    
def search_flux(obs, allflux_file=None, strick_mode=False, debug=False):
    obs_filename = os.path.basename(obs)
    if allflux_file is None:
        if 'ALMACAL_NUMBERCOUNTS_HOME' in os.environ.keys():
            root_path = os.environ['ALMACAL_NUMBERCOUNTS_HOME']
        else:
            root_path = os.path.join(os.path.expanduser('~'), 'projects/almacal/number_counts')
        allflux_file = os.path.join(root_path, 'code/data/allcal.fluxval')
    if debug:
        print(allflux_file)
    # allflux_file = './allcal.fluxval'
    try:
        allflux = np.loadtxt(allflux_file, dtype=str)
    except:
        if strick_mode:
            raise ValueError("Unsupported allflux file!")
        return 0
    obs_select = allflux[:,1] == obs_filename

    if np.sum(obs_select) < 1:
        if strick_mode:
            raise ValueError("No flux can be found!")
        return 0
    if np.sum(obs_select) > 1:
        print("Warning: not an unique flux!")

    return np.float(allflux[:,3][obs_select][0])

def ms_restore(obs_list, allflux_file=None, basedir=None, outdir='./output', tmpdir='./tmp', debug=True):
    """put back the central point source
    """

    # if not os.path.isdir(outdir):
    os.system('mkdir -p {}'.format(outdir))
    # if not os.path.isdir(tmpdir):
    os.system('mkdir -p {}'.format(tmpdir))

    obs_match = re.compile('(?P<obsname>uid___\w*\.ms(\.split\.cal)?\.(?P<objname>[\s\w+-]+)_(?P<band>B\d+))')
    
    if isinstance(obs_list, str):
        if os.path.isfile(obs_list):
            file_list = []
            with open(obs_list) as f:
                file_list_readlines = f.readlines()
            for item in file_list_readlines:
                # remove the '\n' at the end of item
                file_list.append(item.strip())
            obs_list = file_list
        else:
            obs_list = [obs_list,]
    if basedir is not None:
        lamf = lambda x: os.path.join(basedir, x)
        obs_list = map(lamf, obs_list)
    for obs in obs_list:
        if debug:
            print(">>>> {}".format(obs))
        try:
            obs_matched = obs_match.search(obs).groupdict()
        except:
            print("Unsupported filename:".format(obs))

        obsname = obs_matched['obsname']
        objname = obs_matched['objname']
        band = obs_matched['band']

        # read the basic information from ms
        mydirection = read_refdir(obs).encode('utf-8')
        spw_list = read_spw(obs)
        myfreq = str(np.mean(spw_list)) + 'GHz'
        # read the flux from the fluxval file
        try:
            myflux = search_flux(os.path.basename(obs), allflux_file=allflux_file)
        except:
            print("No fluxval found for {}".format(obs))
            continue

        # create a point source
        comp_cal = os.path.join(tmpdir, obsname+'_central_cal.cl')
        os.system('rm -rf {}'.format(comp_cal))
        cl.done()
        cl.addcomponent(dir=mydirection, flux=myflux, fluxunit='Jy', freq=myfreq, shape="point", spectrumtype='constant')
        cl.rename(comp_cal)
        cl.done()


        basedir = os.path.dirname(obs)
        filename = os.path.basename(obs)
        if outdir is None:
            outdir = basedir
        # setup the temperary files
        tmpfile = os.path.join(tmpdir, filename+'.tmp')
        os.system('rm -rf {}'.format(tmpfile))
        os.system('cp -r {} {}'.format(obs, tmpfile))
        # restore the point source to the uv data
        ft(vis=tmpfile, complist=comp_cal)
        uvsub(vis=tmpfile, reverse=True)
        # split the restore data
        outfile = os.path.join(outdir, filename+'.restore')
        os.system('rm -rf {}'.format(outfile))
        split(vis=tmpfile, outputvis=outfile, datacolumn='corrected')
    # remove all the temperary files
    os.system('rm -rf {}'.format(tmpdir))

def copy_ms(basedir=None, outdir=None, selectfile=None, debug=False, time_select=False, 
            start_time='2010-01-01T00:00:00', end_time='2050-01-01T00:00:00', 
            select_band=None, overwrite=False):
    p_obs = re.compile('uid___')
    band_match = re.compile('_(?P<band>B\d{1,2})$')
    selectfiles_list = []
    start_time = Time(start_time)
    end_time = Time(end_time)
    if isinstance(selectfile, str):
        if not os.path.isfile(selectfile):
            raise ValueError("Cannot open the file {}".format(selectfile))
        with open(selectfile) as f:
            selectfiles_readlines = f.readlines()
        for item in selectfiles_readlines:
            # remove the '\n' at the end of item
            selectfiles_list.append(item.strip())
    elif isinstance(selectfile, list):
        selectfiles_list = selectfile

    else:
        selectfiles_list = []
        all_files = os.listdir(basedir)
        for obs in all_files:
            if not p_obs.match(obs):
                if debug:
                    print("Not a valid file: {}".format(obs))
                continue
            if debug:
                print('Find {}'.format(obs))
            if time_select:
                try:
                    tb.open(os.path.join(basedir, obs))
                    obs_time = Time(tb.getcol('TIME').max()/24/3600, format='mjd')
                    tb.close()
                except:
                    if debug:
                        print("Error in opening the visibility!")
                    continue
                if debug:
                    print('> obs_time', obs_time.iso)
                    print('> start_time', start_time.iso)
                    print('> end_time', end_time.iso)
                if obs_time < start_time or obs_time > end_time:
                    if debug:
                        print(">> Skip by wrong observation time: {}".format(obs))
                    continue
            if select_band:
                if band_match.search(obs):
                    obs_band = band_match.search(obs).groupdict()['band']
                else:
                    continue
                if obs_band not in select_band:
                    if debug:
                        print(">> Skip by wrong band: {}".format(obs))
                    continue
            selectfiles_list.append(obs)
        if debug:
            print(selectfiles_list)
        for obs in selectfiles_list:
            print("Copying {}".format(obs))
            if not os.path.isdir(outdir):
                os.system('mkdir -p {}'.format(outdir))
            if os.path.isdir(os.path.join(outdir, obs)) and not overwrite:
                print("File existing: {}".format(obs))
            else:
                os.system('cp -r {} {}'.format(os.path.join(basedir, obs), outdir))

def gaussian(x, u0, amp, std):
    return amp*np.exp(-0.5*((x-u0)/std)**2)

def check_image(img, plot=False, radius=8, debug=False, sigmaclip=True, 
                minimal_fluxval=0.001, outlier_frac=0.05, gaussian_deviation=0.02, 
                central_median_deviation=1.0, central_mean_deviation=1.0, savefig=False, 
                figname=None, outdir=None, strick_mode=False):
    """This program designed to determine the validity of the image after point source subtraction
    
    The recommended img is the fits image, if not, it will be converted into fits using casa exportfits
    
    Args:
        radius: in units of half major axis of sythesized beam
    """
    p_uidimg = re.compile('(?P<uidname>uid___\w+.ms.split.cal.J\d+[+-]\d+_B\d+).cont.auto.fits')
    if p_uidimg.search(img):
        uidname = p_uidimg.search(img).groupdict()['uidname']
    else:
        uidname = None

    with fits.open(img) as hdu:
        header = hdu[0].header
        data = hdu[0].data
    if outdir is None:
        outdir = os.path.dirname(img)
    ny, nx = data.shape[-2:]
    masked_data = np.ma.masked_invalid(data.reshape(ny, nx))
    mask_invalid = masked_data.mask
    data4hist = masked_data[~mask_invalid]
    # Testing the sigma clip
    if sigmaclip:
        # run only once sigma-clip, to exclude the most extreme value
        clipped_data, clip_low, clip_up = scipy.stats.sigmaclip(data4hist, 10, 10)
        data4hist = clipped_data 
        # clipped_data = sigma_clip(masked_data, sigma=5, iters=5)
        # data4hist = clipped_data.compressed()
    # statistics the noisy of the image
    # hist, bins = np.histogram(masked_data.data[~masked_data.mask], bins=100)
    hist, bins = np.histogram(data4hist, bins=100)
    bins_mid = (bins[:-1] + bins[1:])*0.5
    
    # sigmaclip based statistics
    # clip any bad value to get the most accurate statistics
    mean0, median0, std0 = sigma_clipped_stats(masked_data, sigma=5.0,
            iters=5)
    p0 = (mean0, 1.0, std0)
    #p0 = (0, 1, 1e-4) # mean, max and std
    amp_scale = 1.0*np.max(hist) # change to int into float
    try:
        popt, pcov = curve_fit(gaussian, bins_mid, hist/amp_scale, p0=p0)
    except:
        print("`Fitting failed!")
        popt = p0
        
    hist_fit = gaussian(bins_mid, *popt)*amp_scale
    mean, amp, sigma = popt
    lower_1sigma = mean - 1.0*sigma
    lower_2sigma = mean - 2.0*sigma
    lower_3sigma = mean - 3.0*sigma
    upper_1sigma = mean + 1.0*sigma
    upper_2sigma = mean + 2.0*sigma
    upper_3sigma = mean + 3.0*sigma
    ## checking the fitting
    # the fraction of the 3 sigma outlier, theoretical gaussian value is 
    n_3sigma = np.sum(hist[(bins_mid > lower_3sigma) & (bins_mid < upper_3sigma)])
    percent_3sigma = 1.0 * n_3sigma / np.sum(hist)
    percent_2sigma = 1.0 * np.sum(hist[(bins_mid > lower_2sigma) & (bins_mid < upper_2sigma)]) / np.sum(hist) 
    percent_1sigma = 1.0 * np.sum(hist[(bins_mid > lower_1sigma) & (bins_mid < upper_1sigma)]) / np.sum(hist) 
    # calculating the deviation from Gaussian
    deviation_1sigma = np.abs((percent_1sigma - 0.6827)/0.6827)
    deviation_2sigma = np.abs((percent_2sigma - 0.9545)/0.9545)
    deviation_3sigma = np.abs((percent_3sigma - 0.9973)/0.9973)

    # statistics in the central region
    bmaj = header['BMAJ']
    bmin = header['BMIN']
    bpa = header['BPA']
    # print(bmaj, bmin, bpa)
    # radius = 5 
    bmaj_pixel_size = bmaj / np.abs(header['CDELT1'])
    bmin_pixel_size = bmaj / np.abs(header['CDELT2'])
    # select the central region with side length of radius*bmaj 
    x_index = np.arange(0, nx) - nx/2.0
    y_index = np.arange(0, ny) - ny/2.0
    x_map, y_map = np.meshgrid(x_index, y_index)
    mask = np.sqrt(x_map**2 + y_map**2) < radius*bmaj_pixel_size*0.5 #units in half major axis 
    mask = ~masked_data.mask & mask # also masked the central inf and nan

    hist_center, bins_center = np.histogram(masked_data[mask], bins=bins)
    bins_mid_center = (bins_center[:-1] + bins_center[1:])*0.5
    amp_scale_center = 1.0*np.max(hist_center) # change to int into float

    mean_central = np.ma.mean(masked_data[mask])
    median_central = np.ma.median(masked_data[mask])
    
    n_outlier_central = np.sum(hist_center[bins_mid_center>upper_3sigma])
    n_outlier_central_neg = np.sum(hist_center[bins_mid_center<lower_3sigma])
    percent_outlier_central = 1.0*n_outlier_central/np.sum(hist_center+1e-6)
    percent_outlier_central_neg = 1.0*n_outlier_central_neg/np.sum(hist_center+1e-6)
    
    if debug:
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print('>> Checking the fitting of noise:')
        print('mean: {}, sigma:{}'.format(mean, sigma))
        # print('number of 3 sigma outlier: {}'.format(n_outlier_3sigma))
        print('fraction of 1 sigma outlier: {:.4f}%. [theoretical: 68.27%]'.format(percent_1sigma*100.))
        print('deviation of 1 sigma: {:.4f}%.'.format(deviation_1sigma*100.))
        print('fraction of 2 sigma outlier: {:.4f}%. [theoretical: 95.45%]'.format(percent_2sigma*100.))
        print('deviation of 2 sigma: {:.4f}%.'.format(deviation_2sigma*100.))
        print('fraction of 3 sigma outlier: {:.4f}%. [theoretical: 99.73%]'.format(percent_3sigma*100.))
        print('deviation of 3 sigma: {:.4f}%.'.format(deviation_3sigma*100.))
        print('>> Statistics of central region')
        print("central mean: {}".format(mean_central))
        print("central median: {}".format(median_central))
        print('deviation of central median: {:.4f}.'.format(np.abs(median_central)/sigma))
        print('number of 3sigma outlier of central region: {}'.format(n_outlier_central))
        print('fraction of 3sigma outlier of central region: {}'.format(percent_outlier_central))

    if plot:
        # show the image
        fig = plt.figure(figsize=(16, 4.5))
        fig.suptitle(uidname)
        ax = fig.add_subplot(131)
        scale = np.abs(header['CDELT1'])*3600
        x_index = (np.arange(0, nx) - nx/2.0) * scale
        y_index = (np.arange(0, ny) - ny/2.0) * scale
        x_map, y_map = np.meshgrid(x_index, y_index)
        ax.pcolormesh(x_map, y_map, masked_data, vmin=3*lower_1sigma, vmax=5*upper_1sigma)
        ax.text(0, 0, '+', color='r', fontsize=24, fontweight=100, horizontalalignment='center',
                verticalalignment='center')
        circle = patches.Circle((0, 0), radius=bmaj*3600*radius*0.5, facecolor=None, fill=None, 
                                edgecolor='red', linewidth=2, alpha=0.5)
        ellipse = patches.Ellipse((0.8*np.min(x_index), 0.8*np.min(y_index)), width=bmin*3600, 
                                  height=bmaj*3600, angle=bpa, facecolor='orange', edgecolor=None, 
                                  alpha=0.8)
        ax.add_patch(circle)
        ax.add_patch(ellipse)
        ax.set_xlabel('RA [arcsec]')
        ax.set_ylabel('Dec [arcsec]')
        # show the image of the most central region
        ax = fig.add_subplot(132)
        ax.pcolormesh(x_map, y_map, masked_data)
        xlim, ylim = 40*scale, 40*scale
        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(-ylim, ylim)
        ellipse = patches.Ellipse((0, 0), width=bmin*3600, height=bmaj*3600, angle=bpa, fill=None, 
                                  facecolor=None, edgecolor='red', alpha=0.5)
        ax.add_patch(ellipse)
        ax.set_xlabel('RA [arcsec]')
        ax.set_ylabel('Dec [arcsec]')

        # niose statistics
        ax = fig.add_subplot(133)
        ax.step(bins_mid, hist/amp_scale, where='mid', color='b', label='Noise Distribution')
        ax.step(bins_mid, hist_center/amp_scale_center, where='mid', color='orange', 
                label='Central Noise Distribution')
        ax.plot(bins_mid, hist_fit/amp_scale, color='r', label='Gaussian Fitting')
        ax.vlines(upper_3sigma, 0, 2.0, color='k', label=r'3$\sigma$ boundary')
        ax.vlines(lower_3sigma, 0, 2.0, color='k')
        ax.vlines(5.*upper_1sigma, 0, 2.0, color='k', lw=4, label=r'5$\sigma$ upper boundary')
        ax.set_xlabel('Flux density [Jy/beam]')
        ax.set_ylabel('Normalized Pixel numbers')
        ax.tick_params(axis='x', which='major', labelsize=8)
        ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='x')
        ax.legend(loc=2, prop={'size': 6})
        # plt.tight_layout()
        if debug:
            plt.show()
        else:
            plt.close()
        if savefig:
            if not figname:
                figname = os.path.join(outdir, os.path.basename(img) + '.png')
            fig.savefig(figname, dpi=2*fig.dpi)
    
    # return the checking results
    # check the fiiting
    if np.abs(mean-0.0)<1e-8 and np.abs(sigma-0.001)<1e-8: #compare with initial guess
        if debug:
            print("Fitting failed!")
        return False
    if np.abs(median_central) > central_median_deviation*sigma:
        if debug:
            print("Rejected, large central median value!\n")
        return False
    # comparing the noise distribution with Gaussian
    for deviation in [deviation_1sigma, deviation_2sigma]:
        if deviation > gaussian_deviation:
            if debug:
                print("Rjected, non-Gaussian noise")
            return False
    if percent_outlier_central >= outlier_frac:
        if debug:
            print("Rejected, large central residual!\n")
        return False
    if percent_outlier_central_neg >= outlier_frac:
        if debug:
            print("Rejected, negative central residual!\n")
        return False
    if strick_mode:
        if deviation_3sigma > gaussian_deviation:
            if debug:
                print("Rejected, non-Gaussian noise at 3sigma boudary!\n")
            return False
        if mean_central > threshold_mean:
            if debug:
                print("Rejected, large central mean value!\n")
            return False
        if debug:
            print("\n")
    return True

def check_images(imgs, outdir=None, basename='', band=None, debug=False, overwrite=False,
        **kwargs):
    """wraps up check_image to handle multiple images
    """
    if outdir:
        goodfile = os.path.join(outdir, basename+"_good_imgs.txt")
        badfile = os.path.join(outdir, basename+"_bad_imgs.txt")
        if os.path.isfile(goodfile) or os.path.isfile(badfile):
            if not overwrite:
                print("Classification already done!")
                return [], []
    if isinstance(imgs, str):
        all_files = glob.glob(imgs)
    elif isinstance(imgs, list):
        all_files = imgs

    # p_uidimg = re.compile('(?P<uidname>uid___\w+.ms.split.cal.J\d+[+-]\d+_B\d+).*.fits')
    p_uidimg = re.compile('(?P<uidname>uid___\w*\.ms[\.split\.cal]*\.J\d*[+-]+\d*_B\d+).*\.fits')

    good_imgs = []
    bad_imgs = []
    if band:
        image_outdir = os.path.join(outdir, band)
    else:
        image_outdir = outdir
    os.system('mkdir -p {}'.format(image_outdir))
    for img in all_files:
        print("img: {}".format(img))
        # continue
        if p_uidimg.search(img):
            uidname = p_uidimg.search(img).groupdict()['uidname']
        else:
            uidname = None
        try: 
            if check_image(img, debug=debug, outdir=image_outdir, **kwargs):
                if uidname is not None:
                    good_imgs.append(uidname)
                else:
                    if debug:
                        print("Not support filenames: {}".format(img))
            else:
                if uidname is not None:
                    bad_imgs.append(uidname)
                else:
                    if debug:
                        print("{} is not include in the returns!".format(img))
        except:
            print("Faild check: {}".format(img))
            if uidname:
                bad_imgs.append(uidname)
    if outdir:
        if not os.path.isdir(outdir):
            os.system('mkdir -p {}'.format(outdir))
        if len(good_imgs) > 0:
            with open(os.path.join(outdir, basename+"_good_imgs.txt"), "w") as good_outfile:
                good_outfile.write("\n".join(str(item) for item in good_imgs))
        if len(bad_imgs) > 0:
            with open(os.path.join(outdir, basename+"_bad_imgs.txt"), "w") as bad_outfile:
                bad_outfile.write("\n".join(str(item) for item in bad_imgs))
    return good_imgs, bad_imgs

def check_images_manual(imagedir=None, goodfile=None, badfile=None, debug=False, ncol=1, nrow=3, suffix='.updated'):
    '''visual inspection of the classification results from check_images
    
    '''
    if imagedir is None:
        raise ValueError("basedir should be defined along with filelist")
    all_good_files = []
    all_bad_files = []
    has_goodfile = True
    has_badfile = True
    if goodfile:
        # if debug:
        print('imagedir', imagedir)
        print('goodfile', goodfile)
        if isinstance(goodfile, str):
            try:
                flist = []
                with open(goodfile) as f:
                    filelist_lines = f.readlines()
                for line in filelist_lines:
                    pngfile = line.strip()+'.cont.auto.image.fits.png'
                    all_good_files.append(os.path.join(imagedir, pngfile))
            except:
                print("Failed in open {}".format(goodfile))
                has_goodfile = False
                pass
    
    if badfile:
        if debug:
            print('badfile', badfile)
        if isinstance(badfile, str):
            try:
                with open(badfile) as f:
                    filelist_lines = f.readlines()
                for line in filelist_lines:
                    pngfile = line.strip()+'.cont.auto.image.fits.png'
                    all_bad_files.append(os.path.join(imagedir, pngfile))
            except:
                print("Failed in open {}".format(badfile))
                has_badfile = False
                pass
    if (not has_goodfile) and (not has_badfile):
        return

    list_patches = {}
    for desc,all_files in list(zip(['good', 'bad'], [all_good_files, all_bad_files])):
        print(">>>>>>>>>>> {} images".format(desc))
        is_continue = int(raw_input("Continue? File list has changed! \n ") or '1')
        if is_continue == 0:
            break
        total_num = len(all_files)
        select_num = 0
        print("Find {} files".format(total_num))
        all_select = []
        for i in range(0, len(all_files), ncol*nrow):
            fig = plt.figure(figsize=(12*ncol, 4*nrow))
            for j in range(0, nrow*ncol):
                if (i+j)>=len(all_files):
                    continue
                try:
                    imagedata = plt.imread(all_files[i+j])
                except:
                    print("Error in reading: {}".format(all_files[i+j]))
                    continue
                ax = fig.add_subplot(nrow, ncol, j+1)#, projection=wcs2, slices=(0, 0, 'x', 'y'))
                # ax.text(10, 10, str(j), fontsize=20)
                ax.set_title(str(j+1))
                ax.imshow(imagedata, interpolation='none')
                plt.tight_layout()
        
            # show the image and record the 
            plt.show()
            print('Input the index of images, seperate with comma: [{}/{}]'.format(i+3,total_num))
            try:
                find_zero = False
                idx_input = input()
                while idx_input < 0:
                    print("Previous selection: \n    {}".format(all_select[idx_input]))
                    print("Select again:")
                    idx_input = input()

                if idx_input == 0:
                    print("Currently at {}/{}".format(i, len(all_files)))
                    plt.close('all')
                    break
                if isinstance(idx_input, int):
                    idx_input = [idx_input]
                select_num += len(idx_input)
                for ind in idx_input:
                    if ind == 0:
                        print("Currently at {}/{}".format(i, len(all_files)))
                        plt.close('all')
                        find_zero = True
                        break
                    all_select.append(all_files[i+ind-1])
                if find_zero:
                    break
            except:
                plt.close('all')
                continue
            # plt.clf()
            plt.close('all')
        list_patches[desc] = all_select
        if total_num > 0:
            print("Totally {}% of data have been selected.".format(100.*select_num/total_num))
        else:
            print("No data found.")

    list_updated = {}
    list_updated['good'] = (set(all_good_files) - set(list_patches['good'])).union(set(list_patches['bad']))
    list_updated['bad'] = (set(all_bad_files) - set(list_patches['bad'])).union(set(list_patches['good']))
    
    if debug:
        print('list_patches')
        print(list_patches)
        print("list_updated")
        print(list_updated)
    else:
        # obsname_match = re.compile('(?P<obsname>uid___\w*\.ms\.split\.cal\.J\d*[+-]+\d*_B\d+)')
        # fix the problem with naming issue
        obsname_match = re.compile('(?P<obsname>uid___\w*\.ms[\.split\.cal]*\.J\d*[+-]+\d*_B\d+)')
        for desc, f in zip(['good', 'bad'], [goodfile, badfile]):
            # if not os.path.isfile(f):
                # print("Warning! {} doesn't exist!".format(f))
                # continue
            if len(list_updated[desc]) < 1:
                continue
            with open(f+suffix, 'w+') as f:
                for item in list_updated[desc]:
                    try:
                        obsname = obsname_match.search(item).groupdict()['obsname']
                        f.write(obsname+'\n')
                    except:
                        print("Error in matching the obs name for filname: {}".format(item))
                        continue

def check_images_manual_gui(imagedir=None, goodfile=None, badfile=None, debug=False, ncol=1, nrow=3):
    '''visual inspection of the classification results from check_images
    
    '''
        # gui program
    from Tkinter import Tk, Frame, Checkbutton, IntVar, Button, Canvas
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    #
    if imagedir is None:
        raise ValueError("basedir should be defined along with filelist")
    all_good_files = []
    all_bad_files = []
    if goodfile:
        if debug:
            print('goodfile', goodfile)
        if isinstance(goodfile, str):
            try:
                flist = []
                with open(goodfile) as f:
                    filelist_lines = f.readlines()
                for line in filelist_lines:
                    pngfile = line.strip()+'.cont.auto.fits.png'
                    all_good_files.append(os.path.join(imagedir, pngfile))
            except:
                print("Failed in open {}".format(goodfile))
                pass
    
    if badfile:
        if debug:
            print('badfile', badfile)
        if isinstance(badfile, str):
            try:
                with open(badfile) as f:
                    filelist_lines = f.readlines()
                for line in filelist_lines:
                    pngfile = line.strip()+'.cont.auto.fits.png'
                    all_bad_files.append(os.path.join(imagedir, pngfile))
            except:
                print("Failed in open {}".format(badfile))
                pass

    root = Tk()
    root.title('Checking images')
    # frame=Frame(root, width=300, height=300)
    # frame.pack(expand = True, fill=BOTH)
    # canvas = Canvas(frame, bg='white', width=300, height=300)


    # Var1 = IntVar()
    # Var2 = IntVar()
     
    # ChkBttn = Checkbutton(frame, width = 15, variable = Var1)
    # ChkBttn.pack(padx = 5, pady = 5)
     
    # ChkBttn2 = Checkbutton(frame, width = 15, variable = Var2)
    # ChkBttn2.pack(padx = 5, pady = 5)

    # Button = Button(frame, text = "Submit", command = retrieve)
    # Button.pack(padx = 5, pady = 5)

    list_patches = {}
    for desc,all_files in list(zip(['good', 'bad'], [all_good_files, all_bad_files])):
        print(">>>>>>>>>>> {} images".format(desc))
        is_continue = int(raw_input("Continue? File list has changed! \n ") or '1')
        if is_continue == 0:
            break
        total_num = len(all_files)
        select_num = 0
        print("Find {} files".format(total_num))
        all_select = []
        for i in range(0, len(all_files), ncol*nrow):
            fig = plt.figure(figsize=(12*ncol, 4*nrow))
            for j in range(0, nrow*ncol):
                if (i+j)>=len(all_files):
                    continue
                try:
                    imagedata = plt.imread(all_files[i+j])
                except:
                    print("Error in reading: {}".format(all_files[i+j]))
                    continue
                ax = fig.add_subplot(nrow, ncol, j+1)#, projection=wcs2, slices=(0, 0, 'x', 'y'))
                # ax.text(10, 10, str(j), fontsize=20)
                ax.set_title(str(j+1))
                ax.imshow(imagedata, interpolation='none')
                #plt.tight_layout()
                canvas = FigureCanvasTkAgg(f, master=root)
                canvas.show()
                button = Tk.Button(master=root, text='Quit', command=sys.exit)
                button.pack(side=Tk.BOTTOM)
                Tk.mainloop()
        
            # show the image and record the 
            # plt.show()
            # print('Input the index of images, seperate with comma: [{}/{}]'.format(i+3,total_num))
            # try:
                # find_zero = False
                # idx_input = input()
                # while idx_input < 0:
                    # print("Previous selection: \n    {}".format(all_select[idx_input]))
                    # print("Select again:")
                    # idx_input = input()

                # if idx_input == 0:
                    # print("Currently at {}/{}".format(i, len(all_files)))
                    # plt.close('all')
                    # break
                # if isinstance(idx_input, int):
                    # idx_input = [idx_input]
                # select_num += len(idx_input)
                # for ind in idx_input:
                    # if ind == 0:
                        # print("Currently at {}/{}".format(i, len(all_files)))
                        # plt.close('all')
                        # find_zero = True
                        # break
                    # all_select.append(all_files[i+ind-1])
                # if find_zero:
                    # break
            # except:
                # plt.close('all')
                # continue
            # # plt.clf()
            # plt.close('all')
        # list_patches[desc] = all_select
        # if total_num > 0:
            # print("Totally {}% of data have been selected.".format(100.*select_num/total_num))
        # else:
            # print("No data found.")

    list_updated = {}
    list_updated['good'] = (set(all_good_files) - set(list_patches['good'])).union(set(list_patches['bad']))
    list_updated['bad'] = (set(all_bad_files) - set(list_patches['bad'])).union(set(list_patches['good']))
    
    if debug:
        print('list_patches')
        print(list_patches)
        print("list_updated")
        print(list_updated)
    else:
        obsname_match = re.compile('(?P<obsname>uid___\w*\.ms\.split\.cal\.J\d*[+-]+\d*_B\d+)')
        for desc, f in zip(['good', 'bad'], [goodfile, badfile]):
            # if f is None:
                # continue
            with open(f+'.updated', 'w+') as f:
                for item in list_updated[desc]:
                    try:
                        obsname = obsname_match.search(item).groupdict()['obsname']
                        f.write(obsname+'\n')
                    except:
                        print("Error in matching the obs name for filname: {}".format(item))
                        continue

def vis_clone(vis, outdir='./', drop_FDM=True, computwt=True):
    """helper for make_all_goodimages, make clone visibility in parallel mode"""
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
    spw = ''
    if drop_FDM:
        spw_select = []
        tb.open(os.path.join(vis, 'SPECTRAL_WINDOW'))
        chan_num = tb.getcol('NUM_CHAN')
        for s,sn in enumerate(chan_num):
            if sn > 5:
                spw_select.append(str(s))
        spw = ','.join(spw_select)

    v_new = os.path.join(outdir, os.path.basename(vis))
    if os.path.isdir(v_new):
        print('Found splitted file: {}'.format(v_new))
        return v_new
    # Only use valid data
    intent_list = check_intent(vis)
    intent_list_full = []
    for i in intent_list:
        intent_list_full.append('CALIBRATE_{}#ON_SOURCE'.format(i))
    intent = ','.join(intent_list_full)
    # os.system('rm -rf {}*'.format(v_new))
    print('Splitting... {} {}'.format(vis, v_new))
    #os.system('cp -r {} {}'.format(v, v_new))
    if not split(vis=vis, outputvis=v_new, spw=spw, datacolumn='corrected', intent=intent):
        print("No corrected data existing, spliting the data column")
        if not split(vis=vis, outputvis=v_new, spw=spw, datacolumn='data', intent=intent):
            print("Data may be corrupted: {}".format(vis))
            os.system('rm -rf {}'.format(v_new))
            return None
    if computwt:
        statwt(v_new, datacolumn='data')
    return v_new

def make_good_image(vis=None, basename='', basedir=None, outdir='./', concatvis=None, debug=False, 
                    only_fits=True, niter=1000, clean=True, pblimit=-0.01, fov_scale=2.0, 
                    computwt=True, drop_FDM=True, save_psf=True, check_sensitivity=True,
                    uvtaper_list=['0.3arcsec', '0.6arcsec'], make_jackknif=False,
                    uvtaper_scale=None,#[1.5, 2.0], 
                    **kwargs):
    """make the final good image with all the good observations

    uvtaper_list: [['0.3arcsec'], ['0.8arcsec']], outer taper width
    uvtaper_scale: [1.0, 1.7], 1.5 and 2 times and three times worse
    """
    if len(vis) < 1:
        return False
    n_select = len(vis)
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
    if basedir:
        vis_fullpath = []
        for v in vis:
            vis_fullpath.append(os.path.join(basedir, v))
        vis = vis_fullpath
    if debug:
        print(vis)
    if concatvis is None:
        concatvis = os.path.join(outdir, basename+'.ms')
    # map_input = list(zip(vis, [outdir]*n_select))
    # print(map_input)
    # p = mp.Pool(processes=n_core)
    # vis_cloned = p.map(vis_clone, map_input)
    # p.close()
    # print(vis_cloned)
    # return
    vis_cloned = []
    for v in vis:
        try:
            v_cloned = vis_clone(v, outdir=outdir, computwt=computwt, drop_FDM=drop_FDM)
            if v_cloned is not None:
                vis_cloned.append(v_cloned)
        except:
            if debug:
                print("Error in clone: {}".format(v))
            continue
    vis = vis_cloned
        
    if check_sensitivity:
        vis_valid = []
        vis_excluded = []
        # if check_sensitivity:
            # p = mp.Pool(processes=n_core)
            # print('vis:', vis)
            # check_vis2image(vis[0])
            # p.map(check_vis2image, vis)
            # print(is_usable)
            # p.close()
        for v in vis:
            if check_sensitivity:
                is_usable = False
                imagefile = gen_image(vis=v, outdir=os.path.join(outdir, 'images'), suffix='.cont.auto')
                print("imagefile", imagefile)
                print('vis', v)
                try:
                    is_usable = check_vis2image(vis=v, imagefile=imagefile)
                except:
                    print("Warning: error in check {}".format(v))
                if is_usable:
                   vis_valid.append(v)
                else:
                   vis_excluded.append(v)
        if only_fits:
            os.system('rm -rf {}'.format(os.path.join(outdir, 'images')))
        print('After sensitivity check, v:', vis_valid)
        vis = vis_valid
        # print(vis_valid)
        with open(concatvis+'.included.txt', 'w') as f:
            f.write('\n'.join(str(item) for item in vis))
        with open(concatvis+'.excluded.txt', 'w') as f:
            f.write('\n'.join(str(item) for item in vis_excluded))

    if len(vis) < 1:
        print("No valid visibility left!")
        os.system('touch {}'.format(os.path.join(outdir, 'Done')))
        # clean the folder
        rmtables(os.path.join(outdir,'uid___*'))
        os.system('rm -rf {}'.format(os.path.join(outdir,'uid___*.flagversions')))
        return 0

    if os.path.isdir(concatvis):
        print("Skip concating, file exists....")
    else:
        concat(vis=vis, concatvis=concatvis)

    if uvtaper_list:
        for uvtaper in uvtaper_list:
            if uvtaper == '':
                uvtaper_name = ''
                myuvtaper = None
            else:
                uvtaper_name = '.'+uvtaper
                myuvtaper = uvtaper
            myimagename = concatvis+'.auto.cont{}'.format(uvtaper_name)
            make_cont_img(vis=concatvis, myimagename=myimagename, clean=clean, niter=niter, 
                          pblimit=pblimit, fov_scale=fov_scale, myuvtaper=myuvtaper, debug=debug, 
                          **kwargs)
    if uvtaper_scale:
            myimagename = concatvis+'.auto.cont'
            make_cont_img(vis=concatvis, myimagename=myimagename, clean=clean, niter=niter, 
                          pblimit=pblimit, fov_scale=fov_scale, uvtaper_scale=uvtaper_scale, 
                          debug=debug, **kwargs)
    if make_jackknif:
        vis_jackkniffed = vis_jackknif(concatvis)
        if uvtaper_scale:
                myimagename = concatvis+'.auto.cont.jackknif'
                make_cont_img(vis=vis_jackkniffed, myimagename=myimagename, clean=clean, niter=niter, 
                              pblimit=pblimit, fov_scale=fov_scale, uvtaper_scale=uvtaper_scale, 
                              debug=debug, **kwargs)

        if uvtaper_list:
            for uvtaper in uvtaper_list:
                if uvtaper == '':
                    uvtaper_name = ''
                    myuvtaper = None
                else:
                    uvtaper_name = '.'+uvtaper
                    myuvtaper = uvtaper
                myimagename = concatvis+'.auto.cont{}.jackknif'.format(uvtaper_name)
                make_cont_img(vis=vis_jackkniffed, myimagename=myimagename, clean=clean, niter=niter, 
                              pblimit=pblimit, fov_scale=fov_scale, myuvtaper=myuvtaper, debug=debug, 
                              **kwargs)
    if only_fits:
        for i in glob.glob(concatvis+'*.image'):
            exportfits(imagename=i, fitsimage=i+'.fits')
        if save_psf:
            for psf in glob.glob(concatvis+'*.psf'):
                exportfits(imagename=psf, fitsimage=psf+'.fits')

        rmtables(concatvis+'*')
        rmtables(os.path.join(outdir,'uid___*'))
        os.system('rm -rf {}'.format(os.path.join(outdir,'uid___*.flagversions')))

def calculate_effectivearea(flux=np.linspace(0.1, 1, 10), snr_threshold=5.0, images=None, 
        images_pbcorr=None, fovscale=1.5, central_mask_radius=1.0):
    """calculate the effective area for given images

    rlimit: in arcsec
    """
    if len(images) != len(images_pbcorr):
        raise ValueError('images and images_pbcorr should contain same sources!')
    n_fields = len(images)
    effarea = np.zeros_like(flux)
    for i in range(n_fields):
        with fits.open(images[i]) as hdu:
            data = hdu[0].data
            header = hdu[0].header
        with fits.open(images_pbcorr[i]) as hdu_pbcor:
            data_pbcor = hdu_pbcor[0].data
        pixel2arcsec = np.abs(header['CDELT1'])*3600
        pix2area = pixel2arcsec**2  # pixel to arcsec^2
        freq = header['CRVAL3']
        lam = (const.c/(freq*u.Hz)).decompose().to(u.um)
        fov = 1.02 * (lam / (12*u.m)).decompose()* 206264.806
        ny, nx = data.shape[-2:]
        x = (np.arange(0, nx) - nx/2.0) * pixel2arcsec
        y = (np.arange(0, ny) - ny/2.0) * pixel2arcsec
        r = np.sqrt(x**2 + y**2)
        a, b = header['BMAJ']*3600, header['BMIN']*3600
        # print('pixel2arcsec', pixel2arcsec)
        # print('r range', np.min(r), np.max(r))
        # print('a',a,'b',b)
        data_masked = np.ma.masked_invalid(data.reshape(ny, nx))
        # mean, median, std = sigma_clipped_stats(data_masked, sigma=5.0, iters=5)  
        mean, median, std = calculate_image_sensitivity(images[i])
        # std = imstat(images[sigma['sigma'][0]
        pbcor = (data / data_pbcor).reshape(ny, nx)
        if std < 1e-7:
            print("Suspicious image: {}".format(images[i]))
        mask = data_masked.mask & (data_masked.data >= 5.0*std)
        pbcor[mask] = 0.0
        if fovscale:
            rlimit = 0.5 * fovscale * fov.value
            print("rlimit={}, maximal radius={}".format(rlimit, ny/2.0*pixel2arcsec))
            if rlimit < 0.8*ny/2.0*pixel2arcsec:
                print("Warning: checking the fovscale!.....")
            pbcor[r > rlimit] = 0.0
        if central_mask_radius > 1e-8:
            pbcor[r < central_mask_radius*a] = 0.0
        effarea_list = []
        for f in flux:
            snr = f / (std * 1000) # from Jy to mJy
            snr_map = snr * pbcor
            #print('selected number of pixels', np.sum(snr_map > snr_threshold))
            area = np.sum(snr_map > snr_threshold) * pix2area
            effarea_list.append(area)
        effarea = effarea + np.array(effarea_list) / 3600. # convert to arcmin^2
    return np.array([flux, effarea])

def read_fluxsummary(basedir, obj, band, resolution,):
    """read flux from the obj summary file
    """
    summary_file = os.path.join(basedir, obj, obj+'.sources_found.txt')
    with open(summary_file) as sf:
        lines = sf.readlines()
        entry = '# {} {} {}\n'.format(obj, band, resolution)
        if entry in lines:
            idx = lines.index(entry)
            flux = (np.array(map(np.float, lines[idx+1].split()))[2:])
            flux = flux.reshape((len(flux)/2, 2))
        else:
            flux = None
    return flux

def search_band_detection(basedir=None, band='B3', outdir='./', debug=False, **kwargs):
    """This function used to search source detections in other bands

    basedir: the obj directory contains the data of all the bands
    band: the target band
    """
    if os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
    copy_ms(basedir=basedir, outdir=outdir, select_band=band, debug=debug)
    image_dir = os.path.join(outdir, 'images')
    goodimage_dir = os.path.join(outdir, 'goodimages')
    gen_images(dirname=outdir, debug=debug, outdir=image_dir)
    basename = os.path.basename(basedir)+'_{}'.format(band)
    
    is_visual_checking = str(raw_input("Checking the images now? [n/y]: ") or 'n')
    if is_visual_checking == 'y':
        # selected_files = os.path.join(image_dir, 'seleted.txt')
        # show_images(fileglob=image_dir+'/*.fits', savefile=selected_files)
        check_images(os.path.join(image_dir,'*.fits'), outdir=goodimage_dir, band=band, 
                basename=basename, plot=True, savefig=True)
        band_imagedir = os.path.join(goodimage_dir, band)
        goodfile = os.path.join(goodimage_dir, "{}_good_imgs.txt".format(basename))
        badfile = os.path.join(goodimage_dir, "{}_bad_imgs.txt".format(basename))
        check_images_manual(imagedir=band_imagedir, goodfile=goodfile, badfile=badfile, 
                ncol=1, nrow=3)
    else:
        print("See you next time...")
        return
    goodfile_updated = os.path.join(goodimage_dir, "{}_good_imgs.txt.updated".format(basename))
    print(goodfile_updated)
    if os.path.isfile(goodfile_updated):
        vis_list = gen_filenames(listfile=goodfile_updated, basedir=outdir)
        make_good_image(vis=vis_list, basename=basename+'_combined', outdir=goodimage_dir, **kwargs)

def combine_sources(coords_list, units=u.arcsec, tolerance=0.3*u.arcsec):
    """This function used to combine the position of the detections,
       which come out from the images with different resolution,
       thus their coordinates could have small offset

       The program will prefer the first one when multiple coordinates 
    """
    coords_unique = [coords_list[0]]
    for coord1 in coords_list[1:]:
        is_unique = 1
        for coord2 in np.array(coords_unique):
            coord_difference = np.array(coord1)*units - np.array(coord2)*units
            # print(np.sqrt(coord_difference[0]**2 + coord_difference[1]**2))
            if coord_difference[0]**2 + coord_difference[1]**2 < tolerance**2:
                is_unique = False
        if is_unique:
            coords_unique.append(coord1)
    return coords_unique

def mock_observation(image=None, image_pbcorr=None, radius=None, max_radius=16, step=1.0,
        fNN=None, snr_threshold=5.0, fovscale=2.0, central_mask_radius=1.0):
    """
    images: the real images accounting primary beam response
    fNN: callable, the cumulative function give the number of detections at a give sensitivity and area
    """
    if radius is None:
        radius = np.arange(central_mask_radius, max_radius, step)
    with fits.open(image) as hdu:
        data = hdu[0].data
        header = hdu[0].header
    with fits.open(image_pbcorr) as hdu_pbcor:
        data_pbcor = hdu_pbcor[0].data
    pixel2arcsec = np.abs(header['CDELT1'])*3600
    pix2area = pixel2arcsec**2  # pixel to arcsec^2
    deg2pixel = 1.0/np.abs(header['CDELT1'])
    #print('pixel2arcsec', pixel2arcsec)
    freq = header['CRVAL3']
    lam = (const.c/(freq*u.Hz)).decompose().to(u.um)
    fov = 1.02 * (lam / (12*u.m)).decompose()* 206264.806 # to arcsec2
    ny, nx = data.shape[-2:]
    x = (np.arange(0, nx) - nx/2.0) * pixel2arcsec
    y = (np.arange(0, ny) - ny/2.0) * pixel2arcsec
    x_map, y_map = np.meshgrid(x, y)
    #print(np.min(x), np.max(x))
    r = np.sqrt(x_map**2 + y_map**2)
    #a, b = header['BMAJ']*3600, header['BMIN']*3600
    a, b = header['BMAJ']*deg2pixel, header['BMIN']*deg2pixel
    beamsize = np.pi*a*b/(4*np.log(2))
 
    data_masked = np.ma.masked_invalid(data.reshape(ny, nx))
    mean, median, std = calculate_image_sensitivity(image)
    pbcor = (data / data_pbcor).reshape(ny, nx)
    if std < 1e-7:
        print("Suspicious image: {}".format(image))
    mask = data_masked.mask #& (data_masked.data >= 5.0*std)
    pbcor[mask] = 0.0
    if fovscale:
        rlimit = 0.5 * fovscale * fov.value
        pbcor[r > rlimit] = 0.0
    if central_mask_radius > 1e-8:
        pbcor[r < central_mask_radius*a] = 0.0
    Nr = []
    # std_map = std / (pbcor + 1e-8)
    # print('r map shape:', r.shape)
    # print('totally pixels:', np.sum(r > 0))
    n_radii = len(radius)
    for i in range(n_radii - 1):
        r_lower = radius[i]
        r_upper = radius[i+1]
        r_select = (r > r_lower) & (r < r_upper) 
        #r_select = r < r_upper 
        area = np.pi*(r_upper**2-r_lower**2)
        area2 = np.sum(r_select) * pix2area
        # print('area1:', area, "area2:", area2)
        #print('pixels:', np.sum(r_select))
        pb_select = r_select & (pbcor > 1e-8)
        if len(pbcor[pb_select]) > 0:
            pb_mean = np.mean(pbcor[pb_select])
        else:
            pb_mean = 1e-8
        # print("pb_mean", pb_mean)
        # print('beam size', beamsize)
        # flux_sensitivity = snr_threshold * std * 1000 / beamsize / pb_mean #change to mJy
        flux_sensitivity = snr_threshold * std * 1000 / pb_mean #change to mJy
        Nr.append(fNN(flux_sensitivity, area))
        # print('std', std, 'pb_median', pb_mean)
        # print('flux', flux_sensitivity, 'r_lower', r_lower, 'r_upper', r_upper, 'area', area, 'area2', np.pi*(r_upper**2-r_lower**2))
    return Nr

def flux_scaling(wavelength, target_wavelength, flux=1, SED_template='composite_sed_Lir_01_09.dat', z=2.0,
                      debug=False):
    """interpolate the flux at a specified wavelength
    
    Params:
        wavelength : in um
        flux : any units, self consistent
    """
    if 'ALMACAL_NUMBERCOUNTS_HOME' in os.environ.keys():
        root_path = os.environ['ALMACAL_NUMBERCOUNTS_HOME']
    else:
        root_path = os.path.join(os.path.expanduser('~'), 'work/projects/almacal/number_counts')
    SED_template_file = os.path.join(root_path, 'code/data/', SED_template)
    SED = np.loadtxt(SED_template_file)
    f_SED = interpolate.interp1d(SED[:,0]*(1.0+z), SED[:,1]/(1+z)**4)
    scale_factor = flux / f_SED(wavelength)
    if debug:
        print("flux at {}: {}".format(wavelength, f_SED(wavelength)))
        print("flux at {}: {}".format(target_wavelength, f_SED(target_wavelength)))
        print("scale_factor: {}".format(scale_factor))
        print("flux difference: {}".format(f_SED(wavelength)/f_SED(target_wavelength)))
    return scale_factor * f_SED(target_wavelength)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# The ALMA run automatic pipeline section #
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def run_line_search(basedir=None, almacal_z=None, zrange=None, lines=None, debug=False, savefile=None, 
                    v_tol=500, strick_mode=False):
    """
    Args:
        lines (dict):in units of GHz, example: {'CO1-0':115.3, 'CO2-1':230.5}
        v_tol: velocity tolerance, units in km/s
    """
    c = 299792.458 # speed of light, in km/s
    if basedir is None:
        raise ValueError("ALMACAL database should be specified with basedir!")
    if almacal_z is None:
        almacal_z = os.path.join(os.path.expanduser('~'), 'Documents/projects/almacal/data/almacal_z.txt')
    almacal_coords = Table.read(almacal_z, format='csv')

    if zrange is None:
        zrange = [0., np.inf]
    zselect = np.bitwise_and(almacal_coords['zCal']>zrange[0], almacal_coords['zCal']<zrange[1]) 

    almacal_zselect = almacal_coords[zselect]

    source_list = []
    base_dir = ''

    #ALMA band information, in GHz
    band_list = {'B3':[84, 116], 'B4':[125, 163], 'B5':[163, 211], 
                 'B6':[211, 275], 'B7':[275, 373], 'B8':[385, 500],
                 'B9':[602, 720], 'B10':[787, 950]}


    n_select = len(almacal_zselect)

    p_obs = re.compile('uid___')
    band_match = re.compile('_(?P<band>B\d{1,2})')

    searching_info = {}
    for name, freq in lines.items():
        searching_info[name] = []
        for i in range(n_select):
            obj_info = almacal_zselect[i]
            freq_observed = freq / (1.+obj_info['zCal'])
            freq_tol = v_tol/c * freq_observed
            freq_observed_range = [freq_observed-freq_tol, freq_observed+freq_tol]
            # print('Observed freq:', freq_observed)
            for band, band_freq in band_list.items():
                if freq_observed>band_freq[0] and freq_observed<band_freq[1]:
                    target_band = band
            obj = obj_info['#calName']
            obj_basedir = os.path.join(basedir, obj)
            if not os.path.isdir(obj_basedir):
                if debug:
                    print("skip {}".format(obj))
                continue
            if debug:
                print('>>>>>>>>>\n{}\n'.format(obj))
                print('z={}'.format(obj_info['zCal']))
            obj_searching_info = {'obj':obj}
            obj_total_time = 0
            for obs in os.listdir(obj_basedir):
                if p_obs.search(obs):
                    if band_match.search(obs):
                        obsband = band_match.search(obs).groupdict()['band']
                        if obsband != target_band:
                            continue
                        obs_fullname = os.path.join(obj_basedir, obs)
                        spw_list = read_spw(obs_fullname)
                        for spw in spw_list:
                            if strick_mode:
                                is_coveraged = (freq_observed_range[0]>spw[0]) and (freq_observed_range[1]<spw[1])
                            else:
                                is_coveraged = (freq_observed_range[1]>spw[0]) or (freq_observed_range[0]<spw[1])
                            if is_coveraged:
                                try:
                                    time_on_source = au.timeOnSource(obs_fullname, verbose=False, debug=False)
                                except:
                                    print("Faild when reading {}".format(obs_fullname))
                                    timeOnSource = 0
                                # total_time += time_on_source[0]['minutes_on_source']
                                obj_searching_info[obs] = time_on_source[0]['minutes_on_source']
                                obj_total_time += time_on_source[0]['minutes_on_source']
                                break
            if obj_total_time > 1e-8:
                obj_searching_info['total_time'] = obj_total_time
                searching_info[name].append(obj_searching_info)
            if debug:
                print(obj, obj_searching_info)
                print('>>>>>>>>>>>>\n\n')
    if savefile:
        with open(savefile, 'w+') as jf:
            json.dump(searching_info, jf)
    return searching_info

def run_gen_all_obstime(basedir=None, objs=None, outdir=None, bad_obs=None, 
        bands=['B3','B4','B5','B6','B7','B8','B9','B10'], info_file=None, 
        exclude_aca=True, 
        time_select=False, start_time='2010-01-01T00:00:00', end_time='2050-01-01T00:00:00', 
        debug=False, **kwargs):
    """generate the on-source time and spw distribution for the whole almacal
       dataset
    
    Params:
        base_dir: the root folder contains all the measurements
        output_dir: the root folder for placing the results of each calibrator,
                    including the json file and the plots
        bad_obs: the file contains unusable observation
        info_file: the file contains the general information for all the 
                   calibrators
        **kwargs: support the addition parameters of `spw_stat`
    

    default run:
    obj_output_dir = output_dir + '/' + obj + '/'
    os.system('mkdir {}'.format(obj_output_dir))
    spw_stat('/', plot=True, showfig=False, figname=obj_output_dir + obj+'.pdf', savedata=True, filename=obj_output_dir + obj+'.json')
    """

    p_obj = re.compile('J\d+[+-]\d')
    p_obs = re.compile('uid___')

    os.system('mkdir -p {}'.format(outdir))
    if bad_obs is not None:
        with open(bad_obs) as f:
            all_bad_obs = f.readlines()
        for i in range(len(all_bad_obs)):
            all_bad_obs[i] = all_bad_obs[i].strip()

    band_match = re.compile('_(?P<band>B\d{1,2})$')
    obj_match = re.compile('J\d{3,4}[-+]\d{3,4}')

    if objs is None:
        objs = []
        for item in os.listdir(basedir):
            if p_obj.match(item):
                objs.append(item)
            else:
                if debug:
                    print('Error load obj:', item)

    for obj in objs:
        obj_dir = os.path.join(basedir, obj)
        obj_outdir = os.path.join(outdir, obj)
        os.system('mkdir -p {}'.format(obj_outdir))
        vis_list = gen_filenames(dirname=obj_dir)
        obj_stat = spw_stat(vis=vis_list, debug=debug, exclude_aca=exclude_aca,
                            bands=bands, savefile=obj_outdir+'/'+ obj+'.json', 
                            time_select=time_select, start_time=start_time, end_time=end_time, 
                            **kwargs)
        
        if info_file is not None:
            string_stat = '{:<12s}'.format(obj)
            for band in obj_stat.keys():
                string_stat += ' {:>8.2f}'.format(np.sum(obj_stat[band]['time']))
            with open(info_file, 'a+') as f_info:
                f_info.write(string_stat + '\n')

def run_gen_stats_old(basedir=None, listfile_dir=None, objs=None, bands=['B6','B7'], 
        suffix=['good_imgs.txt.updated', 'bad_imgs.txt.updated'],
        exclude_aca=True, debug=False, pickle_file=None,
        func_list=[],
        **kwargs):
    """The function to make some statistics like gen_all_obstime in a more general way
    """

    p_obj = re.compile('J\d+[+-]\d')
    p_obs = re.compile('uid___')
    band_match = re.compile('_(?P<band>B\d{1,2})$')
    obj_match = re.compile('J\d{3,4}[-+]\d{3,4}')

    if objs is None:
        objs = []
        for item in os.listdir(basedir):
            if obj_match.match(item):
                objs.append(item)
            else:
                if debug:
                    print('Error load obj:', item)

    all_stats = {}
    for i,obj in enumerate(objs):
        obj_stat = {}
        print('index=', i, "obj:", obj)
            
        obj_dirname = os.path.join(basedir, obj)
        # obj_output_dir = os.path.join(output_dir, obj)
        # os.system('mkdir -p {}'.format(obj_output_dir))
        if listfile_dir:
            vis_list = []
            for suf in suffix:
                vis_list_suf = []
                for band in bands:
                    obj_band_selectfile = os.path.join(listfile_dir, obj, 
                            "{}_{}_{}".format(obj, band, suf))
                    if not os.path.isfile(obj_band_selectfile):
                        continue
                    obj_band_list = gen_filenames(listfile=obj_band_selectfile, basedir=obj_dirname, exclude_aca=exclude_aca)
                    vis_list_suf.append(obj_band_list)
                vis_list.append(flatten(vis_list_suf)) 
        else:
            vis_list = gen_filenames(dirname=obj_dirname)

        for func in func_list:
            obj_stat[func.func_name] = func(vis_list)
        all_stats[obj] = obj_stat
    if pickle_file:
        with open(pickle_file, 'wb') as f_pickle:
                pickle.dump(all_stats, f_pickle)

def run_gen_oteo2016_data(basedir, outdir, objs=None, select_band=['B6', 'B7'], debug=False):
    """copy the data used in oteo2016

    default run: run_gen_oteo2016_data('/science_ALMACAL/data', outdir='./')
    """
    p_obj = re.compile('J\d+[+-]\d+')
    if objs is None:
        objs = []
        for item in os.listdir(basedir):
            if p_obj.match(item):
                objs.append(item)
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
    for obj in objs:
        obj_input_folder = os.path.join(basedir, obj)
        obj_output_folder = os.path.join(outdir, obj)
        print("Copying {}".format(obj))
        copy_ms(obj_input_folder, obj_output_folder, time_select=True, 
                end_time='2015-07-01T00:00:00',
                select_band=select_band, debug=debug)

def run_gen_all_images(basedir, objlist=None, outdir='./', bands=['B6','B7'], exclude_aca=True, 
                  debug=False, niter=0, **kwargs):
    """fix the missing and wrong images for gen_all_image output
    

    default run: run_gen_all_image('basedir', outdir='', bands=[], )
    """
    filelist = []
    obj_match = re.compile('^J\d*[+-]\d*$')
    obs_match = re.compile('(?P<obsname>uid___\w*\.ms(\.split\.cal)?\.(?P<objname>[\s\w+-]+)_(?P<band>B\d+))')
    if objlist is None:
        objlist = []
        for obj in os.listdir(basedir):
            if obj_match.match(obj):
                objlist.append(obj)
    for obj in objlist:
        print(obj)
        # start to go through each obj

        for obs in os.listdir(os.path.join(basedir, obj)):
            if obs_match.match(obs):
                obs_band = obs_match.search(obs).groupdict()['band']
                infile = os.path.join(basedir, obj, obs)
                if debug:
                    print(obs)
                if obs_band not in bands:
                    continue
                outfile_fullpath = os.path.join(outdir, obs_band, obj)
                outfile_fullpath_invalid = os.path.join(outfile_fullpath, 'invalid')
                os.system('mkdir -p {}'.format(outfile_fullpath_invalid))
                
                outfile_fullname = os.path.join(outfile_fullpath, obs+'.cont.auto.fits') 
                try:
                    fitsimage = gen_image(infile, outdir=outfile_fullpath, check_validity=True, 
                                          debug=debug, niter=niter, **kwargs)
                    # if is_valid:
                        # os.system('mv {} {}'.format(outfile_fullname, os.path.join(outfile_fullpath, 'valid')))
                    # if not is_valid:
                        # os.system('mv {} {}/'.format(outfile_fullname, outfile_fullpath_invalid))
                except:
                    print("Warning:\n Error found in imaging {}\n".format(infile))
                    continue

def run_auto_classify_goodimags(imagedir=None, objlist=None, outdir='./', bands=[], 
        debug=False, suffix='_good_imgs.txt', plot=True, savefig=True, **kwargs):
    """generate the good image list for all the calibrators

    default run: run_make_all_goodimags(imgs_dir='all_img_dir', basedir='science_ALMACAL', outdir='./') 
    """
    if imagedir is None:
        raise ValueError("imagedir is required!")
    obj_match = re.compile('^J\d*[+-]\d*$')
    for band in bands:
        band_imagedir = os.path.join(imagedir, band)
        if not os.path.isdir(band_imagedir):
            print("Warning: no such derectory: {}".format(band_imagedir))
            continue
        band_outdir = os.path.join(outdir, band)
        os.system('mkdir -p {}'.format(band_outdir))
        if objlist is None:
            objlist = []
            objs_indir = os.listdir(band_imagedir)
            if len(objs_indir) < 1:
                contibue
            for obj in objs_indir:
                if obj_match.match(obj):
                    objlist.append(obj)
        for obj in objlist:
            print('obj', obj)
            obj_dir = os.path.join(band_outdir, obj)
            if not os.path.isdir(obj_dir):
                os.system('mkdir -p {}'.format(os.path.join(band_outdir, obj)))
            obj_path = os.path.join(band_imagedir, obj)
            # print(obj_path+'/*.fits')
            good_imgs, bad_imgs = check_images(obj_path+'/*.fits', 
                    outdir=os.path.join(band_outdir, obj), plot=plot, savefig=savefig, 
                    band=band, basename=obj+'_'+band, debug=debug, **kwargs)
            if debug:
                print(good_imgs)

def run_manual_inspection(classifiedfolder=None, objlist=None, bands=['B6','B7'], 
        suffix='imgs.txt'):
    """manually checking all the images
    """
    if isinstance(objlist, str):
        if os.path.isfile(objlist):
            objlist = gen_filenames(listfile=objlist)
    elif isinstance(objlist, (list, np.ndarray)):
        pass
    else:
        objlist = []
        obj_match = re.compile('^J\d*[+-]\d*$')
        for band in bands:
            for obj in os.listdir(os.path.join(classifiedfolder, band)):
                if obj_match.match(obj):
                    objlist.append(obj)
        objlist = np.unique(objlist)
        #raise ValueError("Unsupported objlist!")
    for obj in objlist:
        for band in bands:
            obj_classified_folder = os.path.join(classifiedfolder, band, obj)
            obj_classified_imagedir = os.path.join(obj_classified_folder, band)
            goodfile = os.path.join(obj_classified_folder, "{}_{}_good_{}".format(obj, band, suffix))
            badfile = os.path.join(obj_classified_folder, "{}_{}_bad_{}".format(obj, band, suffix))

            if os.path.isfile(goodfile+'.updated') or os.path.isfile(badfile+'.updated'):
                print('{} of {} already done'.format(band, obj))
                continue
            else:
                print(obj_classified_imagedir)
                print(">goodfile: {}\n>badfile: {}".format(goodfile, badfile))
                check_images_manual(imagedir=obj_classified_imagedir, goodfile=goodfile, 
                                    badfile=badfile, debug=False, ncol=1, nrow=3)

def run_make_all_goodimages(classifiedfolder=None, objlist=None, bands=['B6','B7'], basedir=None, 
        outdir='./', debug=False, only_fits=True, update=True, imagefile_suffix='combine',
        computwt=True, listfile_suffix='good_imgs.txt.updated', overwrite=False, 
        make_jackknif=False, **kwargs):
    """generate the good image with updated list

    default run: run_make_all_goodimags(imgs_dir='all_img_dir', basedir='science_ALMACAL', outdir='./make_good_image') 
    """
    obj_match = re.compile('^J\d*[+-]\d*$')
    for band in bands:
        band_indir = os.path.join(classifiedfolder, band)
        band_outdir = os.path.join(outdir, band)
        if objlist is None:
            objlist = []
            for obj in os.listdir(band_indir):
                if obj_match.match(obj):
                    objlist.append(obj)
        for obj in objlist:
            obj_indir = os.path.join(band_indir, obj)
            obj_outdir = os.path.join(band_outdir, obj)
            print(obj)
            print(obj_outdir)
            if not os.path.isdir(obj_outdir):
                os.system('mkdir -p {}'.format(obj_outdir))
            good_image_file = os.path.join(obj_indir, "{}_{}_{}".format(obj, band, listfile_suffix))
            concatvis = os.path.join(obj_outdir, "{}_{}_{}.ms".format(obj, band, imagefile_suffix))
            if debug:
                print("good image file: {}\n concatvis_name: {}".format(good_image_file, concatvis))
            
            if os.path.isfile(good_image_file):
                good_image_fitsfile = concatvis+'.auto.cont.*image.fits'
                print('good_image_fitsfile: {}'.format(good_image_fitsfile))
                if os.path.isdir(obj_outdir):
                    uid_files = glob.glob(os.path.join(obj_outdir,'uid___*'))
                    if len(uid_files) > 0:
                        print("Found existing measurements!")
                        print(uid_files)
                        print("Pruning folder....")
                        rmtables(os.path.join(obj_outdir,'uid___*'))
                        os.system('rm -rf {}'.format(os.path.join(obj_outdir,'uid___*.flagversions')))
                        os.system('rm -rf {}'.format(os.path.join(obj_outdir,'uid___*.cfg')))
                if not overwrite:
                    if os.path.isfile(os.path.join(obj_outdir, 'Done')):
                        print("Skip {} of {}, delete the 'Done' file to recreate...".format(band,obj))
                        continue
                    if len(glob.glob(good_image_fitsfile)) > 0:
                        print("Skip {} of {}, delete the fits file to recreate...".format(band,obj))
                        continue
                combined_vis = gen_filenames(listfile=good_image_file)
                make_good_image(combined_vis, concatvis=concatvis, only_fits=only_fits, 
                        outdir=obj_outdir, computwt=computwt, make_jackknif=make_jackknif, 
                        basedir=os.path.join(basedir, obj), **kwargs)

def run_make_all_jeckknifimages():
    pass
    if true:
        # make_jackknif = True 
        # if jackknif:
            # good_image_jackknig_fitsfile = concatvis+'.auto.cont.*jackknif.image.fits'
            # if not overwrite and os.path.isfile(os.path.join(obj_outdir, 'Jackknif_Done')):
                # print("Skip {} of {}, delete the 'Jackknif_Done' file to continue...".format(band,obj))
                # make_jackknif = False
            # elif not overwrite and (len(glob.glob(good_image_jackknig_fitsfile)) > 0):
                # print("Skip {} of {}, delete the jackknif fits file to continue...".format(band,obj))
                # make_jackknif = False
        pass

def run_check_SMGs(basedir, objs=None, bands=['B6','B7'], suffix='combine.ms.auto.cont', 
                   resolutions=['0.3arcsec', '0.6arcsec'], ncol=1,
                   summary_file=None, view_mode='multiple', focus_bands=None,
                   interactive=True, outdir=None, continue_mode=True, 
                   methods=['aperture','peak'], **kwargs):

    """finding sources
    Adding a dot in the string: 
        resolutions = ['0.3arcsec','0.6arcsec']
    """
    if outdir is not None:
        if not os.path.isdir(outdir):
            os.system('mkdir -p {}'.format(outdir))
    objs_finished = []
    if focus_bands == None:
        focus_bands = bands
    if summary_file is not None:
        if continue_mode:
            if os.path.isfile(summary_file):
                summary = Table.read(summary_file, format='ascii')
                objs_finished = summary['obj']
            else:
                os.system('mv {} {}.old'.format(summary_file, summary_file))
        if not os.path.isfile(summary_file):
            print('Initializing the output file')
            with open(summary_file, 'w+') as f:
                f.write("obj")
                for band in focus_bands:
                    for res in resolutions:
                        f.write(' ')
                        f.write(band+'_'+res)
                for band in bands:
                    f.write(' detection_{} goodfield_{}'.format(band, band))
                f.write(' is_SMG')
                f.write(' is_RG')
                f.write(' is_Jet')
                f.write('\n')

    obj_match = re.compile('^J\d*[+-]\d*$')
    if objs is None:
        objs = []
        for band in focus_bands:
            for item in os.listdir(os.path.join(basedir, band)):
                if obj_match.match(item):
                    objs.append(item)
        objs = np.unique(objs).tolist()

    # star to loop over all the objs
    failed_files = []
    try:
        for obj in objs:
            if obj in objs_finished:
                print("{} already done!".format(obj))
                continue
            print(obj)
            if obj_match.match(obj):
                print('>>>>> {}'.format(obj))
                # write into summary file
                if outdir is not None:
                    obj_outdir = os.path.join(outdir, obj)
                    obj_summary_file = os.path.join(obj_outdir, '{}.sources_found.txt'.format(obj))
                    if not os.path.isdir(obj_outdir):
                        os.system('mkdir -p {}'.format(obj_outdir))
                    # open the summary file
                    obj_summary = open(obj_summary_file, 'w+')
                    summary_plot = os.path.join(obj_outdir, '{}.summary.png'.format(obj))
                else:
                    obj_summary = None
                    obj_outdir = None
                    summary_plot = None
                obj_sourcefound = {}
                obj_validbands = {}
                for band in focus_bands:
                    obj_validbands[band] = True
                    for res in resolutions:
                        obj_sourcefound[band+'_'+res] = []
                # make a summary plot
                nrow = int(1.0 * len(bands) / ncol)
                fig, ax = plt.subplots(nrow, len(resolutions)*ncol, 
                                       figsize=(ncol*4.2*len(resolutions),4*nrow)) 
                fig.suptitle(obj)
                #for img in imgs:
                # sources_number = {}
                # obj_sourcefound
                for i,band in enumerate(bands):
                    i_ax = i // ncol
                    obj_band_dir = os.path.join(basedir, band, obj)
                    for j,res in enumerate(resolutions):
                        j_ax = j + i%ncol*ncol
                        if len(bands) > 1:
                            ax_select = ax[i_ax,j_ax]
                        else:
                            ax_select = ax[max(i_ax,j_ax)]
                        if res == '':
                            res_string = ''
                        else:
                            res_string = res+'.'
                        image_name = "{}_{}_{}.{}image.fits".format(obj, band, suffix, res_string)
                        image_fullpath = os.path.join(obj_band_dir, image_name)
                        # print(image_fullpath)
                        #print('Finding source in:', image_fullpath)
                        if not os.path.isfile(image_fullpath):
                            if band in focus_bands:
                                obj_validbands['{}'.format(band)] = False
                            continue
                        savefile = image_name + '.source_found.txt'
                        figname = image_name + '.png'
                        sources_found = source_finder(image_fullpath, outdir=obj_outdir, 
                                ax=ax_select, pbcor=True, **kwargs)
                        ax_select.set_title('{} {}'.format(band, res))
                        #try:
                        #    sources_found = source_finder(image_fullpath, outdir=obj_outdir, 
                        #            ax=ax[i,j], pbcor=True)
                        #except:
                        #    print("Error found for {}".format(image_name))
                        #    failed_files.append(image_name)
                        if len(sources_found) > 0 and obj_summary is not None:
                            obj_summary.write('# {} {} {}\n'.format(obj, band, res))
                            obj_sourcefound['{}_{}'.format(band, res)] = sources_found
                            source_idx = 0
                            for ra, dec, flux, snr in sources_found:
                                obj_summary.write('{} {:.6f} {:.6f} '.format(source_idx, ra, dec))
                                for f_m, f_snr in zip(flux, snr):
                                    obj_summary.write(" {:.4f} {:.2f} ".format(f_m, f_snr))
                                obj_summary.write('\n')
                                source_idx += 1
                # write into files
                if obj_summary is not None:
                    obj_summary.close()
                found_string = obj
                for band in focus_bands:
                    for res in resolutions:
                        # print(obj_sourcefound[band+'_'+res])
                        found_string += ' '+str(len(obj_sourcefound[band+'_'+res]))
                # print(found_string)
                # save figure
                fig.subplots_adjust(wspace=0.2, hspace=0.2)
                if summary_plot:
                    fig.savefig(summary_plot, bbox_inches='tight', dpi=400)
                if summary_file: 
                    SMG_input = 0
                    RG_input = 0
                    Jet_input = 0
                    detections = {}
                    goodfields = {}
                    for band in bands:
                        goodfields[band] = 0
                        detections[band] = 0
                    if interactive:
                        plt.show()
                        print("Single Band:\n") 
                        print("Detection: 0)None +n)N Detections -1)Not Sure")
                        print("Usable: 0)No 1)Full -n)exclude central n*FWHM region")
                        print("General Classification")
                        print("Is SMG: 0)No 1)Yes +n)Multiple -1)Not Sure")
                        print("Is RG: 0)No 1)Yes +n)Multiple -1)Not Sure")
                        print("Is Jet: 0)No 1)Compact +n)Enlongated -1)Not Sure")
                        for band in focus_bands:
                            if not obj_validbands[band]:
                                continue
                            detection_input=int(raw_input("Dection in Band:{} (integer, 0,1,2,3,-1) [0]?: ".format(band)) or 0)
                            goodfield_input=int(raw_input("Usable for Band:{} (integer, 0,1,2,-1,-2) [0]?: ".format(band)) or 0)
                            detections[band] = detection_input
                            goodfields[band] = goodfield_input
                        if bands > 0:
                            SMG_input = int(raw_input("Is SMG? (integer, 0,1,2,-1) [0]: ") or 0)
                            RG_input = int(raw_input("Is RG? (integer, 0,1,2,-1) [0]: ") or 0)
                            Jet_input = int(raw_input("Is Jet? (integer, 0,1,2,-1) [0]: ") or 0)
                        plt.close()
                    for band in focus_bands:
                        found_string += ' {} {}'.format(detections[band], goodfields[band])
                    found_string += ' {}'.format(SMG_input)
                    found_string += ' {}'.format(RG_input)
                    found_string += ' {}'.format(Jet_input)
                    with open(summary_file, 'a+') as f:
                        f.write("{}\n".format(found_string)) 
                else:
                    next_one = int(raw_input("Next one [1/0] [1]: ") or 1)
                    if next_one == 1:
                        if view_mode == 'single':
                            plt.close()
                        elif view_mode == 'multiple':
                            continue
                    else:
                        return 0
                if (outdir is None) and (view_mode == 'multiple'):
                    is_close = int(raw_input("Close all windows? (1/0) [1]") or 1)
                    if is_close == 0:
                        return 0
                plt.close()
    except KeyboardInterrupt:
        return 0

def run_gen_stats(basedir, imagedir=None, outdir=None, savefile='gen_stat.txt',
        suffix='combine.ms.included.txt', debug=False,
        bands=['B3','B4','B5','B6','B7','B8','B9','B10'], **kwargs):
    """This function used to generate the spw statistics of final samle
    """
    if not os.path.isdir(outdir):
        os.system('mkdir {}'.format(outdir))

    band_list = {'B3':[84, 116], 'B4':[125, 163], 'B5':[163, 211], 
            'B6':[211, 275], 'B7':[275, 373], 'B8':[385, 500], \
            'B9':[602, 720], 'B10':[787, 950]}

    p_obj = re.compile('J\d+[+-]\d')
    p_obs = re.compile('uid___')

    band_match = re.compile('_(?P<band>B\d{1,2})$')
    obj_match = re.compile('J\d{3,4}[-+]\d{3,4}')

    for band in bands:
        print("\nWorking on {}...\n".format(band))
        band_imagedir = os.path.join(imagedir, band)
        band_outdir = os.path.join(outdir, band)
        os.system('mkdir -p {}'.format(band_outdir))
        band_summary_file = os.path.join(band_outdir, '{}_obstime.txt'.format(band))
        
        if os.path.isfile(band_summary_file):
            print("Found summary_file, delete it for overwriting!")
            continue

        band_objs = []
        for item in os.listdir(band_imagedir):
            if obj_match.match(item):
                band_objs.append(item)
            else:
                if debug:
                    print('Error load obj:', item)

        for i,obj in enumerate(band_objs):
            obj_basedir = os.path.join(basedir, obj)
            obj_exptime = {}
            print('index=', i, "obj:", obj)
                
            obj_imagedir = os.path.join(band_imagedir, obj)
            
            obj_band_selectfile = os.path.join(obj_imagedir, "{}_{}_{}".format(obj, band, suffix))
            if os.path.isfile(obj_band_selectfile):
                band_obj_list = gen_filenames(listfile=obj_band_selectfile, basedir=obj_basedir)
                obj_savefile = os.path.join(band_outdir, obj+'.json')
                if os.path.isfile(obj_savefile):
                    print("gen_stats file already exists, delete it for overwrite!")
                    obj_stat = spw_stat(jsonfile=obj_savefile, bands=[band,])
                else:
                    obj_stat = spw_stat(vis=band_obj_list, debug=debug,
                                bands=[band,], savefile=obj_savefile, 
                                **kwargs)
                if band_summary_file is not None: 
                    string_stat = '{:<12s} {:>8.2f}'.format(obj, np.sum(obj_stat[band]['time']))
                    with open(band_summary_file, 'a+') as f_info:
                        f_info.write(string_stat + '\n')

def run_measure_flux(basedir, objs=None, bands=['B6','B7'], focus_band=None,
                   suffix='combine.ms.auto.cont', target_wave=None,
                   resolutions=['0.3arcsec', '0.6arcsec'], 
                   selected_resolution='0.3arcsec', ncol=1, 
                   summary_file=None, view_mode='multiple',
                   continue_mode=True, **kwargs):
    """finding sources
    Adding a dot in the string: 
        resolutions = ['0.3arcsec','0.6arcsec']
    """
    objs_finished = []
    if summary_file is not None:
        if continue_mode:
            if os.path.isfile(summary_file):
                summary = Table.read(summary_file, format='ascii')
                objs_finished = summary['obj']
            else:
                os.system('mv {} {}.old'.format(summary_file, summary_file))
        if not os.path.isfile(summary_file):
            print('Initializing the output file')
            with open(summary_file, 'w+') as f:
                f.write("obj idx type ra dec flux_aperture flux_snr_aperture flux_gaussian flux_snr_gaussian flux_peak flux_snr_peak radial_distance")
                f.write('\n')

    obj_match = re.compile('^J\d*[+-]\d*$')
    if objs is None:
        objs = []
        for item in os.listdir(os.path.join(basedir, focus_band)):
            if obj_match.match(obj):
                objs.append(item)
        objs = np.unique(objs).tolist()

    # star to loop over all the objs
    failed_files = []
    try:
        for obj in objs:
            if obj in objs_finished:
                print("{} already done!".format(obj))
                continue
            print(obj)
            obj_focus_band_dir = os.path.join(basedir, focus_band, obj)
            obj_sourcefound = {}
            if obj_match.match(obj):
                print('>>>>> {}'.format(obj))

                # make a summary plot
                nrow = int(1.0 * len(bands) / ncol)
                fig, ax = plt.subplots(nrow, len(resolutions)*ncol, 
                                       figsize=(ncol*4.2*len(resolutions),4*nrow)) 
                fig.suptitle(obj)
                #for img in imgs:
                # sources_number = {}
                for i,band in enumerate(bands):
                    i_ax = i // ncol
                    obj_band_dir = os.path.join(basedir, band, obj)
                    for j,res in enumerate(resolutions):
                        j_ax = j + i%ncol*ncol
                        if len(bands) > 1:
                            ax_select = ax[i_ax,j_ax]
                        else:
                            ax_select = ax[max(i_ax,j_ax)]
                        if res == '':
                            res_string = ''
                        else:
                            res_string = res+'.'
                        image_name = "{}_{}_{}.{}image.fits".format(obj, band, suffix, res_string)
                        image_fullpath = os.path.join(obj_band_dir, image_name)
                        print(image_fullpath)
                        #print('Finding source in:', image_fullpath)
                        if not os.path.isfile(image_fullpath):
                            continue
                        ax_select.set_title('{}'.format(res))
                        sources_found = source_finder(image_fullpath,
                                methods=['single_aperture', 'peak'],
                                ax=ax_select, pbcor=True, central_mask_radius=2.0, **kwargs)
                        if band == focus_band:
                            obj_sourcefound['{}'.format(res)] = sources_found

                # interactive selecting sources
                if summary_file:
                    detections_summary = open(summary_file, 'a+')
                    coords_list = []
                    is_detection = int(raw_input("Dection in band: {} (integer, 0/1) [1]?: ".format(focus_band)) or 1)
                    # is_detection = 1
                    if is_detection > 0:
                        for res in resolutions:
                            detection_idx = str(
                                    raw_input("The index for the true detections of band {} in reselection:{}\nSeperate with comma: ".format(
                                               focus_band, res)))
                            for idx in detection_idx.split(','):
                                if idx == '':
                                    continue
                                ra, dec, flux, snr = obj_sourcefound['{}'.format(res)][int(idx)]
                                coords_list.append([ra, dec])
                        coords_unique = combine_sources(coords_list)
                        selected_image = "{}_{}_{}.{}.image.fits".format(obj, focus_band, suffix, selected_resolution)
                        seleted_image_fullpath = os.path.join(obj_focus_band_dir, selected_image)
                        print("Using {} for flux measurements.\n".format(seleted_image_fullpath))
                        sources_flux, sources_flux_snr, sources_radial_distance = flux_measure(
                                seleted_image_fullpath, coords_unique, target_wave=target_wave,
                                methods=['adaptive_aperture', 'gaussian', 'peak'], 
                                calculate_radial_distance=True)
                        # sources_flux, sources_flux_snr = sources_flux, sources_flux_snr
                        # print("sources_flux", sources_flux)
                        # print("sources_flux_snr", sources_flux_snr)
                        #print(sources_flux)
                        
                        # classify the detections: 0)not sure 1)SMG; 2)radio;
                        n_uniq = len(coords_unique)
                        source_type = []
                        print("Classify the detections: 0)not sure 1)SMG; 2)radio:")
                        for si in range(n_uniq):
                            print("For selected source: {}: flux={:.2f} snr={:.2f} dist={:.2f}".format(si, sources_flux[si][0], sources_flux_snr[si][0], 
                                                                                       sources_radial_distance[si]))
                        for si in range(n_uniq):
                            s_type = int(raw_input("Selected source: {} type:(integer, 0/1/2) [1]?: ".format(si)) or 1)
                            source_type.append(s_type)
                        print("Totoal {} sources".format(n_uniq))
                        for i in range(n_uniq):
                            detections_summary.write('{} {} {} {:.6f} {:.6f} {} {} {} {} {} {} {}'.format(obj,
                                                     i, source_type[i], coords_unique[i][0], coords_unique[i][1],
                                                     sources_flux[i][0], sources_flux_snr[i][0],
                                                     sources_flux[i][1], sources_flux_snr[i][1],
                                                     sources_flux[i][2], sources_flux_snr[i][2],
                                                     sources_radial_distance[i]))
                            detections_summary.write('\n')
                    detections_summary.close()
                    plt.close()

                else:
                    next_one = int(raw_input("Next one [1/0] [1]: ") or 1)
                    if next_one == 1:
                        if view_mode == 'single':
                            pass
                        elif view_mode == 'multiple':
                            continue
                    else:
                        return 0
                if not summary_file and (view_mode == 'multiple'):
                    is_close = int(raw_input("Close all windows? (1/0) [1]") or 1)
                    if is_close == 0:
                        return 0
                plt.close()
    except KeyboardInterrupt:
        return 0
    print("Failed files:", failed_files)

def run_flux_measure(detections=None, imagedir=None, savefile=None,
        resolution='0.3arcsec', suffix='combine.ms.auto.cont', target_wave=None,
        methods=['adaptive_aperture', 'gaussian', 'peak'], 
        bands=['B3','B4','B5','B6','B7','B8','B9'], ):
    """generate the catalogue base on the detections at different bands.

    It will read the source coordinates from flux measurement files and complete 
    their flux density at different bands.

    It also can be used to apply the flux correction by extrapolate the flux density 
    to a different wavelength.

    """
    detections = Table.read(detections, format='ascii')

    if resolution == '':
        res_string = ''
    else:
        res_string = resolution+'.'
    if savefile:
        with open(savefile, 'w+') as fp:
            fp.write('obj idx type ra dec radial_distance ')
            for band in bands:
                for m in methods:
                    fp.write("flux_{0}_{1} flux_snr_{0}_{1} ".format(band, m))
            fp.write('\n')
    for det in detections:
        obj = det['obj']
        ra, dec = det['ra'], det['dec'] # arcsec
        det_bands = {}
        for band in bands:
            obj_band_dir = os.path.join(imagedir, band, obj)
            image_name = "{}_{}_{}.{}image.fits".format(obj, band, suffix, res_string)
            image_fullpath = os.path.join(obj_band_dir, image_name)
            # print("flux measurement in:")
            # print(image_fullpath)
            if os.path.isfile(image_fullpath):
                source_flux, source_flux_snr = flux_measure(
                        image_fullpath, [[ra,dec],], target_wave=target_wave,
                        methods=methods, 
                        calculate_radial_distance=False)
                det_bands[band] = np.array([source_flux[0], source_flux_snr[0]])
            else:
                det_bands[band] = np.zeros((2, len(methods)))
        if savefile:
            with open(savefile, 'a+') as fp:
                fp.write('{obj} {idx} {type} {ra} {dec} {radial_distance} '.format(
                    obj=obj, idx=det['idx'], type=det['type'], ra=ra, dec=dec, \
                    radial_distance=det['radial_distance']))
                for band in bands:
                    for i in range(len(methods)):
                        fp.write("{} {} ".format(*det_bands[band][:,i]))
                fp.write('\n')

def run_gen_catalogue(fluxfiles=[], savefile=None,
        bands=['B3','B4','B5','B6','B7','B8','B9'], ):
    """generate the catalogue base on the flux measurements at different bands.
    """
    cat_base = Table.read(fluxfiles[0], format='ascii')
    cat_base['id'] = int(bands[0][1:])*10 + cat_base['id']
    cat_base_skcoords = SkyCoord(ra=cat_base['ra'], dec=cat_base['dec'])
    for ff,band in zip(fluxfiles[1:], bands[1:]):
        cat_new = Table.read(ff, format='ascii')
        cat_new['id'] = int(band[1:])*10 + cat_new['id']
        cat_new_skcoords = SkyCoord(cat_new['ra'], dec=cat_new['dec'])
        idx_new, idx_base, d2d, d3d = cat_base_skcoords.search_around_sky(cat_new_skcoords, 2*u.arcsec)
        idx_new_compset = []
        idx_base_compset = []
        cat_new[idx_new]['id'] = cat_base[idx_base]['id']
        matched_cat = join(cat_base[idx_base], cat_new[idx_new], key=['obj', 'id'], 
             table_names=[], uniq_col_name='{table_names}_{col_name}', join_type='inner')
        unique_cat = join(cat_base[idx_base_compset], cat_new[idx_new_compset], key=['obj', 'id'], 
             table_names=[], uniq_col_name='{table_names}_{col_name}', join_type='outer')
        cat_base = vstack(matched_cat, unique_cat)
    return cat_base

def combine_catalogues(fluxfiles=[], savefile=None,
        bands=['B3','B4','B5','B6','B7','B8','B9'], ):
    """generate the catalogue base on the flux measurements at different bands.
    """
    cat_base = Table.read(fluxfiles[0], format='ascii')
    band_base = bands[0]
    cat_base['idx'] = int(bands[0][1:])*10 + cat_base['idx']
    cat_base.rename_columns(cat_base.colnames[2:], [band_base+'_'+str(m) for m in cat_base.colnames[2:]])
    # add recomment ra and dec
    cat_base.add_column(cat_base[band_base+'_ra'], name='ra')
    cat_base.add_column(cat_base[band_base+'_dec'], name='dec')
    bands_done = [bands[0],]
    
    cat_base_skcoords = SkyCoord(ra=cat_base['ra'], dec=cat_base['dec'], unit='arcsec')
    for ff,band in zip(fluxfiles[1:], bands[1:]):
        cat_new = Table.read(ff, format='ascii')
        cat_new['idx'] = int(band[1:])*10 + cat_new['idx']
        cat_new.rename_columns(cat_new.colnames[2:], [band+'_'+str(m) for m in cat_new.colnames[2:]])
        cat_new_skcoords = SkyCoord(ra=cat_new[band+'_ra'], dec=cat_new[band+'_dec'], unit='arcsec')
        idx_new, idx_base, d2d, d3d = cat_base_skcoords.search_around_sky(cat_new_skcoords, 1*u.arcsec)
        if len(idx_new) > 1:
            cat_new['idx'][idx_new] = cat_base['idx'][idx_base]

            idx_new_compset = list(set(np.arange(0, len(cat_new))) - set(idx_new))
            idx_base_compset = list(set(np.arange(0, len(cat_base))) - set(idx_base))

            matched_cat = join(cat_base[idx_base], cat_new[idx_new], keys=['obj', 'idx'], 
                               join_type='inner')
            unique_cat = join(cat_base[idx_base_compset], cat_new[idx_new_compset], keys=['obj', 'idx'], 
                              join_type='outer')
            cat_base = vstack([matched_cat, unique_cat])
        else:
            #cat_base = vstack([cat_base, cat_new])
            cat_base = join(cat_base, cat_new, keys=['obj', 'idx'], 
                            join_type='outer')
        bands_done.append(band)
        #print(np.ma.masked_invalid([cat_base[b+'_ra'] for b in bands_done]))
        #print(">>>>>")
        #print(np.mean(np.ma.masked_invalid([cat_base[b+'_ra'] for b in bands_done]), axis=0))
        cat_base['ra'][:] = np.ma.mean([cat_base[b+'_ra'] for b in bands_done], axis=0)
        cat_base['dec'][:] = np.ma.mean([cat_base[b+'_dec'] for b in bands_done], axis=0)
        cat_base_skcoords = SkyCoord(ra=cat_base['ra'], dec=cat_base['dec'], unit='arcsec')
    return cat_base

def run_make_simulations(imagedir=None, objs=None, band=None, outdir='./', 
        resolution='0.3arcsec', repeat=1000, suffix='image.fits'):
    """make simulations for for the seleted calibrators and their bands

    imagedir: the image derectories generated by make_all_goodimages
    objs: should be the fields with detections
    bands: the wanted bands for each field
    outdir: the output directory, should be a place with enough free space
    resolution: the resolution of the image produced in make_all_goodimages, that is the base image
    repeat: how many simulations should be made for one field
    """
    if resolution == '':
        resolution_string = ''
    else:
        resolution_string = resolution+'.'
    for obj in objs:
        obj_imagedir = os.path.join(imagedir, band, obj)
        imagefile = os.path.join(obj_imagedir,
                '{}_{}_combine.ms.auto.cont.{}{}'.format(obj, band, resolution_string,
                                                         suffix))
        obj_outdir = os.path.join(outdir, obj)
        print(imagefile)
        if os.path.isfile(imagefile):
            if os.path.isfile(os.path.join(obj_outdir, '{}_simulation.txt'.format(obj))):
                print("Simulation file already exits, delete the {}_simulation.txt to overwrite.".format(obj))
                continue
            gen_sim_images(mode='image',imagefile=imagefile, outdir=obj_outdir, 
                           snr=(0.1,20), repeat=repeat)
            calculate_sim_images(obj_outdir, baseimage=imagefile, savefile=os.path.join(
                obj_outdir,'{}_simulation.txt'.format(obj)), plot=False, repeat=repeat, threshold=5, 
                second_check=False, snr_mode='peak', mode='image')

def run_calculate_effarea(imagedir=None, flux=np.linspace(0.01, 10, 500),  objs=None, band=None, 
        suffix='combine.ms.auto.cont', resolution='0.3arcsec', objs_nocenter=None, 
        almacal_catalogue=None,
        snr_threshold=5.0, savefile=None, **kwargs):
    """calculate the effecitve area for all the usable fields
    """
    obj_match = re.compile('^J\d*[+-]\d*$')
    
    if objs is None:
        objs = []
        for item in os.listdir(imagedir):
            if obj_match.match(obj):
                objs.append(item)
    if almacal_catalogue is not None:
        almacal_cat = Table.read(almacal_catalogue, format='ascii')
    effarea = Table()
    effarea['flux'] = flux
    effarea[band] = np.zeros_like(flux)
    for obj in objs:
        obj_dir = os.path.join(imagedir, obj)
        image = "{}_{}_{}.{}.image.fits".format(obj, band, suffix, resolution)
        image_fullpath = os.path.join(obj_dir, image)
        image_pbcorr = image.replace('image', 'pbcor.image')
        image_pbcorr_fullpath = os.path.join(obj_dir, image_pbcorr)
        if os.path.isfile(image_fullpath):
            mask_radius = 0.0
            if almacal_catalogue is not None:
                mask_code = almacal_cat[almacal_cat['obj'] == obj]['goodfield_{}'.format(band)].data[0]
                if mask_code < 0:
                    mask_radius = -1.0 * mask_code
            elif objs_nocenter is not None:
                mask_radius = 2.0
            _, objarea = calculate_effectivearea(flux, snr_threshold=snr_threshold, 
                        images=[image_fullpath,], images_pbcorr=[image_pbcorr_fullpath,],
                        central_mask_radius=mask_radius,**kwargs)
            effarea[band] += objarea
    if savefile:
        effarea.write(savefile, format='ascii')

def run_calculate_effarea2(imagedir=None, flux=np.linspace(0.01, 10, 500),  objs=None, band=None, 
        stats_dir=None, basedir=None, included_suffix='combine.ms.included.txt',
        suffix='combine.ms.auto.cont', resolution='0.3arcsec', objs_nocenter=None, 
        almacal_catalogue=None, 
        snr_threshold=5.0, savefile=None, debug=False, **kwargs):
    """calculate the effecitve area for all the usable fields

    This new version can also calculate the effective wavelength at different sensitivity
    
    imagedir : make_good_image/band
    stats_dir : gen_stats/band
    basedir : the ALMACAL base directory, where all the calibrated file is
    """
    obj_match = re.compile('^J\d*[+-]\d*$')
    
    if objs is None:
        objs = []
        for item in os.listdir(imagedir):
            if obj_match.match(obj):
                objs.append(item)
    if almacal_catalogue is not None:
        almacal_cat = Table.read(almacal_catalogue, format='ascii')
    flux_mean = 0.5*(flux[:-1] + flux[1:])
    effarea = Table()
    effarea['total_time'] = np.zeros_like(flux_mean)
    effarea['flux'] = flux_mean
    effarea[band] = np.zeros_like(flux_mean)
    effarea['effecitve_wavelength'] = np.zeros_like(flux_mean)

    for obj in objs:
        obj_basedir = os.path.join(basedir, obj)
        obj_dir = os.path.join(imagedir, obj)
        # find the selected files
        obj_band_selectfile = os.path.join(obj_dir, "{}_{}_{}".format(obj, band, included_suffix))
        band_obj_list = gen_filenames(listfile=obj_band_selectfile, basedir=obj_basedir)
        obj_savefile = os.path.join(stats_dir, obj+'.json')
        if os.path.isfile(obj_savefile):
            print("gen_stats file already exists, delete it for overwrite!")
            obj_spw_list = spw_stat(jsonfile=obj_savefile, bands=[band,])
        else:
            obj_spw_list = spw_stat(vis=band_obj_list, debug=debug,
                        bands=[band,], savefile=obj_savefile, 
                        **kwargs)
        obj_spwstat = discrete_spw(obj_spw_list, bands=[band,])

        # find the images
        image = "{}_{}_{}.{}.image.fits".format(obj, band, suffix, resolution)
        image_fullpath = os.path.join(obj_dir, image)
        image_pbcorr = image.replace('image', 'pbcor.image')
        image_pbcorr_fullpath = os.path.join(obj_dir, image_pbcorr)
        if os.path.isfile(image_fullpath):
            mask_radius = 0.0
            if almacal_catalogue is not None:
                mask_code = almacal_cat[almacal_cat['obj'] == obj]['goodfield_{}'.format(band)].data[0]
                if mask_code < 0:
                    mask_radius = -1.0 * mask_code
            elif objs_nocenter is not None:
                mask_radius = 2.0
            _, objarea = calculate_effectivearea(flux, snr_threshold=snr_threshold, 
                        images=[image_fullpath,], images_pbcorr=[image_pbcorr_fullpath,],
                        central_mask_radius=mask_radius,**kwargs)
            effarea[band] += objarea[:-1]
            objarea_diff = np.diff(objarea)
            effarea['effecitve_wavelength'][objarea_diff > 1e-8] += np.sum((obj_spwstat[band][0] * obj_spwstat[band][2]))
            effarea['total_time'][objarea_diff > 1e-8] += np.sum(obj_spwstat[band][2])
            
    effarea['effecitve_wavelength'] = effarea['effecitve_wavelength'] / effarea['total_time']
    if savefile:
        effarea.write(savefile, format='ascii')

def run_calculate_effarea(imagedir=None, flux=np.linspace(0.01, 10, 500),  objs=None, band=None, 
        suffix='combine.ms.auto.cont', resolution='0.3arcsec', objs_nocenter=None, 
        almacal_catalogue=None,
        snr_threshold=5.0, savefile=None, **kwargs):
    """calculate the effecitve area for all the usable fields
    """
    obj_match = re.compile('^J\d*[+-]\d*$')
    
    if objs is None:
        objs = []
        for item in os.listdir(imagedir):
            if obj_match.match(obj):
                objs.append(item)
    if almacal_catalogue is not None:
        almacal_cat = Table.read(almacal_catalogue, format='ascii')
    effarea = Table()
    effarea['flux'] = flux
    effarea[band] = np.zeros_like(flux)
    for obj in objs:
        obj_dir = os.path.join(imagedir, obj)
        image = "{}_{}_{}.{}.image.fits".format(obj, band, suffix, resolution)
        image_fullpath = os.path.join(obj_dir, image)
        image_pbcorr = image.replace('image', 'pbcor.image')
        image_pbcorr_fullpath = os.path.join(obj_dir, image_pbcorr)
        if os.path.isfile(image_fullpath):
            mask_radius = 0.0
            if almacal_catalogue is not None:
                mask_code = almacal_cat[almacal_cat['obj'] == obj]['goodfield_{}'.format(band)].data[0]
                if mask_code < 0:
                    mask_radius = -1.0 * mask_code
            elif objs_nocenter is not None:
                mask_radius = 2.0
            _, objarea = calculate_effectivearea(flux, snr_threshold=snr_threshold, 
                        images=[image_fullpath,], images_pbcorr=[image_pbcorr_fullpath,],
                        central_mask_radius=mask_radius,**kwargs)
            effarea[band] += objarea
    if savefile:
        effarea.write(savefile, format='ascii')

def run_mock_observation(radius=np.arange(1.0, 22.,1.5), imagedir=None, fNN=None, objs=None, band=None,
        suffix='combine.ms.auto.cont', resolution='0.3arcsec', objs_nocenter=None, 
        snr_threshold=5.0, savefile=None):
    """
    images: the real images accounting primary beam response
    fNN: callable, the cumulative function give the number of detections at a give sensitivity and area
    """
    # max_radius = np.max(radius)
    # step = radius[1] - radius[0]
    # central_mask_radius = radius[0]
    obj_match = re.compile('^J\d*[+-]\d*$')
    if objs is None:
        objs = []
        for item in os.listdir(imagedir):
            if obj_match.match(obj):
                objs.append(item)
    radius_mean = 0.5*(radius[:-1] + radius[1:])
    Nr_array = np.zeros_like(radius_mean)
    for obj in objs:
        obj_dir = os.path.join(imagedir, obj)
        image = "{}_{}_{}.{}.image.fits".format(obj, band, suffix, resolution)
        image_fullpath = os.path.join(obj_dir, image)
        image_pbcorr = image.replace('image', 'pbcor.image')
        image_pbcorr_fullpath = os.path.join(obj_dir, image_pbcorr)
        if os.path.isfile(image_fullpath):
            Nr = mock_observation(radius=radius, image=image_fullpath, 
                    image_pbcorr=image_pbcorr_fullpath, fNN=fNN, snr_threshold=snr_threshold,)
            Nr_array = Nr_array + np.array(Nr) 
    if savefile:
        with open(savefile, 'w+') as f:
            for i in range(len(radius_mean)):
                f.write("{} {}\n".format(radius_mean[i], Nr_array[i]))
    else:
        return radius_mean, Nr_array.tolist()

def run_number_counts_old(flist, detections_file=None, effective_area_file=None, band='B6',
        simulation_folder=None, ax=None, flux_mode=['aperture', 'gaussian'], completeness_mode='peak',
        default_simulation=None, objs_withdeafultsim=[], n_bootstrap=1000, savefile=None,
        plot=True, mode='cumulative'):
    """

    flist: flux list, [0.2, 0.6, 1.0]
    detections_file: the file with detection objs and their flux and flux snr, file format:
        obj B6 B6_SNR B7 B&_SRN
    band: the selected band
    effective_area_file: the file include two columns: flux and its effective area
    simulation_folder: the folder generated by run_make_simulations
    ax: the plot axes
    mode: can be 'cumulative' and 'differential'
    """
    flist = np.array(flist)
    cat = Table.read(detections_file, format='ascii')
    is_SMG = cat['type'] == 1
    is_Jet = cat['type'] == 2
    is_Unknown = cat['type'] == 0
    # tab = cat[is_SMG | is_Unknown]
    tab = cat[is_SMG]
    n_sources = len(tab)
    Ni_list = []
    Ni_bootstrap_array = np.zeros((n_sources, 2, n_bootstrap))

    # effective area
    effarea = np.loadtxt(effective_area_file, skiprows=1)
    cs_effarea = interpolate.interp1d(effarea[:, 0], effarea[:,1], fill_value="extrapolate")
    
    for i in range(n_sources):
        item = tab[i]
        if isinstance(flux_mode, (list, tuple)):
            flux1 = item['flux_'+flux_mode[0]]
            flux2 = item['flux_'+flux_mode[1]]
            # print(flux1, flux2)
            flux = 0.5*(flux1 + flux2)
            flux_err = 0.5*(item['flux_'+flux_mode[0]]/item['flux_snr_'+flux_mode[0]]
                            + item['flux_'+flux_mode[1]]/item['flux_snr_'+flux_mode[1]]) \
                            + np.diff([flux1, flux2])
        else:
            flux = item['flux_'+flux_mode]
            flux_err = 1.2*flux/item['flux_snr_'+flux_mode] 
        obj = item['obj']
        completeness_snr = item['flux_snr_'+completeness_mode]
        
        flux_bootstrap = flux + np.random.randn(n_bootstrap)*flux_err
        # define the completeness function
        if obj in objs_withdeafultsim:
            sim_jsonfile = default_simulation
        else:
            sim_jsonfile = os.path.join(simulation_folder, obj, obj+'_simulation.txt')
        if not os.path.isfile(sim_jsonfile):
            # raise ValueError('No simulation could be found for {}'.format(obj))
            print('Warning: using the default simulation results!')
            sim_jsonfile = default_simulation
        try:
            snr_list, apert_boost_list, comp_list, fake_rate_list = plot_sim_results(
                jsonfile=sim_jsonfile, snr=np.arange(0.2, 11, 0.2), plot=False)
        except:
            snr_list, apert_boost_list, comp_list, fake_rate_list = plot_sim_results(
                jsonfile=default_simulation, snr=np.arange(0.2, 11, 0.2), plot=False)
            print("Using default simulations: {}".format(default_simulation))

        cs_comp = scipy.interpolate.interp1d(snr_list, comp_list)
        def cs_comp2(snr):
            if snr<10.0:
                return cs_comp(snr)
            else:
                return 1.0
        # print('effarea:', cs_effarea(flux)/3600.)
        # print('completeness', cs_comp2(snr))
        Ni = 1 / (cs_effarea(flux)/3600.) / cs_comp2(completeness_snr)
        Ni_error = 1/(cs_effarea(flux)**2/3600)/cs_comp2(completeness_snr)*np.sqrt(cs_effarea(flux))
        Ni_bootstrap = 1/(cs_effarea(flux_bootstrap)/3600.)/cs_comp2(completeness_snr)
        #print(Si, Ni)
        Ni_list.append([flux, Ni, Ni_error, completeness_snr])
        Ni_bootstrap_array[i] = np.array([flux_bootstrap, Ni_bootstrap])
    Ni_array = np.array(Ni_list)
    
    if mode == 'differential':
        # calculation the deritive number counts
        dN = []
        dN_err = []
        dN_number = []
        for i in range(len(flist)-1):
            flist_select = (Ni_array[:,0] > flist[i]) & (Ni_array[:,0] < flist[i+1])
            dN.append(np.sum(Ni_array[:,1][flist_select]) / (0.5*(flist[i]+flist[i+1])))
            dN_number.append(np.sum(flist_select))
        dN_flist = 0.5 * (flist[1:] + flist[:-1])

        # print(dN_flist, dN, dN_number)
        if savefile:
            with open(savefile, 'w+') as sf:
                sf.write('flist dN_number dN\n')
                for i in range(len(dN_flist)):
                    sf.write("{} {} {}\n".format(dN_flist[i], dN_number[i], dN[i]))
        if plot:
            plt.loglog(dN_flist, dN)
        
    if mode == 'cumulative':
        # calculation the cumulative number counts
        NN = []
        NN_number = []
        NN_err = []
        for f in flist:
            NNi = np.sum(Ni_array[:,1][Ni_array[:,0] >= f])
            NNi_err2 = np.sum(Ni_array[:,2][Ni_array[:,0] >= f]**2)
            NN.append(NNi)
            NN_number.append(np.sum(Ni_array[:,0] >= f))
            NN_bootstrap = []
            for j in range(n_bootstrap):
                NN_bootstrap.append(np.sum(Ni_bootstrap_array[:,1,j][Ni_bootstrap_array[:,0,j] >= f]))
            NN_err.append(np.sqrt(np.std(NN_bootstrap)**2 + NNi_err2 + NNi))
        NN_err = (2*np.sqrt(np.array(NN_err)**2 + np.array(NN))).tolist()
        print(flist, NN_number, NN, NN_err)
        # print(flist, NN, np.sqrt(NN).tolist())
        NN_err_upper = []
        NN_err_lower = []
        poisson_error_table = Table.read(os.path.join(ROOT_DIR, 'code/data/Poisson_error.txt'), 
                format='ascii')
        for err, n in zip(NN_err, NN_number):
            if n<=30:
                poisson_select = poisson_error_table['n'] == n
                poisson_err_upper = poisson_error_table['upper'][poisson_select].data[0] - n
                poisson_err_lower = poisson_error_table['lower'][poisson_select].data[0] - n
            else:
                poisson_err_lower = np.sqrt(n)
                poisson_err_upper = np.sqrt(n)
            NN_err_upper.append(np.sqrt(err**2 + poisson_err_upper**2)) 
            NN_err_lower.append(np.sqrt(err**2 + poisson_err_lower**2))
        if savefile:
            print("flist", flist)
            print("NN_number", NN_number)
            print('NN', NN)
            print('NN_err_lower', NN_err_lower)
            print('NN_err_upper', NN_err_upper)
            with open(savefile, 'w+') as sf:
                sf.write('flist NN_number NN NN_err_lower NN_err_upper\n')
                for i in range(len(flist)):
                    sf.write("{} {} {} {} {}\n".format(flist[i], NN_number[i], NN[i], 
                                                       NN_err_lower[i], NN_err_upper[i]))
        if plot:
            plot_number_counts(flist=flist, NN=NN, NN_err=[NN_err_lower, NN_err_upper], ax=ax, band=band)

def run_number_counts(flist=None, detections_file=None, effective_area_file=None, band='B6',
        simulation_folder=None, ax=None, flux_mode='aperture', completeness_mode='peak',
        default_simulation=None, objs_withdeafultsim=[], n_bootstrap=1000, savefile=None,
        target_wavelength=None,
        plot=False, mode='cumulative', perturbation=0.1, n_points=None):
    """

    flist: flux list, [0.2, 0.6, 1.0]
    detections_file: the file with detection objs and their flux and flux snr, file format:
        obj B6 B6_SNR B7 B&_SRN
    band: the selected band
    effective_area_file: the file include two columns: flux and its effective area
    simulation_folder: the folder generated by run_make_simulations
    ax: the plot axes
    mode: can be 'cumulative' and 'differential'
    """
    cat = Table.read(detections_file, format='ascii')
    # cat.sort(['flux_'+flux_mode])
    n_sources = len(cat)
    cat_index = np.arange(n_sources)
    is_SMG = cat['type'] == 1
    is_Jet = cat['type'] == 2
    is_Unknown = cat['type'] == 0
    n_Unknown = int(np.sum(is_Unknown))
    n_SMG = int(np.sum(is_SMG))
    n_Jet = np.sum(is_Jet)
    SMG_frac = 1.0 * np.sum(is_SMG) / (np.sum(is_SMG) + np.sum(is_Jet))
    index_SMG = cat_index[is_SMG].tolist()
    index_Unknow = cat_index[is_Unknown].tolist()
    print("n_sources" ,n_sources)
    print("n_SMG" ,n_SMG)
    print("n_Unknown" ,n_Unknown)
    print("n_Jet" ,n_Jet)
    print("SMG fraction", SMG_frac)
    # tab = cat[is_SMG | is_Unknown]
    # tab = cat[is_SMG]

    if flist is None:
        flux_sorted = np.sort(cat['flux_'+flux_mode][index_SMG + index_Unknow])
        if not n_points:
            n_points = int(len(flux_sorted)/10.0)
        if mode == 'cumulative':
            idx_first = int(0.*len(flux_sorted))
            idx_last = int(0.95*len(flux_sorted))
            flist = np.logspace(np.log10(flux_sorted[idx_first]), 
                                np.log10(flux_sorted[idx_last]),n_points)
        if mode == 'differential':
            idx_first = int(0.0*len(flux_sorted))
            idx_last = int(0.95*len(flux_sorted))
            flist = np.logspace(np.log10(flux_sorted[idx_first]), np.log10(flux_sorted[idx_last]),
                                n_points)
    # effective area
    effarea = Table.read(effective_area_file, format='ascii')
    cs_effarea = interpolate.interp1d(effarea['flux'], effarea[band], fill_value="extrapolate")
    if 'effective_freq' in effarea.colnames:
        print("Correcting the effective frequency on the fly.")
        effective_wave = (const.c/(effarea['effective_freq'] * u.GHz)).to(u.um).value
        wavelength_scaling = flux_scaling(effective_wave, target_wavelength)
        func_scale = interpolate.interp1d(effarea['flux'], wavelength_scaling)
   
    # completeness function
    Ni_comp = np.zeros(n_sources)
    Si_deboosted = np.zeros(n_sources)
    for i in range(n_sources):
        item = cat[i]
        obj = item['obj']
        # define the completeness function
        if obj in objs_withdeafultsim:
            sim_jsonfile = default_simulation
        else:
            sim_jsonfile = os.path.join(simulation_folder, obj, obj+'_simulation.txt')
        if not os.path.isfile(sim_jsonfile):
            # raise ValueError('No simulation could be found for {}'.format(obj))
            print('Warning: using the default simulation results!')
            sim_jsonfile = default_simulation
        try:
            snr_list, apert_boost_list, comp_list, fake_rate_list = plot_sim_results(
                jsonfile=sim_jsonfile, snr=np.arange(0.2, 11, 0.2), plot=False)
        except:
            snr_list, apert_boost_list, comp_list, fake_rate_list = plot_sim_results(
                jsonfile=default_simulation, snr=np.arange(0.2, 11, 0.2), plot=False)
            print("Using default simulations: {}".format(default_simulation))
        cs_comp = interpolate.interp1d(snr_list, comp_list, fill_value='extrapolate')
        if flux_mode == 'aperture':
            flux_boosting = np.array(apert_boost_list[0])
        elif flux_mode == 'gaussian':
            flux_boosting = np.array(apert_boost_list[1])
        # print(apert_boost_list)
        # print(len(snr_list), len(flux_boosting[:,0]))
        func_deboosting = interpolate.interp1d(snr_list, flux_boosting[:,0], fill_value='extrapolate')
        def cs_comp2(snr):
            comp_return = cs_comp(snr)
            high_fedelity = (snr>10.0) #to avoid the truncating issue of the simulations
            comp_return[high_fedelity] = 1
            return comp_return
        Ni_comp[i] = cs_comp2(item['flux_snr_'+completeness_mode])
        Si_deboosted[i] = item['flux_'+flux_mode] / func_deboosting(
                          item['flux_snr_'+completeness_mode])

        if 'effective_freq' in effarea.colnames:
            Si_deboosted[i] = Si_deboosted[i] * func_scale(item['flux_'+flux_mode])

    # sim_jsonfile = os.path.join(simulation_folder, obj, obj+'_simulation.txt')
    ##
    #print("Using default simulations: {}".format(default_simulation))
    #sim_jsonfile = default_simulation
    #snr_list, apert_boost_list, comp_list, fake_rate_list = plot_sim_results(
    #    jsonfile=default_simulation, snr=np.arange(0.2, 11, 0.2), plot=False)
    #cs_comp = interpolate.interp1d(snr_list, comp_list, fill_value='extrapolate')
    #def cs_comp2(snr):
    #    comp_return = cs_comp(snr)
    #    high_fedelity = (snr>10.0)
    #    comp_return[high_fedelity] = 1
    #    return comp_return
    # print('effarea:', cs_effarea(flux)/3600.)
    # print('completeness', cs_comp2(snr))
    n_points = len(flist)
    N = np.zeros_like(flist)
    N_err = np.zeros_like(flist)
    N_number = np.zeros_like(flist)
    mode_simu = np.zeros((3, n_bootstrap, n_points))
    for i in range(n_bootstrap):
        N = N * 0.0 # reset the N
        N_err = N_err * 0.0 # reset the N
        N_number = N_number * 0.0 # reset N_number
        n_smgs_expect = SMG_frac*n_Unknown
        if n_Unknown > 0:
            # pearson possion error:
            index_unknown_bootstrap = np.random.choice(index_Unknow, int(n_smgs_expect), 
                                                       replace=False).tolist()
        else:
            index_unknown_bootstrap = []
        index_bootstrap = index_unknown_bootstrap + index_SMG 
        n_select = len(index_bootstrap)
        if perturbation:
            index_SMG_bootstrap = np.random.choice(index_bootstrap, int((1-perturbation)*n_select), 
                                                   replace=False).tolist()
        else:
            index_SMG_bootstrap = index_bootstrap
        # print("two index:", index_unknown_bootstrap, index_SMG_bootstrap)
        cat_bootstrap = cat[index_SMG_bootstrap]
        Ni_comp_boostrap = Ni_comp[index_SMG_bootstrap]
        Si_deboosted_boostrap = Si_deboosted[index_SMG_bootstrap]
        # print('n select: unknown:{}; SMG:{}'.format(len(index_unknown_bootstrap), len(index_SMG_bootstrap)))
        # flux = cat_bootstrap['flux_'+flux_mode]
        flux = Si_deboosted_boostrap
        completeness_snr = cat_bootstrap['flux_snr_'+completeness_mode]
        flux_err = flux/cat_bootstrap['flux_snr_'+flux_mode]
        flux_bootstrap = flux + np.random.randn(len(index_SMG_bootstrap)) * 1.1 * flux_err # add additional 10% error
        Ni = 1 / (cs_effarea(flux_bootstrap)/3600.) / Ni_comp_boostrap#cs_comp2(completeness_snr)
        Ni_err = 1/(cs_effarea(flux_bootstrap)**2/3600)*np.sqrt(cs_effarea(flux_bootstrap))/Ni_comp_boostrap #cs_comp2(completeness_snr)
        
        # testing code
        # testing the each Ni, please change the n_bootstrap to a small number
        # for itt in zip(flux_bootstrap, flux_err, completeness_snr, Ni):
            # if itt[0] > 2.0:
                # print(itt)
        # print("Ni error:", Ni_err/Ni)

        if mode == 'cumulative':
            for j in range(n_points):
                fi_select = (flux >= flist[j])
                N[j] = np.sum(Ni[fi_select])
                N_err[j] = np.sqrt(np.sum(Ni_err[fi_select]**2))
                # N_err[j] = np.std(Ni[fi_select])
                N_number[j] = np.sum(fi_select)
        if mode == 'differential':
            for j in range(n_points-1):
                fi_select = (flux > flist[j]) & (flux <= flist[j+1])
                N[j] = np.sum(Ni[fi_select]) / (flist[j+1] - flist[j])
                N_err[j] = np.sqrt(np.sum(Ni_err[fi_select]**2)) / (flist[j+1] - flist[j])
                N_number[j] = np.sum(fi_select)
        mode_simu[0,i,:] = N
        mode_simu[1,i,:] = N_err
        mode_simu[2,i,:] = N_number
        # print("N err:", N_err/N)
        # print("N=", N)
        # print('N_number=', N_number)

    # print(mode_simu[0])
    if mode == 'differential':
        flist = 0.5*(flist[:-1] + flist[1:])
        mode_simu = mode_simu[:,:,:-1] # remove the last column
    mode_NN_number = np.mean(mode_simu[2], axis=0)
    mode_NN_number_err = np.std(mode_simu[2], axis=0)
    mode_NN = np.mean(mode_simu[0], axis=0)
    # if 
    mode_NN_err = np.std(mode_simu[0], axis=0)
    #mode_NN_err = np.sqrt(np.sum(mode_simu[1]**2, axis=0))/mode_NN_number # another way to calculate the err
    # print('mode_NN', mode_NN)
    # print('mode_NN_err', mode_NN_err)
    # print('mode_NN_number', mode_NN_number)
                
    if True: # add Poisson errors
        NN_err_upper = np.zeros_like(mode_NN_err)
        NN_err_lower = np.zeros_like(mode_NN_err)
        poisson_error_table = Table.read(os.path.join(ROOT_DIR, 'code/data/Poisson_error.txt'), 
                format='ascii')

        # add Poisson Errors
        for i in range(len(flist)):
            n_float, n_float_err = mode_NN_number[i], mode_NN_number_err[i]
            n = int(mode_NN_number[i])
            if n<=30: # read the Poisson error table
                poisson_select = (poisson_error_table['n'] == n)
                poisson_err_upper = (poisson_error_table['upper'][poisson_select].data[0]-n)/n*mode_NN[i]
                poisson_err_lower = (poisson_error_table['lower'][poisson_select].data[0]-n)/n*mode_NN[i]
            else: # direct use sqrt(n) as a proxmimation
                poisson_err_lower = np.sqrt(n_float+0.5*n_float_err/np.sqrt(n_float))/n*mode_NN[i]
                poisson_err_upper = np.sqrt(n_float+0.5*n_float_err/np.sqrt(n_float))/n*mode_NN[i]
            NN_err_upper[i] = (np.sqrt(mode_NN_err[i]**2 + poisson_err_upper**2)) 
            NN_err_lower[i] = (np.sqrt(mode_NN_err[i]**2 + poisson_err_lower**2))
    if savefile:
        print("flist", flist)
        print("NN_number", mode_NN_number)
        print('NN', mode_NN)
        print('NN_err_lower', NN_err_lower)
        print('NN_err_upper', NN_err_upper)
        with open(savefile, 'w+') as sf:
            if mode=='cumulative':
                sf.write('flist N_number N Nerr_lower Nerr_upper\n')
            if mode=='differential':
                sf.write('flist dN_number dN dNerr_lower dNerr_upper\n')
            for i in range(len(flist)):
                sf.write("{} {} {} {} {}\n".format(flist[i], mode_NN_number[i], mode_NN[i], 
                                               NN_err_lower[i], NN_err_upper[i]))

    if plot:
        plot_number_counts(flist=flist, NN=mode_NN, NN_err=[NN_err_lower, NN_err_upper], ax=ax, band=band)

def plot_number_counts(datafile=None, flist=None, NN=None, NN_err=None, ax=None, band=None,
        model_dir=None, show_models=True):
    if ax is None:
        fig = plt.figure(figsize=(6, 8))
        ax = fig.add_subplot(111)
        ax.set_xscale("log")
        ax.set_yscale("log")
    if datafile:
        if isinstance(datafile, str):
            tab = Table.read(datafile, format='ascii')
            flist = tab['flist']
            NN = tab['NN']
            NN_err_upper = tab['NN_err_upper']
            NN_err_lower = tab['NN_err_lower']
            NN_err = [NN_err_lower, NN_err_upper]
        elif isinstance(datafile, (list, tuple)):
            for dfile, dband in zip(datafile, band):
                plot_number_counts(dfile, ax=ax, band=dband, show_models=False,model_dir=model_dir)
    if True:
        ax.plot(flist, NN, 'ko',label='ALMACAL')
        ax.plot(flist, NN, 'k--')
        ax.errorbar(flist, NN, yerr=NN_err, fmt='k')
        ax.set_xlabel('Flux [mJy]')
        ax.set_ylabel('N [deg$^{-2}$]')
        ax.set_title(band)

        if band == 'B7':
            flist_oteo2016 = [0.4, 1.0]
            NN_oteo2016 = [17e3, 3.8e3]
            NN_oteo2016_error_upper = [14.3e3, 3.2e3]
            NN_oteo2016_error_lower = [14.0e3, 3.1e3]
            ax.plot(flist_oteo2016, NN_oteo2016, 'yo', label='Oteo et al (2016) B7')
            ax.errorbar(flist_oteo2016, NN_oteo2016, yerr=NN_oteo2016_error_upper, fmt='y')

        if band == 'B6':
            flist_oteo2016 = [0.2, 0.8]
            NN_oteo2016 = [5.6e3, 0.9e3]
            NN_oteo2016_error_upper = [4.7e3, 0.8e3]
            NN_oteo2016_error_lower = [4.6e3, 0.7e3]
            ax.plot(flist_oteo2016, NN_oteo2016, 'yo', label='Oteo et al. (2016) B6')
            ax.errorbar(flist_oteo2016, NN_oteo2016, yerr=NN_oteo2016_error_upper, fmt='y')

        if band == 'B8':
            flist_anne2020 = [0.67, 2.50]
            NN_anne2020_log = [4.8, 4.1]
            NN_anne2020_log_error = [4.7, 3.9]
            NN_anne2020 = np.array(10)**np.array(NN_anne2020_log)
            NN_anne2020_error = np.array(NN_anne2020_log_error)*np.array(NN_anne2020)*np.log(10)
            ax.plot(flist_anne2020, NN_anne2020, 'yo', label='Klitsch et al. (2020)')
            ax.errorbar(flist_anne2020, NN_anne2020, yerr=NN_anne2020_error, fmt='y')
            # Plot the model of Lagos
            # flist_Lagos_B8 = 10**np.linspace(-1.0, 1.0, 9)
            # NN_lagos_B8 = np.array([102458.04, 70664.18, 44428.86, 24411.87, 11827.86, 5120.97, 1919.28, 602.10, 164.67])
            # plt.plot(flist_Lagos_B8, NN_lagos_B8)

        if band == 'B6':
            # show ASPECS
            flist_aspecs_data  = np.array([[31.6, 39.8], [39.8, 50.1], [50.1, 63.1],[63.1, 79.4], [79.4,100.], 
                             [100., 125.9], [125.9,158.5], [158.5,199.5], [199.5, 251.2], 
                             [251.2,316.2], [316.2, 398.1], [398.1,501.2], [501.2,631.0],
                             [631.0,794.3], [794.3,1000.], [1000.,1258.9], [1258.9,1584.9]])
            flist_aspecs = np.mean(flist_aspecs_data, axis=1) * 1e-3
            NN_aspecs = [47400, 35100, 30700,26500,26500,19300,14800,10600,8700,7800,5300,4300,1720,1720,
                         860,860,1600]
            NN_aspecs_error_lower = [8200,6200,5500,4900,5100,4400,3500,2900,2400,2300,1900,1600,840,840,
                                     490,490,-1]
            NN_aspecs_error_upper = [8900, 6900,6200,5600,5600,4700,4000,3400,3000,2700,2400,2100,1250,
                                     1250,820,820,-1]
            ax.plot(flist_aspecs, NN_aspecs, 'bs', label='ASPECS-LP 1.2mm')
            ax.errorbar(flist_aspecs, NN_aspecs, yerr=NN_aspecs_error_upper, fmt='b', linestyle='..')

            # show Fujimoto et al. 2016
            flist_fujimoto = [0.002, 0.015, 0.027, 0.047, 0.084, 0.150, 0.267, 0.474, 0.843]
            NN_fujimoto_log = [6.6, 5.5, 5.3, 5.0, 4.8, 4.6, 4.2, 3.7, 2.8]
            NN_fujimoto_log_error_upper = [0.5, 0.3, 0.2, 0.2, 0.1, 0.1, 0.2, 0.2, 0.6]
            NN_fujimoto_log_error_lower = [1.1, 0.4, 0.3, 0.2, 0.2, 0.1, 0.2, 0.3, 0.5]
            NN_fujimoto = np.array(10)**NN_fujimoto_log
            NN_fujimoto_error_upper = np.array(NN_fujimoto_log_error_upper)*np.array(NN_fujimoto)*np.log(10)
            NN_fujimoto_error_lower = np.array(10)**NN_fujimoto_log_error_lower
            ax.plot(flist_fujimoto[1:], NN_fujimoto[1:], 'm>', label='Fujimoto et al. (2016) 1.2mm')
            ax.errorbar(flist_fujimoto[1:], NN_fujimoto[1:], yerr=NN_fujimoto_error_upper[1:], fmt='m', linestyle='-.')

            # # show Stach 2018
            # flist_stach = [4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5]
            # NN_stach = [385.3, 216.7, 106.2, 53.6, 29.6, 20.0, 10.5, 5.3, 2.1,2.1,1.0]
            # NN_stach_error_upper = [21.1, 17.3, 11.4, 8.4, 6.5, 5.7, 4.4,3.5, 2.8, 2.8, 2.4]
            # NN_stach_error_lower = [7.7, 6.6, 3.5, 2.5, 1.9, 1.8, 1.2, 0.9, 0.6, 0.6, 0.5]
            # ax.plot(flist_stach, NN_stach, 'co', label='Stach et al. (2018)')
            # ax.errorbar(flist_stach, NN_stach, yerr=NN_stach_error_upper, fmt='c')
        if False:
            # show the theoretical model
            Lagos_model_dir = os.path.join(model_dir, 'ModelsLagos')
            Lagos_B6 = np.loadtxt(os.path.join(Lagos_model_dir,'Band6.txt'))
            Lagos_B7 = np.loadtxt(os.path.join(Lagos_model_dir,'Band7.txt'))
            Lagos_B8 = np.loadtxt(os.path.join(Lagos_model_dir,'Band8.txt'))
            Lagos_B9 = np.loadtxt(os.path.join(Lagos_model_dir,'Band9.txt'))
            ax.plot(10**Lagos_B6[:,0], Lagos_B6[:,1], 'r--', label='Lagos et al. (2019)')
            ax.plot(10**Lagos_B7[:,0], Lagos_B7[:,1], 'r--', )
            ax.plot(10**Lagos_B8[:,0], Lagos_B8[:,1], 'r--', )
            ax.plot(10**Lagos_B9[:,0], Lagos_B9[:,1], 'r--', )


        ax.legend()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#   The pipeeline for number counts       #
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''
run_gen_all_obstime('/scratch/ALMACAL', outdir='gen_all_obstime', bands=['B3','B4','B5','B6','B7','B8'], debug=False)
run_gen_all_images('/jchen/work/ALMACAL', outdir='gen_all_images', bands=['B8',])
run_auto_classify_goodimags(imagedir='gen_all_images', outdir='make_classifications', bands=['B8',])
B8_list = np.loadtxt('B8_ordered.txt', dtype=str)
run_manual_inspection(classifiedfolder='make_classifications', objlist=B8_list, bands=['B8',])
run_make_all_goodimags(classifiedfolder='make_classifications', objlist=B8_list, 
                       basedir='/jchen/work/ALMACAL', outdir='make_all_good_images')
run_check_SMGs('make_all_good_images', objs=None, bands=['B3','B4','B5','B6','B7','B8'], ncol=2,
                summary_file='check_SMGs.txt', outdir='./check_SMGs')
run_gen_stats('/scratch/ALMACAL', imagedir='make_all_good_images', outdir='gen_stats', 
               bands=['B3','B4','B5','B6','B7','B8'], debug=False)
run_measure_flux('make_all_good_images', objs=None, bands=['B6','B7','B8'], focus_band='B8', 
                 selected_resolution='0.3arcsec', summary_file='B8_SMGs.txt', target_wave=650)
run_calculate_effarea(imagedir='make_good_images', flux=np.linspace(0.1, 10, 50), objs=B8_list, 
                      band='B8', savefile='B8_effarea.txt')
run_make_simulations(imagedir='make_good_images', objs=B8_dets_uniq, band='B8', outdir='simulations/B8')
run_number_counts(flist, detections_file='B8_SMGs.txt', effective_area_file='B8_effarea.txt', band='B8'
                  default_simulation='', simulation_folder='./simulations/B8')

run_number_counts(flist=np.logspace(-1.1, 0.4, 15),
  detections_file='B6/B6_SMGs_run3.0.3arcsec.flux.txt',
  effective_area_file='B6/B6_effarea2.txt',
  band='B6',
  default_simulation='../main_run2/simulations/B8/J0120-2701/J0120-2701_simulation.txt',
  simulation_folder='simulations/B8',
  flux_mode='aperture',
  plot=True,
  savefile='B6/B6_number_counts.txt')

# flist used for different bands:
    B5: flist=np.array([0.3,0.7])

'''
