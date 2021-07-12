# This file include all the tasks for ALMACAL number counts projects

# The task list include:
# 1. gen_obstime: generate the observational time for the whole dataset
# 2. gen_dirty_image: generate the dirty image for all the observations of one calibrator
# 3. show_image: interative way to inspect the images manually 
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
from astropy.stats import sigma_clipped_stats, sigma_clip

import analysisUtils as au
# from analysisUtils import tbtool
# from imaging_utils import make_cont_img

tb = au.tbtool()

def flatten(t):
    flat_list = [item for sublist in t for item in sublist]
    return flat_list

def gen_filenames(listfile=None, dirname=None, basedir=None, debug=False, exclude_aca=True):
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
            line = l.strip()
            if basedir:
                line = os.path.join(basedir, line)
            file_list.append(line)
    return file_list

def savelist(l, filename=None, outdir='./', file_mode='a+'):
    with open(os.path.join(outdir, filename), file_mode) as f:
        for item in l:
            f.write(item+'\n')

def gen_image(vis=None, band=None, outdir='./', niter=0, exclude_aca=False, check_calibrator=False, debug=False, update_raw=False, overwrite=False, **kwargs):
    """make images for one calibrator on all or specific band

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
    if band is not None:
        band_match = re.compile('_(?P<band>B\d{1,2})$')
        if band_match.search(vis):
            obs_band = band_match.search(vis).groupdict()['band']
            if obs_band != band:
                print("Invalid band!")
                return False

    basename = os.path.basename(vis)
    myimagename = os.path.join(outdir, basename + '.cont.auto')
    if os.path.isfile(myimagename+'.fits') and not overwrite:
        print('image already existing, return without making new.')
        return 0

    if exclude_aca:
        try:
            tb.open(vis + '/ANTENNA')
            antenna_diameter = np.mean(tb.getcol('DISH_DIAMETER'))
            tb.close()
        except:
            return False
        if antenna_diameter < 12.0:
            if debug:
                print("Excuding data from {}".format(antenna_diameter))
            return False
    try:
        make_cont_img(vis=vis, clean=True, myimagename=myimagename, outdir=outdir, niter=niter, **kwargs)
    except:
        print("Error in imaging {}".format(vis))
    if check_calibrator:
        imstat_info = imstat(myimagename+'.image')
        if imstat_info['flux'][0] > 0.05: #flux/flux density large than 0.05 Jy
            print("Maybe point source subtraction is failed")
            outfile = os.path.join(outdir, os.path.basename(vis)+".point.cl")
            uvmodelfit(vis=vis, niter=5, comptype="P", sourcepar=[1.0, 0.0, 0.0], varypar=[True,False,False], outfile=outfile)
            # uvmodelfit(vis=vis, niter=5, comptype="G", outfile=outfile,
                       # sourcepar=[1.0, 0.0, 0.0, 1.0, 0.9, 0], 
                       # varypar=[True, False, False, True, True, True])
            ft(vis=vis, complist=outfile)
            uvsub(vis=vis,reverse=False)
            outputvis = os.path.join(outdir, 'updated_data', os.path.basename(vis))
            os.system('mkdir -p {}'.format(os.path.dirname(outputvis)))
            split(vis=vis, outputvis=outputvis)
            if update_raw:
                print('removing', vis)
                rmtables(vis)
                print('copying', outputvis, vis)
                os.system('mv {} {}'.format(outputvis, vis))
                outputvis = vis
            rmtables(myimagename+'.*')
            make_cont_img(vis=outputvis, clean=True, myimagename=myimagename, outdir=outdir, niter=1000, only_fits=True, **kwargs)
            return 0
            
    exportfits(imagename=myimagename+'.image', fitsimage=myimagename+'.fits')
    rmtables(tablenames=myimagename+'.*')
    return 0

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
    for infile in filelist:
        if bands is not None:
            for band in bands:
                outdir_band = os.path.join(outdir, band)
                os.system('mkdir -p {}'.format(outdir_band))
                gen_image(infile, band=band, outdir=outdir_band, exclude_aca=exclude_aca, **kwargs)
        else:
            gen_image(infile, outdir=outdir, exclude_aca=exclude_aca, **kwargs)

def show_images(fileglob=None, filelist=None, basedir=None, mode='auto', nrow=3, ncol=3, savefile=None, debug=False):
    """show images in an interative ways, and record the input from inspector
    
    Parameters:
        fileglob: the match patter of the image filenames, like: './*.image'
    """
    
    if fileglob:
        if debug:
            print(fileglob)
        all_files = glob.glob(fileglob)
    elif filelist:
        if debug:
            print(filelist)
        if not os.path.isfile(filelist):
            raise ValueError("{} is not found".format(filelist))
        all_files = []
        if basedir is None:
            raise ValueError("basedir should be defined along with filelist")
        if isinstance(filelist, str):
            try:
                flist = []
                with open(filelist) as f:
                    filelist_lines = f.readlines()
                for line in filelist_lines:
                    flist.append(line.strip()+'.cont.auto.fits')
            except:
                raise ValueError("file {} cannot be open".format(filelist))
        elif isinstance(filelist, list):
            flist = filelist
        else:
            raise ValueError("Wrong file type of filelist!")
        for item in flist:
            all_files.append(os.path.join(basedir, item))

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

def check_image(img, plot=False, radius=8, debug=False, sigmaclip=True, check_flux=True, 
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
        print('mean: {}, std:{}'.format(mean, sigma))
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
    # check the fluxval
    if check_flux:
        # check the flux value of the subtracted point source
        if uidname:
            fluxval = search_flux(uidname)
            if debug:
                print(uidname)
                print("fluxval: {}".format(fluxval))
            if fluxval > 0:
                if fluxval/minimal_fluxval < 1:
                    if debug:
                        print("Rejected, wrong fluxval!")
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
                return
    if isinstance(imgs, str):
        all_files = glob.glob(imgs)
    elif isinstance(imgs, list):
        all_files = imgs

    p_uidimg = re.compile('(?P<uidname>uid___\w+.ms.split.cal.J\d+[+-]\d+_B\d+).cont.auto.fits')

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
        # 
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

def check_images_manual(imagedir=None, goodfile=None, badfile=None, debug=False, ncol=1, nrow=3):
    '''visual inspection of the classification results from check_images
    
    '''
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

def make_good_image(vis=None, basename='', basedir=None, outdir='./', concatvis=None, debug=False, 
                    only_fits=True, niter=1000, clean=True, pblimit=-0.01, fov_scale=2.0, 
                    computwt=True, save_psf=True,
                    uvtaper_list=[['0.3arcsec'], ['0.6arcsec']], 
                    uvtaper_scale=None,#[1.5, 2.0], 
                    **kwargs):
    """make the final good image with all the good observations

    uvtaper_list: [['0.3arcsec'], ['0.8arcsec']], outer taper width
    uvtaper_scale: [1.0, 1.7], 1.5 and 2 times and three times worse
    """
    if len(vis) < 1:
        return False
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
    if computwt:
        vis_new = []
        for v in vis:
            v_new = os.path.join(outdir, os.path.basename(v))
            if os.path.isdir(v_new):
                print('Found splitted file: {}'.format(v_new))
                continue
            # os.system('rm -rf {}*'.format(v_new))
            print('Splitting... {} {}'.format(v, v_new))
            #os.system('cp -r {} {}'.format(v, v_new))
            if not split(vis=v, outputvis=v_new, datacolumn='corrected'):
                print("No corrected data existing, spliting the data column")
                if not split(vis=v, outputvis=v_new, datacolumn='data'):
                    print("Data may be corrupted: {}".format(v))
                    os.system('rm -rf {}'.format(v_new))
                    continue

            statwt(v_new, datacolumn='data')
            vis_new.append(v_new)
        print(vis_new)
        vis = vis_new
    if os.path.isdir(concatvis):
        print("Skip concating, file exists....")
    else:
        concat(vis=vis, concatvis=concatvis)
    if uvtaper_list:
        for uvtaper in uvtaper_list:
            if uvtaper == []:
                uvtaper_name = ''
            else:
                uvtaper_name = '.'+uvtaper[0]
            myimagename = concatvis+'.auto.cont{}'.format(uvtaper_name)
            make_cont_img(vis=concatvis, myimagename=myimagename, clean=clean, niter=niter, 
                          pblimit=pblimit, fov_scale=fov_scale, myuvtaper=uvtaper, debug=debug, 
                          **kwargs)
    if uvtaper_scale:
            myimagename = concatvis+'.auto.cont'
            make_cont_img(vis=concatvis, myimagename=myimagename, clean=clean, niter=niter, 
                          pblimit=pblimit, fov_scale=fov_scale, uvtaper_scale=uvtaper_scale, 
                          debug=debug, **kwargs)
    if only_fits:
        for i in glob.glob(concatvis+'*.image'):
            exportfits(imagename=i, fitsimage=i+'.fits')
        if save_psf:
            for psf in glob.glob(concatvis+'*.psf'):
                exportfits(imagename=psf, fitsimage=psf+'.fits')

        rmtables(concatvis+'*')
        rmtables(os.path.join(outdir,'uid___*'))
        os.system('rm -rf {}'.format(os.path.join(outdir,'uid___*.flagversions')))

def calculate_completeness(objfolder, vis=None, baseimage=None, n=20, repeat=10, snr=np.arange(1,20,0.5), 
        suffix='.cont.auto', known_file=None, obj=None, band=None, basename=None, savefile=None, 
        threshold=5.0, plot=False, snr_mode='integrated', **kwargs):
    """simulation the completeness of source finding algorithm

    mode:
        peak: snr is the peak value
        integrated: snr is the integrated value
    """
    # one time function
    f_mean = lambda x: np.mean(x)
    
    # image statistics
    im_head = imhead(baseimage)
    im_beam = im_head['restoringbeam']
    im_incr = im_head['incr']
    im_info = imstat(baseimage)
    # beamsize = np.pi*a*b/(4*np.log(2))
    beamsize = np.pi/(4*np.log(2))* im_beam['major']['value'] * im_beam['minor']['value'] / (im_incr[0]/np.pi*180*3600)**2
    rms = im_info['rms']
    rms_flux = rms * 1000 # convert to mJy/pixel
    # sensitivity = 1000 * calculate_sensitivity(vis) # convert into mJy
    print('rms',rms)
    print('beamsize',beamsize) 
    print('rms_flux',rms_flux)
    
    flux_match = re.compile('(?P<obj>J\d*[+-]\d*)_(?P<band>B\d+)\w*.snr(?P<snr>\d+.\d+).run(?P<run>\d+)')
    if basename is None:
        basename = os.path.join(objfolder, '{obj}_{band}_combine.ms'.format(obj=obj, band=band))
        print('basename', basename)
    all_fake_images = glob.glob(basename+'*{}.fits'.format(suffix))

    flux_input_theory = []
    flux_input_list = []
    flux_input_autolist = []
    flux_input_foundlist = []
    flux_found_autolist = []
    snr_input_list = []
    # snr_inputfound_comp = []
    detection_input_array = np.array([[0., 0.]])
    detection_found_array = np.array([[0., 0.]])
    for s in snr:
        print('SNR: {}'.format(s))
        flux_input_theory.append(snr*s)
        n_input = 0
        n_found = 0
        # flux = s*sensitivity
        # flux_input_list.append(flux)
        for run in np.arange(repeat):
            #print('SNR: {}, RUN:{}'.format(s, run))
            #print('flux:', flux)
            img = "{basename}.snr{snr}.run{run}{suffix}.fits".format(basename=basename, snr=s, run=run, suffix=suffix)
            sf_return = source_finder(img, sources_file=basename+'.snr{}.run{}.txt'.format(s, run), known_file=known_file, **kwargs)
            # if sf_return == 0:
                # continue
            flux_input, flux_input_auto, flux_found_auto, idxs = sf_return 
            # print('flux_input', flux_input)
            # print('flux_input_auto', flux_input_auto)
            # print('flux_found_auto', flux_found_auto)
            # print('idxs', idxs)
            
            if len(idxs[0])<1:
                print("Skip snr={}, run{}".format(s, run))
                continue
            flux_input_list.append(flux_input)
            flux_input_autolist.append(flux_input_auto)
            flux_found_autolist.append(flux_found_auto)

            #calculate the snr
            if mode == 'peak':
                snr_input = flux_input_auto[:,2] / rms_flux
                snr_input_list.append(snr_input)
            elif mode == 'integrated':
                snr_input = flux_input / rms_flux
                snr_input_list.append(snr_input)
            # the snr of the detection
            snr_input_found = flux_input_auto[:,2][idxs[0]] / rms_flux
            snr_input_failed = flux_input_auto[:,2][idxs[2]] / rms_flux
            snr_found_input = flux_found_auto[:,2][idxs[1]] / rms_flux
            snr_found_fake = flux_found_auto[:,2][idxs[3]] / rms_flux
            # snr_input_foundlist.append(snr_input_found)
            # print('snr_inputfound',snr_inputfound)
            # print('snr_inputfound_faild',snr_inputfound_failed)
            # print('detection_array', detection_array)
            if len(snr_input_found) > 0:
                detection_input_array = np.vstack([detection_input_array, np.array(zip(snr_input_found, 
                                                                      np.ones_like(snr_input_found)))])
            if len(snr_input_failed) > 0:
                detection_input_array = np.vstack([detection_input_array, np.array(zip(snr_input_failed, 
                                                                      np.zeros_like(snr_input_failed)))])
            if len(snr_found_input) > 0:
                detection_found_array = np.vstack([detection_found_array, np.array(zip(snr_found_input, 
                                                                      np.ones_like(snr_found_input)))])
            if len(snr_found_fake) > 0:
                detection_found_array = np.vstack([detection_found_array, np.array(zip(snr_found_fake, 
                                                                      np.zeros_like(snr_found_fake)))])

    # save data into json
    data_saved = {}
    snr_flat = [item for sublist in snr_input_list for item in sublist]
    flux_input_flat = [item for sublist in flux_input_list for item in sublist]
    flux_input_aperture_flat = [item for sublist in flux_input_autolist for item in sublist[:,0]]
    flux_input_gaussian_flat = [item for sublist in flux_input_autolist for item in sublist[:,1]]
    flux_aperture_flat = [item for sublist in flux_found_autolist for item in sublist[:,0]]
    flux_gaussian_flat = [item for sublist in flux_found_autolist for item in sublist[:,1]]
    data_saved['snr_input'] = snr_flat
    data_saved['flux_input'] = flux_input_flat
    data_saved['flux_input_aperture'] = flux_input_aperture_flat
    data_saved['flux_input_gaussian'] = flux_input_gaussian_flat
    data_saved['flux_gaussian'] = flux_gaussian_flat
    data_saved['flux_aperture'] = flux_aperture_flat
    data_saved['detection_snr'] = detection_input_array[:,0].tolist()
    data_saved['detection_input_array'] = detection_input_array.tolist()
    data_saved['detection_found_array'] = detection_found_array.tolist()
    if savefile:
        with open(savefile, 'w') as fp:
           json.dump(data_saved, fp)
    if plot:
        plot_completeness(data_saved)
    # return data_saved
    #flux_list, flux_peak_list, flux_found_list, completeness_list

def calculate_effectivearea(flux=np.linspace(0.1, 1, 10), snr_threshold=5.0, images=None, 
        images_pbcorr=None, fovscale=2.0):
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
        pixel_scale = header['CDELT1']*3600
        pix2area = pixel_scale**2  # pixel to arcsec^2
        freq = header['CRVAL3']
        lam = (const.c/(freq*u.Hz)).decompose().to(u.um)
        fov = 1.02 * (lam / (12*u.m)).decompose()* 206264.806
        ny, nx = data.shape[-2:]
        data_masked = np.ma.masked_invalid(data.reshape(ny, nx))
        mean, median, std = sigma_clipped_stats(data_masked, sigma=5.0, iters=5)  
        mask = data_masked.mask & (data_masked.data >= 5.0*std)
        pbcor = (data / data_pbcor).reshape(ny, nx)
        pbcor[mask] = 0.0
        if fovscale:
            rlimit = 0.5 * fovscale * fov.value
            # print('rlimit', rlimit)
            x = (np.arange(0, nx) - nx/2.0) * pixel_scale
            y = (np.arange(0, ny) - ny/2.0) * pixel_scale
            r = np.sqrt(x**2 + y**2)
            # print('max r:', np.max(r))
            pbcor[r > rlimit] = 0.0
        effarea_list = []
        for f in flux:
            snr = f / (std * 1000) # from Jy to mJy
            snr_map = snr * pbcor
            area = np.sum(snr_map > snr_threshold) * pix2area
            effarea_list.append(area)
        effarea = effarea + np.array(effarea_list) / 3600. # convert to deg^2
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

def run_gen_all_obstime(basedir=None, listfile_dir=None, objs=None, output_dir=None, bad_obs=None, 
        bands=['B3','B4','B5','B6','B7','B8','B9','B10'], info_file=None, 
        suffix='good_imgs.txt.updated', exclude_aca=True,
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

    if bad_obs is not None:
        with open(bad_obs) as f:
            all_bad_obs = f.readlines()
        for i in range(len(all_bad_obs)):
            all_bad_obs[i] = all_bad_obs[i].strip()

    band_match = re.compile('_(?P<band>B\d{1,2})$')
    obj_match = re.compile('J\d{4}[-+]\d{4}')

    if objs is None:
        objs = []
        for item in os.listdir(basedir):
            if obj_match.match(item):
                objs.append(item)
            else:
                if debug:
                    print('Error load obj:', item)

    for i,obj in enumerate(objs):
        obj_exptime = {}
        for band in bands:
            obj_exptime[band] = 0
        print('index=', i, "obj:", obj)
            
        obj_dirname = os.path.join(basedir, obj)
        obj_output_dir = os.path.join(output_dir, obj)
        os.system('mkdir -p {}'.format(obj_output_dir))
        if listfile_dir:
            vis_list = []
            for band in bands:
                obj_band_selectfile = os.path.join(listfile_dir, obj, "{}_{}_{}".format(obj, band, suffix))
                if not os.path.isfile(obj_band_selectfile):
                    continue
                obj_band_list = gen_filenames(listfile=obj_band_selectfile, basedir=obj_dirname, exclude_aca=exclude_aca)
                vis_list.append(obj_band_list)
            vis_list = flatten(vis_list) 
        else:
            vis_list = gen_filenames(dirname=obj_dirname)
        obj_stat = spw_stat(vis=vis_list, debug=True,
                            savedata=True,
                            bands=bands,
                            filename=obj_output_dir+'/'+ obj+'.json', 
                            **kwargs)
    
        if info_file is not None:
            string_stat = '{:<12s}'.format(obj)
            for band in bands:
                string_stat += ' {:>8.2f}'.format(np.sum(obj_stat[band]['time']))
            with open(info_file, 'a+') as f_info:
                f_info.write(string_stat + '\n')

def run_gen_stats(basedir=None, listfile_dir=None, objs=None, bands=['B6','B7'], 
        suffix=['good_imgs.txt.updated', 'bad_imgs.txt.updated'],
        exclude_aca=True, debug=False, pickle_file=None,
        func_list=[],
        **kwargs):
    """The function to make some statistics like gen_all_obstime in a more general way
    """

    p_obj = re.compile('J\d+[+-]\d')
    p_obs = re.compile('uid___')
    band_match = re.compile('_(?P<band>B\d{1,2})$')
    obj_match = re.compile('J\d{4}[-+]\d{4}')

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

def run_gen_all_images(allcal_dir, obj_list=None, outdir='./', bands=['B6','B7'], exclude_aca=True, 
                  debug=False, niter=100, **kwargs):
    """fix the missing and wrong images for gen_all_image output
    

    default run: run_gen_all_image('all_image_dir', outdir='gen_all_image_run1')
    """
    filelist = []
    obj_match = re.compile('^J\d*[+-]\d*$')
    obs_match = re.compile('(?P<obsname>uid___\w*\.ms(\.split\.cal)?\.(?P<objname>[\s\w+-]+)_(?P<band>B\d+))')
    if obj_list is None:
        obj_list = []
        for obj in os.listdir(allcal_dir):
            if obj_match.match(obj):
                obj_list.append(obj)
    for obj in obj_list:
        print(obj)
        # start to go through each obj

        for obs in os.listdir(os.path.join(allcal_dir, obj)):
            if obs_match.match(obs):
                obs_band = obs_match.search(obs).groupdict()['band']
                infile = os.path.join(allcal_dir, obj, obs)
                if debug:
                    print(obs)
                if obs_band in bands:
                    outfile_fullpath = os.path.join(outdir, obj, obs_band)
                    os.system('mkdir -p {}'.format(outfile_fullpath))
                    
                    outfile_fullname = os.path.join(outfile_fullpath, obs+'.cont.auto.fits') 
                    if os.path.isfile(outfile_fullname):
                        if exclude_aca:
                            tb.open(infile + '/ANTENNA')
                            antenna_diameter = np.mean(tb.getcol('DISH_DIAMETER'))
                            tb.close()
                            if antenna_diameter < 12.0:
                                print("Removing aca image: {}".format(outfile_fullname))
                                os.system('rm -f {}'.format(outfile_fullname))
                        else:
                            if debug:
                                print(">> {} already exists.".format(outfile_fullname))
                    else:
                        if gen_image(infile, band=obs_band, outdir=outfile_fullpath, exclude_aca=exclude_aca, 
                                     debug=debug, niter=niter, **kwargs):
                            if debug:
                                print("Adding new image: {}".format(outfile_fullname))

def run_auto_classify_goodimags(imgsdir=None, objlist=None, outdir='./', bands=['B6','B7'], 
        debug=False, suffix='_good_imgs.txt', plot=True, savefig=True, **kwargs):
    """generate the good image list for all the calibrators

    default run: run_make_all_goodimags(imgs_dir='all_img_dir', basedir='science_ALMACAL', outdir='./') 
    """
    if imgsdir:
        obj_match = re.compile('^J\d*[+-]\d*$')
        for obj in os.listdir(imgsdir):
            if obj_match.match(obj):
                if objlist is not None:
                    if obj not in objlist:
                        continue
                    print(obj)
                obj_dir = os.path.join(outdir, obj)
                if not os.path.isdir(obj_dir):
                    os.system('mkdir -p {}'.format(os.path.join(outdir, obj)))
                for band in bands:
                    obj_band_path = os.path.join(imgsdir, obj, band)
                    good_imgs, bad_imgs = check_images(obj_band_path+'/*.fits', 
                            outdir=os.path.join(outdir, obj), plot=plot, savefig=savefig, 
                            band=band, basename=obj+'_'+band, debug=debug, **kwargs)
                    if debug:
                        print(good_imgs)

def run_manual_inspection(imagedir=None, outdir=None, objlist=None, bands=['B6','B7'], 
        suffix='imgs.txt'):
    """manually checking all the images
    """
    if isinstance(objlist, str):
        if os.path.isfile(objlist):
            objlist = gen_filenames(listfile=objlist)
    if not isinstance(objlist, (list, np.ndarray)):
        raise ValueError("Unsupported objlist!")
    for obj in objlist:
        obj_imagedir = os.path.join(imagedir, obj)
        for band in bands:
            obj_band_imagedir = os.path.join(obj_imagedir, band)
            obj_outdir = os.path.join(outdir, obj)
            goodfile = os.path.join(obj_outdir, "{}_{}_good_{}".format(obj, band, suffix))
            badfile = os.path.join(obj_outdir, "{}_{}_bad_{}".format(obj, band, suffix))
            if os.path.isfile(goodfile+'.updated'):
                print('{} of {} already done'.format(band, obj))
                continue
            else:
                print(obj_band_imagedir)
                print(">goodfile: {}\n>badfile: {}".format(goodfile, badfile))
                check_images_manual(imagedir=obj_band_imagedir, goodfile=goodfile, badfile=badfile, debug=False, ncol=1, nrow=3)

def run_make_all_goodimags(imgs_dir=None, objlist=None, bands=['B6','B7'], basedir=None, 
        outdir='./', debug=False, only_fits=True, update=True, 
        computwt=True, suffix='good_imgs.txt.updated', **kwargs):
    """generate the good image with updated list

    default run: run_make_all_goodimags(imgs_dir='all_img_dir', basedir='science_ALMACAL', outdir='./make_good_image') 
    """
    obj_match = re.compile('^J\d*[+-]\d*$')
    if objlist is None:
        objlist = []
        for obj in os.listdir(imgs_dir):
            if obj_match.match(obj):
                objlist.append(obj)
    for obj in objlist:
        obj_outdir = os.path.join(outdir, obj)
        print(obj)
        print(obj_outdir)
        if not os.path.isdir(obj_outdir):
            os.system('mkdir -p {}'.format(obj_outdir))
        for band in bands:
            good_image_file = os.path.join(obj_outdir, "{}_{}_{}".format(obj, band, suffix))
            concatvis = os.path.join(obj_outdir,"{}_{}_combine.ms".format(obj, band))
            if debug:
                print("good image file: {}\n concatvis_name: {}".format(good_image_file, concatvis))
            if os.path.isfile(good_image_file):
                good_image_fitsfile = concatvis+'.*.fits'
                # print('good_image_fitsfile: {}'.format(good_image_fitsfile))
                if len(glob.glob(good_image_fitsfile)) > 0:
                    print("Skip {} of {}, delete the fits file to continue...".format(band,obj))
                    continue
                else:
                    combined_vis = gen_filenames(listfile=good_image_file)
                    make_good_image(combined_vis, concatvis=concatvis, only_fits=only_fits, 
                            outdir=os.path.join(outdir, obj), computwt=computwt, 
                            basedir=os.path.join(basedir, obj), **kwargs)

def run_check_SMGs(basedir, objs=None, bands=['B6','B7'], suffix='combine.ms.auto.cont', 
                   resolutions=['0.3arcsec', '0.6arcsec'],
                   summary_file='summary.txt',
                   interative=False, outdir=None, continue_mode=True):
    """finding sources
    Adding a dot in the string: 
        resolutions = ['0.3arcsec','0.8arcsec']
    """
    if outdir is not None:
        if not os.path.isdir(outdir):
            os.system('mkdir -p {}'.format(outdir))
        summary_file = os.path.join(outdir, summary_file)
    else:
        summary_file = None
    if summary_file is not None:
        if os.path.isfile(summary_file):
            if continue_mode:
                summary = Table.read(summary_file, format='ascii')
                objs_finished = summary['obj']
            else:
                os.system('mv {} {}.old'.format(summary_file, summary_file))
                objs_finished = []
        else os.path.isfile(summary_file):
            print('Initializing the output file')
            with open(summary_file, 'w+') as f:
                f.write("obj")
                for band in bands:
                    for res in resolutions:
                        f.write(' ')
                        f.write(band+'_'+res)
                for band in bands:
                    f.write(' detection_{} goodfield_{}'.format(band, band))
                f.write(' is_SMG')
                f.write(' is_Jet')
                f.write('\n')

    obj_match = re.compile('^J\d*[+-]\d*$')
    if objs is None:
        objs = []
        for item in os.listdir(basedir):
            if obj_match.match(obj):
                objs.append(item)
    failed_files = []
    try:
        for obj in objs:
            if obj in objs_finished:
                print("{} already done!".format(obj))
                continue
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
                obj_sourcefound = {}
                for band in bands:
                    for res in resolutions:
                        obj_sourcefound[band+'_'+res] = 0
                # make a summary plot
                fig, ax = plt.subplots(len(bands), len(resolutions), 
                                       figsize=(4.2*len(resolutions),4*len(bands))) 
                fig.suptitle(obj)
                #for img in imgs:
                # sources_number = {}
                for i,band in enumerate(bands):
                    obj_band_dir = os.path.join(basedir, band, 'make_good_images', obj)
                    for j,res in enumerate(resolutions):
                        if len(bands) > 1:
                            ax_select = ax[i,j]
                        else:
                            ax_select = ax[max(i,j)]
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
                        savefile = image_name + '.source_found.txt'
                        figname = image_name + '.png'
                        sources_found = source_finder(image_fullpath, outdir=obj_outdir, 
                                ax=ax_select, pbcor=True)
                        ax_select.set_title('{} {}'.format(band, res))
                        #try:
                        #    sources_found = source_finder(image_fullpath, outdir=obj_outdir, 
                        #            ax=ax[i,j], pbcor=True)
                        #except:
                        #    print("Error found for {}".format(image_name))
                        #    failed_files.append(image_name)

                        if len(sources_found) > 0 and obj_summary is not None:
                            obj_summary.write('# {} {} {}\n'.format(obj, band, res))
                            obj_sourcefound['{}_{}'.format(band, res)] = len(sources_found)
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
                for band in bands:
                    for res in resolutions:
                        found_string += ' '+str(obj_sourcefound[band+'_'+res])
                print(found_string)
                # save figure
                fig.subplots_adjust(wspace=0.2, hspace=0.2)
                fig.savefig(summary_plot, bbox_inches='tight', dpi=400)
                if summary_file: 
                    SMG_input = 0
                    Jet_input = 0
                    detections = {}
                    goodfields = {}
                    for band in bands:
                        goodfields[band] = 0
                        detections[band] = 0
                    if interative:
                        plt.show()
                        print("Single Band:\n") 
                        print("Detection: 0)None 1)point source 2)Partially 3)enlongated jet -1)Not Sure")
                        print("Usable: 0)No 1)Full 2)exclude center -1)re-imaging -2)images reselection")
                        print("General Classification")
                        print("Is SMG: 0)No 1)Yes 2)Multiple 3)Radio Source -1)Not Sure")
                        print("Is Jet: 0)No 1)Compact 2)Enlongated -1)Not Sure")
                        for band in bands:
                            detection_input=int(raw_input("Dection in Band:{} (integer, 0,1,2,3,-1) [0]?: ".format(band)) or 0)
                            goodfield_input=int(raw_input("Usable for Band:{} (integer, 0,1,2,-1,-2) [0]?: ".format(band)) or 0)
                            detections[band] = detection_input
                            goodfields[band] = goodfield_input
                        if len(bands) > 1:
                            SMG_input = int(raw_input("Is SMG? (integer, 0,1,2,3,-1) [0]: ") or 0)
                            Jet_input = int(raw_input("Is Jet? (integer, 0,1,2,-1) [0]: ") or 0)
                        plt.close()
                    for band in bands:
                        found_string += ' {} {}'.format(detections[band], goodfields[band])
                    found_string += ' {}'.format(SMG_input)
                    found_string += ' {}'.format(Jet_input)
                    with open(summary_file, 'a+') as f:
                        f.write("{}\n".format(found_string)) 
                else:
                    next_one = int(raw_input("Next one [1/0] [1]") or 1)
                    if next_one == 1:
                        continue
                    else:
                        return 0
                plt.close()
    except KeyboardInterrupt:
        return 0

def run_search_bands_detections(basedir=None, objs=None, bands=['B6','B7'], outdir='./', debug=False):
    for obj in objs:
        for band in bands:
            search_band_detection(basedir=os.path.join(basedir, obj), band=band, 
                                  outdir=os.path.join(outdir, obj, band), debug=debug)

def run_make_simulations(imagedir=None, objs=None, bands=['B6','B7'], outdir='./', 
        resolution=None, repeat=1000, ):
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
        obj_outdir = os.path.join(outdir, obj)
        for band in bands:
            imagefile = os.path.join(obj_outdir,
                    '{}_{}_combine.ms.auto.cont{}.image.fits'.format(obj, band, resolution_string))
            obj_outdir = os.path.join(outdir, obj)
            if os.path.isfile(imagefile):
                gen_sim_images(imagefile=imagefile, outdir=outdir, snr=(0.1,20), repeat=repeat)
                calculate_sim_images(outdir, baseimage=imagefile, savefile=os.path.join(
                    outdir,'{}_simulation.txt'.format(obj)), plot=False, repeat=repeat, threshold=5, 
                    second_check=False, snr_mode='peak')

def run_calculate_effarea(imagedir, flux=np.linspace(0.1, 1, 10),  objs=None, bands=['B6','B7'], 
        suffix='combine.ms.auto.cont', resolutions=['', 'uvscale1.5', 'uvscale2.0'], 
        outdir='./', snr_threshold=5.0, savefile=None):
    """calculate the effecitve area for all the usable fields
    """
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
    obj_match = re.compile('^J\d*[+-]\d*$')
    if objs is None:
        objs = []
        for item in os.listdir(basedir):
            if obj_match.match(obj):
                objs.append(item)

    effarea = Table()
    effarea['flux'] = flux
    for band in bands:
        for res in resolutions:
            effarea[band+'_'+res] = np.zeros_like(flux)
    for obj in objs:
        obj_dir = os.path.join(basedir, obj)
        for band in bands:
            for res in resolutions:
                if res == '':
                    res_string = ''
                else:
                    res_string = res+'.'
                image = "{}_{}_{}.{}image.fits".format(obj, band, suffix, res_string)
                image_fullpath = os.path.join(obj_dir, image)
                image_pbcorr = image.replace('image', 'pbcor.image')
                image_pbcorr_fullpath = os.path.join(obj_dir, image_pbcorr)
                if os.path.isfile(image_fullpath):
                    _, objarea = calculate_effectivearea(flux, snr_threshold=snr_threshold, 
                                    images=[image_fullpath,], images_pbcorr=[image_pbcorr_fullpath,])
                    effarea[band+'_'+res] += objarea
    if savefile:
        savefile_fullpath = os.path.join(outdir, savefile)
        effarea.write(savefile_fullpath, format='ascii')

def run_differential_number_counts(flist, detections_file=None, band='B6', effective_area_file=None, 
        simulation_folder=None, ax=None):
    """

    flist: flux list, [0.2, 0.6, 1.0]
    detections_file: the file with detection objs and their flux and flux snr, file format:
        obj B6 B6_SNR B7 B&_SRN
    band: the selected band
    effective_area_file: the file include two columns: flux and its effective area
    simulation_folder: the folder generated by run_make_simulations
    ax: the plot axes
    """
    Ni_list = []
    tab = Table.read(detections_file, format='ascii')
    effarea = np.loadtxt(effective_area_file, skiprows=1)
    cs_effarea = CubicSpline(effarea[:, 0], effarea[:,1])
    for item in tab:
        flux = item[band]
        if flux < 0:
            continue
        obj = item['obj']
        snr = item[band+'_SNR']
        print(obj, flux, snr)
        snr_list, apert_boost_list, comp_list, fake_rate_list = plot_sim_results(
                jsonfile=os.path.join(simulation_folder, obj, obj+'_simulation.txt'), 
                snr=np.arange(0.2, 10, 0.2), plot=False)
        cs_comp = scipy.interpolate.interp1d(snr_list, comp_list) 
        def cs_comp2(snr):
            if snr<10:
                return cs_comp(snr)
            else:
                return 1.0
        print('effarea:', cs_effarea(flux)/3600.)
        print('completeness', cs_comp2(snr))
        Ni = 1 / (cs_effarea(flux)/3600.) / cs_comp2(snr)
        #print(Si, Ni)
        Ni_list.append([flux, Ni])
    Ni_array = np.array(Ni_list)
    print(Ni_array)

    NN = []
    for f in flist:
        print('flux:', f)
        print(flux)
        NN.append(np.sum(Ni_array[:,1][Ni_array[:,0] >= f]))
    print(flist, NN)
    if ax is None:
        fig = plt.figure(figsize=(6, 8))
        ax = fig.add_subplot(111)
    ax.loglog(flist, NN, 'ro')
    ax.loglog(flist, NN, 'r--')
