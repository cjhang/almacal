# This file include all the tasks for ALMACAL number counts projects

# The task list include:
# 1. gen_obstime: generate the observational time for the whole dataset
# 2. gen_dirty_image: generate the dirty image for all the observations of one calibrator
# 3. show_image: interative way to inspect the images manually 
import os
import glob
import re
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

import analysisUtils as au
# from analysisUtils import tbtool
# from imaging_utils import make_cont_img

tb = au.tbtool()

def gen_filenames(dirname=None, listfile=None, basedir=None, debug=False):
    """generate all the valid files names

    """
    if dirname:
        obj_match = re.compile('^J\d*[+-]\d*$')
        obs_match = re.compile('(?P<obsname>uid___\w*\.ms(\.split\.cal)?\.(?P<objname>[\s\w+-]+)_(?P<band>B\d+))')
        file_list = []
        if not os.path.isdir(dirname):
            raise ValueError('Invalid directory!')
        for item in os.listdir(dirname):
            if debug:
                print(item)
            if obj_match.match(item):
                obj_path = os.path.join(dirname, item)
                for obs in os.listdir(obj_path):
                    if obs_match.match(obs):
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

def gen_image(vis=None, band=None, outdir='./', exclude_aca=False, check_calibrator=False, debug=False, update_raw=False, **kwargs):
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

    if exclude_aca:
        tb.open(vis + '/ANTENNA')
        antenna_diameter = np.mean(tb.getcol('DISH_DIAMETER'))
        tb.close()
        if antenna_diameter < 12.0:
            if debug:
                print("Excuding data from {}".format(antenna_diameter))
            return False
    try:
        make_cont_img(vis=vis, clean=True, myimagename=myimagename, outdir=outdir, niter=1000, **kwargs)
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

def gen_images(allcal_dir=None, vis=None, outdir='./', bands=None, exclude_aca=True, 
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
    filelist = []
    if allcal_dir:
        for obj in os.listdir(allcal_dir):
            if debug:
                print(obj)
            obj_match = re.compile('^J\d*[+-]\d*$')
            if obj_match.match(obj):
                filelist.append(os.path.join(allcal_dir, obj))
                continue
    elif vis:
        if isinstance(vis, str):
            filelist = [vis,]
        elif isinstance(vis, list):
            filelist = vis

    for infile in filelist:
        if bands is not None:
            for band in bands:
                outdir = os.path.join(outdir, band)
                os.system('mkdir -p {}'.format(outdir))
                gen_image(infile, band=band, outdir=outdir, exclude_aca=exclude_aca, **kwargs)
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
    
    print("Totally {}% of data have been selected.".format(100.*select_num/total_num))
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
    
def read_flux(obs, allflux_file=None, strick_mode=False):
    obs_filename = os.path.basename(obs)
    if allflux_file is None:
        allflux_file = os.path.join(os.path.expanduser('~'), 'Documents/projects/almacal/data/allcal.fluxval')
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
            myflux = read_flux(obs, allflux_file=allflux_file)
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

def copy_ms(basedir, outdir, selectfile=None, debug=False, time_select=False, 
            start_time='2010-01-01T00:00:00', end_time='2050-01-01T00:00:00', 
            select_band=None):
    p_obs = re.compile('uid___')
    band_match = re.compile('_(?P<band>B\d{1,2})$')
    selectfiles_list = []
    start_time = Time(start_time)
    end_time = Time(end_time)
    if selectfile is not None:
        with open(selectfile) as f:
            selectfiles_readlines = f.readlines()
        for item in selectfiles_readlines:
            # remove the '\n' at the end of item
            selectfiles_list.append(item.strip())

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
            os.system('cp -r {} {}'.format(os.path.join(basedir, obs), outdir))

def gaussian(x, u0, amp, std):
    return amp*np.exp(-0.5*((x-u0)/std)**2)

def check_image(img, plot=False, radius=6, debug=False, sigmaclip=True, check_flux=True, minimal_fluxval=0.001, outlier_frac=0.02, 
                gaussian_deviation=0.25, central_median_deviation=1.0, central_mean_deviation=2.0):
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

    ny, nx = data.shape[-2:]
    masked_data = np.ma.masked_invalid(data.reshape(ny, nx))
    mask_invalid = masked_data.mask
    data4hist = masked_data[~mask_invalid]
    # Testing the sigma clip
    if sigmaclip:
        clipped_data, clip_low, clip_up = scipy.stats.sigmaclip(data4hist, 5, 5)
        data4hist = clipped_data    
    # statistics the noisy of the image
    # hist, bins = np.histogram(masked_data.data[~masked_data.mask], bins=100)
    hist, bins = np.histogram(data4hist, bins=100)
    bins_mid = (bins[:-1] + bins[1:])*0.5
    p0 = (0, 1, 1e-4) # mean, max and std
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
    upper_3sigma = mean + 3.0*sigma
    lower_5sigma = mean - 5.0*sigma
    upper_5sigma = mean + 5.0*sigma
    ## checking the fitting
    # the fraction of the 3 sigma outlier, theoretical gaussian value is 
    n_outlier_3sigma = np.sum(hist[bins_mid < lower_3sigma])
    percent_outlier_3sigma = 1.0 * n_outlier_3sigma / np.sum(hist) 
    percent_outlier_1sigma = 1.0 * np.sum(hist[bins_mid < lower_1sigma]) / np.sum(hist) 
    percent_outlier_2sigma = 1.0 * np.sum(hist[bins_mid < lower_2sigma]) / np.sum(hist) 
    # calculating the deviation from Gaussian
    deviation_1sigma = np.abs((percent_outlier_1sigma - 0.1587)/0.1587)
    deviation_2sigma = np.abs((percent_outlier_2sigma - 0.0228)/0.0228)
    deviation_3sigma = np.abs((percent_outlier_3sigma - 0.0014)/0.0014)

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
    
    n_outlier_central = np.sum(hist_center[bins_mid_center>upper_5sigma])
    percent_outlier_central = 1.0*n_outlier_central/np.sum(hist_center)
    
    if debug:
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print('>> Checking the fitting of noise:')
        print('mean: {}, std:{}'.format(mean, sigma))
        # print('number of 3 sigma outlier: {}'.format(n_outlier_3sigma))
        print('fraction of 1 sigma outlier: {:.4f}%. [theoretical: 15.87]'.format(percent_outlier_1sigma*100.))
        print('deviation of 1 sigma: {:.4f}%.'.format(deviation_1sigma*100.))
        print('fraction of 2 sigma outlier: {:.4f}%. [theoretical: 2.28%]'.format(percent_outlier_2sigma*100.))
        print('deviation of 2 sigma: {:.4f}%.'.format(deviation_2sigma*100.))
        print('fraction of 3 sigma outlier: {:.4f}%. [theoretical: 0.14%]'.format(percent_outlier_3sigma*100.))
        print('deviation of 3 sigma: {:.4f}%.'.format(deviation_3sigma*100.))
        print('>> Statistics of central region')
        print("central mean: {}".format(mean_central))
        print("central median: {}".format(median_central))
        print('deviation of central median: {:.4f}.'.format(median_central/sigma))
        print('number of 5sigma outlier of central region: {}'.format(n_outlier_central))
        print('fraction of 5sigma outlier of central region: {}'.format(percent_outlier_central))

    if plot:
        # show the image
        fig = plt.figure(figsize=(14, 6))
        ax = fig.add_subplot(121)
        scale = np.abs(header['CDELT1'])*3600
        x_index = (np.arange(0, nx) - nx/2.0) * scale
        y_index = (np.arange(0, ny) - ny/2.0) * scale
        x_map, y_map = np.meshgrid(x_index, y_index)
        ax.pcolormesh(x_map, y_map, masked_data)
        ax.text(0, 0, '+', color='r', fontsize=24, fontweight=100, horizontalalignment='center',
                verticalalignment='center')
        circle = patches.Circle((0, 0), radius=bmaj*3600*radius*0.5, facecolor=None, fill=None, edgecolor='red', linewidth=2, alpha=0.5)
        ellipse = patches.Ellipse((0.8*np.min(x_index), 0.8*np.min(y_index)), width=bmin*3600, height=bmaj*3600, angle=bpa, facecolor='orange', edgecolor=None, alpha=0.8)
        ax.add_patch(circle)
        ax.add_patch(ellipse)
        ax.set_xlabel('RA [arcsec]')
        ax.set_ylabel('Dec [arcsec]')
        
        # niose statistics
        ax = fig.add_subplot(122)
        ax.step(bins_mid, hist/amp_scale, where='mid', color='b', label='Noise Distribution')
        ax.step(bins_mid, hist_center/amp_scale_center, where='mid', color='orange', label='Central Noise Distribution')
        ax.plot(bins_mid, hist_fit/amp_scale, color='r', label='Gaussian Fitting')
        ax.vlines(upper_3sigma, 0, 2.0, color='k', label=r'3$\sigma$ boundary')
        ax.vlines(lower_3sigma, 0, 2.0, color='k')
        ax.vlines(upper_5sigma, 0, 2.0, color='k', lw=4, label=r'5$\sigma$ upper boundary')
        ax.set_xlabel('Flux density [Jy/beam]')
        ax.set_ylabel('Normalized Pixel numbers')
        ax.legend()
        plt.show()
    
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
            fluxval = read_flux(uidname)
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
    strick_mode = False
    # comparing the noise distribution with Gaussian
    for deviation in [deviation_1sigma, deviation_2sigma]:
        if deviation > gaussian_deviation:
            if debug:
                print("Rjected, non-Gaussian noise")
            return False
    if strick_mode:
        if deviation_3sigma > gaussian_deviation:
            if debug:
                print("Rejected, non-Gaussian noise at 3sigma boudary!\n")
            return False
        if percent_outlier_central >= outlier_frac:
            if debug:
                print("Rejected, large central residual!\n")
            return False
        if mean_central > threshold_mean:
            if debug:
                print("Rejected, large central mean value!\n")
            return False
        if debug:
            print("\n")
    return True

    if 0:
        imghead = imhead(img, mode='list')
        imgstat = imstat(img)

        unit = imghead['bunit']
        rms = imgstat['rms']

        # define the central region
        beam_major = imghead['beammajor']['value']*u.Unit(imghead['beammajor']['unit'])
        radius = 5 
        region_string = 'circle[[{},{}], {}]'.format(str(imghead['crval1'])+imghead['cunit1'],
                                                     str(imghead['crval2'])+imghead['cunit2'],
                                                     radius*beam_major)
        imgstat_central = imstat(img, region=region_string)
        central_max = imgstat_central['max']
        central_mean = imgstat_central['mean']
        central_sum = imgstat_central['sum']

        print("For whole image: RMS={} {}, sum={}".format(rms, unit, imgstat['sum']))
        print("For the central region: max:{},   mean:{},   sum:{}".format(central_max, central_mean, central_sum))
        print("For the ratio: max:{},   mean:{},   sum:{}".format(central_max/rms, central_mean/rms, central_sum/rms))

def check_images(imgs, outdir=None, basename='', debug=False, **kwargs):
    """wraps up check_image to handle multiple images
    """
    if isinstance(imgs, str):
        all_files = glob.glob(imgs)
    elif isinstance(imgs, list):
        all_files = imgs

    p_uidimg = re.compile('(?P<uidname>uid___\w+.ms.split.cal.J\d+[+-]\d+_B\d+).cont.auto.fits')

    good_imgs = []
    bad_imgs = []
    for img in all_files:
        if debug:
            print("img: {}".format(img))
        # continue
        if p_uidimg.search(img):
            uidname = p_uidimg.search(img).groupdict()['uidname']
        else:
            uidname = None
        # 
        if check_image(img, debug=debug, **kwargs):
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
        if len(good_imgs) > 0:
            with open(os.path.join(outdir, basename+"_good_imgs.txt"), "w") as good_outfile:
                good_outfile.write("\n".join(str(item) for item in good_imgs))
        if len(bad_imgs) > 0:
            with open(os.path.join(outdir, basename+"_bad_imgs.txt"), "w") as bad_outfile:
                bad_outfile.write("\n".join(str(item) for item in bad_imgs))
    return good_imgs, bad_imgs

def make_good_image(vis=None, basename='', basedir=None, outdir='./', tmpdir='./', concatvis=None, debug=False, only_fits=False,
                    niter=100, clean=True, pblimit=-0.01, fov_scale=2.0, uvtaper_scale=[1.0, 2.0], **kwargs):
    """make the final good image with all the good observations
    """
    if len(vis) < 1:
        return False
    if basedir:
        vis_fullpath = []
        for v in vis:
            vis_fullpath.append(os.path.join(basedir, v))
        vis = vis_fullpath
    if debug:
        print(vis)
    if concatvis is None:
        concatvis = os.path.join(tmpdir, basename+'combine.ms')
    concat(vis=vis, concatvis=concatvis)
    make_cont_img(vis=concatvis, myimagename=concatvis+'.auto.cont', clean=clean, niter=niter, pblimit=pblimit, 
                  fov_scale=fov_scale, uvtaper_scale=uvtaper_scale, debug=debug, **kwargs)

    if only_fits:
        for i in glob.glob(concatvis+'.*.image'):
            exportfits(imagename=i, fitsimage=i+'.fits')
        rmtables(concatvis+'*')

def gen_fake_images(vis, known_file=None, flux=np.arange(0.1,10,0.5), fov_scale=1.5, outdir='./', basename=None):
    """generate the fake images with man-made sources
    """
    # read information from vis 
    spw_specrange = read_spw(vis)
    freq_mean = np.mean(spw_specrange) # in GHz
    tb.open(vis + '/ANTENNA')
    antenna_diameter_list = tb.getcol('DISH_DIAMETER')
    tb.close()
    antenna_diameter = np.max(antenna_diameter_list) * u.m
    wavelength = const.c / (freq_mean * u.GHz) # in um
    fov = (fov_scale * 1.22 * wavelength / antenna_diameter * 206265).decompose()
    print('fov', fov)
    print('radius', 0.5*fov)
    if basename is None:
        basename = os.path.basename(vis)
    for f in flux:
        add_random_sources(vis, n=20, radius=0.5*fov, outdir=outdir, 
                            basename=basename+'.{}mJy'.format(f), flux=f, known_file=known_file,
                            uvtaper_scale=None)

def calculate_completeness(objfolder, image=None, known_file=None, obj=None, band=None, basename=None,):
    """simulation the completeness of source finding algorithm
    """
    # image statistics
    im_info = imstat(image)
    rms = im_info['rms']
    print(rms)
    
    flux_match = re.compile('(?P<obj>J\d*[+-]\d*)_(?P<band>B\d+)_combine.ms.(?P<flux>\d+.\d+)mJy')
    if basename is None:
        basename = os.path.join(objfolder, '{obj}_{band}_combine.ms'.format(obj=obj, band=band))
        print('basename', basename)
    all_fake_images = glob.glob(basename+'.*.fits')

    flux_list = []
    flux_input_list = []
    flux_peak_list = []
    completeness_list = []
    flux_found_list = []
    for img in all_fake_images:
        search_result = flux_match.search(img)
        if search_result:
            flux = search_result.groupdict()['flux']
            #print('flux:', flux)
            sf_return = source_finder(img, sources_file=basename+'.{}mJy.txt'.format(flux), known_file=known_file)
            if sf_return == 0:
                continue
            flux_input, flux_auto, sources_input_found = sf_return 
            #print('flux_input', flux_input)
            #print('flux_auto', flux_auto)
            #print('sources_input_found', sources_input_found)
            
            flux_list.append(float(flux))#sources_input_found[0])
            flux_input_list.append(sources_input_found[0])
            flux_peak_list.append(sources_input_found[2])
            flux_found_list.append(sources_input_found[1])
            completeness_list.append(1.*len(sources_input_found[1])/len(flux_input))
    # print(flux_list)
    # print(flux_peak_list)
    # print(flux_found_list)
    # print(completeness_list)
    if True:
        fig = plt.figure(figsize=(8, 4))
        ax = fig.add_subplot(1,2,1)
        ax.set_xlabel('SNR')
        ax.set_ylabel(r'$S_{\rm out}/S_{\rm in}$')
        for i in range(len(flux_list)):
            ax.plot(flux_peak_list[i]/rms, flux_found_list[i][:,0]/flux_list[i], 'k.')
            ax.plot(flux_peak_list[i]/rms, flux_found_list[i][:,1]/flux_list[i], 'r.')
        # ax.plot(np.array(flux_peak_list)/rms, np.array(flux_found_list)/np.array(flux_list), 'o')
        ax = fig.add_subplot(1,2,2)
        f_mean = lambda x: np.mean(x)
        ax.plot(map(f_mean, flux_peak_list)/rms, completeness_list, 'o')
        # ax.plot(np.array(flux_peak_list)/rms, completeness_list, 'o')
        ax.set_xlabel('SNR')
        ax.set_ylabel(r'Completeness')
        ax.set_xlim((3., 12))
        ax.set_ylim((0, 1.2))
        plt.show()
    return flux_list, flux_peak_list, flux_found_list, completeness_list

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# The ALMA run automatic pipeline section #
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def run_gen_all_obstime(base_dir=None, output_dir=None, bad_obs=None, 
        info_file=None, **kwargs):
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

    for i,obj in enumerate(os.listdir(base_dir)):
        obj_exptime = {'B3':0, 'B4':0, 'B5':0, 'B6':0, 
                    'B7':0, 'B8':0,  'B9':0, 'B10':0}
        if not obj_match.match(obj):
            print('Error load obj:', obj)
            continue
        print('index=', i, "obj:", obj)
            
        obj_dirname = base_dir +'/'+ obj
        obj_output_dir = output_dir + '/' + obj
        os.system('mkdir {}'.format(obj_output_dir))
        obj_stat = spw_stat(obj_dirname,
                            savedata=True, 
                            filename=obj_output_dir+'/'+ obj+'.json', 
                            **kwargs)
    
        if info_file is not None:
            with open(info_file, 'a+') as f_info:
                f_info.write('{:<12s} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f}\n'.format(
            obj, 
            np.sum(obj_stat['B3']['time']),
            np.sum(obj_stat['B4']['time']),
            np.sum(obj_stat['B5']['time']),
            np.sum(obj_stat['B6']['time']),
            np.sum(obj_stat['B7']['time']),
            np.sum(obj_stat['B8']['time']),
            np.sum(obj_stat['B9']['time']),
            np.sum(obj_stat['B10']['time']),
            ))

def run_gen_oteo2016_data(basedir, outdir, select_band=['B6', 'B7'], debug=False):
    """copy the data used in oteo2016

    default run: run_gen_oteo2016_data('/science_ALMACAL/data', outdir='./')
    """
    p_obj = re.compile('J\d+[+-]\d+')
    for obj in os.listdir(basedir):
        if not p_obj.match(obj):
            print('Error load obj:', obj)
            continue
        obj_input_folder = os.path.join(basedir, obj)
        obj_output_folder = os.path.join(outdir, obj)
        print("Copying {}".format(obj))
        copy_ms(obj_input_folder, obj_output_folder, time_select=True, 
                end_time='2015-07-01T00:00:00',
                select_band=select_band, debug=debug)

def run_gen_all_image(allcal_dir, outdir='./', bands=['B6','B7'], exclude_aca=True, 
                  debug=False, **kwargs):
    """fix the missing and wrong images for gen_all_image output
    

    default run: run_gen_all_image('all_image_dir', outdir='gen_all_image_run1')
    """
    filelist = []
    obj_match = re.compile('^J\d*[+-]\d*$')
    obs_match = re.compile('(?P<obsname>uid___\w*\.ms(\.split\.cal)?\.(?P<objname>[\s\w+-]+)_(?P<band>B\d+))')
    for obj in os.listdir(allcal_dir):
        if obj_match.match(obj):
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
                                         debug=debug, **kwargs):
                                if debug:
                                    print("Adding new image: {}".format(outfile_fullname))

def run_make_all_goodimags(imgs_dir=None, objlist=None, good_imgs_file=None, basedir=None, make_image=False, outdir='./', 
                           debug=False, only_fits=False, update=True, **kwargs):
    """generate the good image list for all the calibrators

    default run: run_make_all_goodimags(imgs_dir='all_img_dir', basedir='science_ALMACAL', make_image=True, outdir='./', only_fits=True) 
    """
    if imgs_dir:
        obj_match = re.compile('^J\d*[+-]\d*$')
        for obj in os.listdir(imgs_dir):
            if obj_match.match(obj):
                if objlist is not None:
                    if obj not in objlist:
                        continue
                    print(obj)
                obj_dir = os.path.join(outdir, obj)
                if os.path.isdir(obj_dir):
                    if len(os.listdir(obj_dir)) > 0:
                        if debug:
                            print("Skip {}".format(obj))
                        continue
                else:
                    os.system('mkdir -p {}'.format(os.path.join(outdir, obj)))
            for band in os.listdir(os.path.join(imgs_dir, obj)):
                if update:
                    combined_file = os.path.join(outdir, obj, obj+'_'+band+'_combine.ms')
                    good_image_file = os.path.join(outdir, obj, obj+'_'+band+'_good_imgs.txt')
                    if os.path.isdir(combined_file):
                        print("\n\n'n>>>>>>>>>>>Find combined file>>>>>>>>> \n\n")
                        make_good_image(combined_file, basename=obj+'_'+band+'_', basedir=os.path.join(basedir,obj), 
                                        tmpdir=os.path.join(outdir,obj), only_fits=only_fits, debug=debug)
                    elif os.path.isfile(good_image_file):
                        print("\n\n'n>>>>>>>>>>>Find combined file>>>>>>>>> \n\n")
                        good_imgs = []
                        with open(good_images_file) as f:
                            good_imags_lines = f.readlines()
                        for gi in good_imags_lines:
                            good_imgs.append(gi.strip())
                        make_good_image(good_imgs, basename=obj+'_'+band+'_', basedir=os.path.join(basedir,obj), 
                                        tmpdir=os.path.join(outdir,obj), only_fits=only_fits, debug=debug)

                obj_band_path = os.path.join(imgs_dir, obj, band)
                good_imgs, bad_imgs = check_images(obj_band_path+'/*.fits', outdir=os.path.join(outdir, obj), 
                                                    basename=obj+'_'+band, debug=debug, **kwargs)
                if debug: 
                    print(good_imgs)
                if make_image:
                        make_good_image(good_imgs, basename=obj+'_'+band+'_', basedir=os.path.join(basedir,obj), 
                                        tmpdir=os.path.join(outdir,obj), only_fits=only_fits, debug=debug)

def run_gen_fake_images(basedir, bands=['B7',], outdir='./tmp', flux=np.arange(0.1,10,0.5)):
    obj_match = re.compile('^J\d*[+-]\d*$')
    for obj in os.listdir(basedir):
        if obj_match.match(obj):
            for band in bands:
                objfolder = os.path.join(basedir, obj)
                vis_combined = os.path.join(objfolder, '{}_{}_combine.ms'.format(obj, band))
                if os.path.isdir(vis_combined):
                    gen_fake_images(vis=vis_combined, outdir=os.path.join(outdir, obj), flux=flux,
                                    known_file=vis_combined+'.auto.cont.image.fits.source_found.txt')

def run_calculate_completeness():
    pass

def run_find_source(basedir, summary_file=None):
    obj_match = re.compile('^J\d*[+-]\d*$')
    for obj in os.listdir(basedir):
        if obj_match.match(obj):
            print('>>>>> {}'.format(obj))
            sources_nums = []
            imgs = glob.glob(os.path.join(basedir, obj, '*.fits'))
            for img in imgs:
                print(img)
                savefile = img+ '.source_found.txt'
                sources_found = source_finder(img, savefile=None)
                if sources_found != 0:
                    sources_nums.append(len(sources_found))
        if summary_file:
            with open(summary_file, 'a+') as f:
                f.write("{}, {}\n".format(obj, sources_nums))


