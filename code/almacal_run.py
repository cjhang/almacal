# This file include all the tasks for ALMACAL number counts projects

# The task list include:
# 1. gen_obstime: generate the observational time for the whole dataset
# 2. gen_dirty_image: generate the dirty image for all the observations of one calibrator
# 3. show_image: interative way to inspect the images manually 

import glob
import re
import analysisUtils as au
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
from matplotlib import pyplot as plt


def gen_obstime(base_dir=None, output_dir=None, bad_obs=None, info_file=None, 
                **kwargs):
    """generate the on-source time and spw distribution for the whole almacal dataset
    
    Params:
        base_dir: the root folder contains all the measurements
        output_dir: the root folder for placing the results of each calibrator,
                    including the json file and the plots
        bad_obs: the file contains unusable observation
        info_file: the file contains the general information for all the calibrators
        **kwargs: support the addition parameters of `spw_stat`
    

    default run:
    obj_output_dir = output_dir + '/' + obj + '/'
    os.system('mkdir {}'.format(obj_output_dir))
    spw_stat('/', plot=True, showfig=False, figname=obj_output_dir + obj+'.pdf', savedata=True, filename=obj_output_dir + obj+'.json')
    """

    p_obj = re.compile('J\d+[+-]\d')
    p_obs = re.compile('uid___')
    
    # bad_obs = '/Users/jchen/Desktop/projects/almacal/data/broken_obs.txt'
    # almacal_info_file = '/Users/jchen/Desktop/projects/almacal/data/almacal_timeOnSource.txt'
    # all_obs = Table.read(almacal_info_file, format='ascii')
    # all_obs.sort(['B6', 'B7'])
    # all_obs.reverse()

    if bad_obs is not None:
        with open(bad_obs) as f:
            all_bad_obs = f.readlines()
        for i in range(len(all_bad_obs)):
            all_bad_obs[i] = all_bad_obs[i].strip()

    band_match = re.compile('_(?P<band>B\d{1,2})$')
    obj_match = re.compile('J\d{4}[-+]\d{4}')

    for i,obj in enumerate(os.listdir(base_dir)):
    # for i,obj in enumerate(all_obs['obj']):
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
 
def gen_image(obj, band=None, outdir='./', **kwargs):
    """make images for one calibrator on all or specific band

    Params:
    obj: the folder that contains all the observation of one calibrator
    band: the band want to imaged
    outdir: the directory where all the fits image will be placed
    **kwargs: the additional parameters supported by make_cont_img

    """
    for obs in os.listdir(obj):
        if band is not None:
            band_match = re.compile('_(?P<band>B\d{1,2})$')
            if band_match.search(obs):
                obs_band = band_match.search(obs).groupdict()['band']
                if obs_band != band:
                    continue
        basename = os.path.basename(obs)
        myimagename = os.path.join(outdir, basename + '.cont.auto')
        try:
            make_cont_img(vis=obj+'/'+obs, dirty_image=True, myimagename=myimagename, outdir=outdir, **kwargs)
        except:
            print("Error in imaging {}".format(obj))
        exportfits(imagename=myimagename+'.image', fitsimage=myimagename+'.fits')
        rmtables(tablenames=myimagename+'.*')

def gen_all_image(allcal_dir, outdir='./', bands=['B6','B7'], **kwargs):
    """generate the images of all calibrators

    Params:
        allcal_dir: the root directory contains all the measurements of calibrators
        outdir: the output directory
        bands: the bands to be imaged
        **kwargs: the additional parameters of make_cont_img


    default run: gen_all_image('/science_ALMACAL/data', '/tmp/all_images')
    """
    for obj in os.listdir(allcal_dir):
        obj_match = re.compile('^J\d*[+-]\d*$')
        if not obj_match.match(obj):
            continue
        for band in bands:
            os.system('mkdir -p {}'.format(os.path.join(outdir, obj, band)))
            gen_image(os.path.join(allcal_dir, obj), band=band, 
                      outdir=os.path.join(outdir, obj, band),)

def show_images(fileglob, mode='auto', nrow=3, ncol=3, savefile=None):
    """show images in an interative ways, and record the input from inspector
    
    Parameters:
        fileglob: the match patter of the image filenames, like: './*.image'
    """
    all_files = glob.glob(fileglob)
    print("Find {} files".format(len(all_files)))
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
            except:
                print("Error in reading: {}".format(all_files[i+j]))
                continue
            ax = fig.add_subplot(nrow, ncol, j+1)#, projection=wcs2, slices=(0, 0, 'x', 'y'))
            # ax.text(10, 10, str(j), fontsize=20)
            ax.set_title(str(j+1))
            ax.imshow(imagedata[0,0,:,:], origin='lower')#, cmap='viridis')
            ax.set_xlabel('RA')
            ax.set_ylabel('Dec')
        # show the image and record the 
        plt.show()
        print('Input the index of images (1-9), seperate with comma:')
        try:
            idx_input = input()
            if idx_input == 0:
                print("Currently at {}/{}".format(i, len(all_files)))
                break
            if isinstance(idx_input, int):
                idx_input = [idx_input]
            for ind in idx_input:
                all_select.append(all_files[i+ind-1])
        except:
            plt.close('all')
            continue
        # plt.clf()
        plt.close('all')
    
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

def make_combine_obs(obj, band=None, badfiles=None):
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
            valid_files.append(obs)

    # print(valid_files)
    return valid_files
    
    
