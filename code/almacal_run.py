# This file include all the tasks for ALMACAL number counts projects

import glob
import analysisUtils as au
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
from matplotlib import pyplot as plt


def gen_obstime(base_dir=None, output_dir=None, bad_obs=None):
    """generate the on-source time and spw distribution for the whole almacal dataset
    """

    p_obj = re.compile('J\d+[+-]\d')
    p_obs = re.compile('uid___')
    
    bad_obs = '/Users/jchen/Desktop/projects/almacal/data/broken_obs.txt'
    almacal_info_file = '/Users/jchen/Desktop/projects/almacal/data/almacal_timeOnSource.txt'
    all_obs = Table.read(almacal_info_file, format='ascii')
    all_obs.sort(['B6', 'B7'])
    all_obs.reverse()

    if bad_obs is not None:
        with open(bad_obs) as f:
            all_bad_obs = f.readlines()
        for i in range(len(all_bad_obs)):
            all_bad_obs[i] = all_bad_obs[i].strip()

    band_match = re.compile('_(?P<band>B\d{1,2})$')
    obj_match = re.compile('J\d{4}[-+]\d{4}')

    for i,obj in enumerate(all_obs['obj']):
        obj_exptime = {'B3':0, 'B4':0, 'B5':0, 'B6':0, 
                    'B7':0, 'B8':0,  'B9':0, 'B10':0}
        if not obj_match.match(obj):
            print('Error load obj:', obj)
            continue
        print('index=', i, "obj:", obj)
            
        obj_dirname = base_dir +'/'+ obj
        obj_output_dir = output_dir + '/' + obj
        os.system('mkdir {}'.format(obj_output_dir))
        spw_stat(obj_dirname, plot=True, showfig=False, \
                figname=obj_output_dir+'/'+obj+'.pdf', \
                plotbands=['B5','B6','B7','B8'], 
                savedata=True, filename=obj_output_dir+'/'+ obj+'.fits')


def gen_image(obj, band=None, outdir='./', **kwargs):
    """make images for one calibrator on all or specific band
    """
    for obs in os.listdir(obj):
        if band is not None:
            band_match = re.compile('_(?P<band>B\d{1,2})')
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

def show_images(fileglob, mode='auto', nrow=3, ncol=3):
    """show images in an interative ways, and record the input from inspector
    
    Parameters:
        fileglob: the match patter of the image filenames, like: './*.image'
    """
    all_files = glob.glob(fileglob)
    if mode == 'random':
        all_files = np.random.choice(all_files, nrow * ncol)
    for i in range(0, len(all_files), ncol*nrow):
        fig = plt.figure(figsize=(3*nrow,3*ncol))
        for j in range(0, nrow*ncol):
            print(all_files[i+j])
            try:
                with fits.open(all_files[i+j]) as fitsimage:
                    imageheader = fitsimage[0].header
                    imagedata = fitsimage[0].data
                    wcs = WCS(imageheader)
            except:
                print("Error in reading: {}".format(all_files[i+j]))
                continue
            ax = fig.add_subplot(nrow, ncol, j)#, projection=wcs, slices=(0, 0, 'x', 'y'))
            ax.imshow(imagedata[0,0,:,:], origin='lower')#, cmap='viridis')
            ax.set_xlabel('RA')
            ax.set_ylabel('Dec')
    plt.show()

