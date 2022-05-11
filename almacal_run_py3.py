# This is a affiliated file for testing and additional runs for python3 library

# This inital design of this library is to use sep to do source finding

import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from astropy.table import Table, vstack
import astropy.units as u
import astropy.constants as const
from astropy.stats import sigma_clipped_stats
from scipy import ndimage, interpolate


sys.path.append('/home/jchen/work/projects/code_snippets/python')
from image_tools import FitsImage, source_finder, measure_flux

ROOT_DIR = os.path.join(os.path.expanduser('~'), 'work/projects/almacal/number_counts')

class CalImage(FitsImage):
    def __init__(self, image_file, **kwargs):
        FitsImage.__init__(self, image_file, **kwargs)
        self.lam = (const.c/self.reffreq).decompose().to(u.um)
        # self.fov = (1.02 * (self.lam / (12*u.m)).decompose()* 206264.806).value #in arcsec
        self.fov = (1.02 * (self.restlam / (12*u.m)).decompose()* 206264.806).value #in arcsec
        self.imstat()
    def make_rmap(self):
        ny, nx = self.imagesize
        pixel_rescale = self.pixel2deg_ra * 3600
        x_index = (np.arange(0, nx) - nx/2.0) * pixel_rescale # to arcsec
        y_index = (np.arange(0, ny) - ny/2.0) * pixel_rescale # to arcsec
        x_map, y_map = np.meshgrid(x_index, y_index)
        rmap = np.sqrt(x_map**2 + y_map**2)
        return rmap
    def mask_center(self, radius=None, radius_scale=3.0):
        rmap = self.make_rmap()
        if radius is not None:
            self.central_mask = rmap < radius
        elif radius_scale is not None:
            # mask central region in units of major beam size
            central_mask = rmap < radius_scale*self.bmaj*3600
        return central_mask
    def imstat(self):
        mask = self.imagemask | self.find_structure() | self.mask_center()
        self.mean, self.median, self.std = sigma_clipped_stats(self.image, sigma=5.0, mask=mask)

       
def combine_detection_py3(dets, units=u.deg, tolerance=0.3*u.arcsec):
    dets_ulist = [dets[0],]
    for idx in range(1, len(dets)):
        det_candidate = dets[idx]
        is_unique = True
        for udet in dets_ulist:
            ra_diff = udet['ra'] - det_candidate['ra']
            dec_diff = udet['dec'] - det_candidate['dec']
            # coord_difference = udet['ra','dec'] - det_candidate['ra','dec']
            # print(np.sqrt(coord_difference[0]**2 + coord_difference[1]**2))
            if (ra_diff**2 + dec_diff**2)*units**2 < tolerance**2:
                is_unique = False
        if is_unique:
            dets_ulist.append(det_candidate)
    return vstack(dets_ulist)

def calculate_effectivearea_py3(flux=np.linspace(0.1, 1, 10), snr_threshold=5.0, 
        calimage=None, image_file=None, image_pbcor_file=None, fov_scale=2.0,
        central_mask_scale=1.0):
    """calculate the effective area for given images

    rlimit: in arcsec
    """
    if calimage is None:
        calimage = CalImage(image_file, pbcor_file=image_pbcor_file)
    pix2area = calimage.pixel2deg_ra*calimage.pixel2deg_dec*3600**2  # pixel to arcsec^2
    lam = (const.c/(calimage.reffreq)).decompose().to(u.um)
    fov = 1.02 * (lam / (12*u.m)).decompose()* 206264.806
    mask = calimage.imagemask
    rmap = calimage.make_rmap()
    pbcor = calimage.image / calimage.image_pbcor

    pbcor[mask] = 0.0
    if fov_scale:
        rlimit = 0.5 * fov_scale * fov.value
        pbcor[rmap > rlimit] = 0.0
    if central_mask_scale > 1e-8:
        pbcor[rmap < central_mask_scale*calimage.bmaj*3600] = 0.0
    effarea_list = []
    for f in flux:
        snr = f / (calimage.std * 1000) # from Jy to mJy
        snr_map = snr * pbcor
        #print('selected number of pixels', np.sum(snr_map > snr_threshold))
        area = np.sum(snr_map > snr_threshold) * pix2area
        effarea_list.append(area)
    effarea = np.array(effarea_list) / 3600. # convert to arcmin^2
 
    return np.array([flux, effarea])


def run_check_SMGs_py3(basedir, objs=None, bands=['B6','B7'], suffix='combine.ms.auto.cont', 
                   resolutions=['0.3arcsec', '0.6arcsec'], ncol=1,
                   save_obj_summary=False, aperture_scale=6.0, show_detections=True,
                   summary_file=None, view_mode='multiple', focus_bands=None,
                   interactive=True, outdir=None, continue_mode=True, 
                   central_mask_scale=2,
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
                # for band in focus_bands:
                    # for res in resolutions:
                        # f.write(' ')
                        # f.write(band+'_'+res)
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
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
    try:
        # make a summary plot
        nrow = int(1.0 * len(bands) / ncol)
        fig = plt.figure(figsize=(ncol*4.2*len(resolutions),4*nrow))
        for obj in objs:
            if obj in objs_finished:
                print("{} already done!".format(obj))
                continue
            print(obj)
            if obj_match.match(obj):
                print('>>>>> {}'.format(obj))
                # write into summary file
                if outdir is not None:
                    # obj_outdir = os.path.join(outdir, obj)
                    # if not os.path.isdir(obj_outdir):
                        # os.system('mkdir -p {}'.format(obj_outdir))
                    # open the summary file
                    if save_obj_summary:
                        obj_summary_file = os.path.join(outdir, '{}.sources_found.txt'.format(obj))
                        obj_summary = open(obj_summary_file, 'w+')
                    summary_plot = os.path.join(outdir, '{}.summary.png'.format(obj))
                else:
                    obj_summary = None
                    # obj_outdir = None
                    summary_plot = None
                obj_sourcefound = {}
                obj_validbands = {}
                for band in focus_bands:
                    obj_validbands[band] = True
                    for res in resolutions:
                        obj_sourcefound[band+'_'+res] = []
                #for img in imgs:
                # sources_number = {}
                # obj_sourcefound
                if os.path.isfile(summary_plot):
                    ax = fig.add_subplot(111)
                    ax.axis('off')
                    print("Using existing summary plot...")
                    imagedata = plt.imread(summary_plot)
                    ax.imshow(imagedata, interpolation='none')
                    for band in bands:
                        obj_band_dir = os.path.join(basedir, band, obj)
                        for res in resolutions:
                            if res == '':
                                res_string = ''
                            else:
                                res_string = res+'.'
                            image_name = "{}_{}_{}.{}image.fits".format(obj, band, suffix, res_string)
                            image_fullpath = os.path.join(obj_band_dir, image_name)
                            if not os.path.isfile(image_fullpath):
                                if band in focus_bands:
                                    obj_validbands['{}'.format(band)] = False
                else:
                    ax = fig.subplots(nrow, len(resolutions)*ncol) 
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
                            #print('Finding source in:', image_fullpath)
                            if not os.path.isfile(image_fullpath):
                                if band in focus_bands:
                                    obj_validbands['{}'.format(band)] = False
                                continue
                            savefile = image_name + '.source_found.txt'
                            fitsimage = CalImage(image_file=image_fullpath)
                            sources_found = source_finder(fitsimage, method='sep', plot=False,
                                                          aperture_scale=6.0, detection_threshold=5.0)
                            nocenter_dets = sources_found[sources_found['radial_distance'] 
                                                           > central_mask_scale*fitsimage.bmaj*3600]
                            fitsimage.plot(ax=ax_select, show_detections=show_detections, detections=nocenter_dets, 
                                           show_rms=True, aperture_scale=aperture_scale, fontsize=10)
                            ax_select.set_title('{} {}'.format(band, res))
                            #try:
                            #    sources_found = source_finder(image_fullpath, outdir=obj_outdir, 
                            #            ax=ax[i,j], pbcor=True)
                            #except:
                            #    print("Error found for {}".format(image_name))
                            #    failed_files.append(image_name)
                            if len(sources_found) > 0 and save_obj_summary:
                                obj_summary.write('# {} {} {}\n'.format(obj, band, res))
                                # obj_sourcefound['{}_{}'.format(band, res)] = sources_found
                                source_idx = 0
                                for ra, dec, flux, snr in sources_found['ra','dec','peak_flux', 'peak_snr']:
                                    obj_summary.write('{} {:.6f} {:.6f} '.format(source_idx, ra, dec))
                                    obj_summary.write(" {:.4f} {:.2f} ".format(flux, snr))
                                    obj_summary.write('\n')
                                    source_idx += 1
                    fig.subplots_adjust(wspace=0.2, hspace=0.2)
                    # fig.suptitle(obj)
                    fig.savefig(summary_plot, bbox_inches='tight', dpi=400)
                # write into files
                if save_obj_summary:
                    obj_summary.close()
                found_string = obj
                # for band in focus_bands:
                    # for res in resolutions:
                        # print(obj_sourcefound[band+'_'+res])
                        # found_string += ' '+str(len(obj_sourcefound[band+'_'+res]))
                # print(found_string)
                # save figure
                if interactive:
                    if summary_file: 
                        SMG_input = 0
                        RG_input = 0
                        Jet_input = 0
                        detections = {}
                        goodfields = {}
                        for band in bands:
                            goodfields[band] = 0
                            detections[band] = 0
                        plt.show()
                        print("Single Band:\n") 
                        print("Detection: 0)None +n)N Detections -1)Not Sure")
                        print("Usable: 0)No 1)Full -n)exclude central n*FWHM region")
                        print("General Classification")
                        print("Is SMG: 0)No 1)Yes +n)Multiple -1)Not Sure")
                        print("Is RG: 0)No 1)Yes +n)Multiple -1)Not Sure")
                        print("Is Jet: 0)No 1)Compact +n)Enlongated [in arcsec] -1)Not Sure")
                        have_obs = False
                        for band in focus_bands:
                            if not obj_validbands[band]:
                                continue
                            have_obs = True
                            detection_input=int(input(
                                "Dection in Band:{} (integer, 0,1,2,3,-1) [0]?: ".format(band)) or 0)
                            goodfield_input=int(input(
                                "Usable for Band:{} (integer, 0,1,2,-1,-2) [0]?: ".format(band)) or 0)
                            detections[band] = detection_input
                            goodfields[band] = goodfield_input
                        if have_obs:
                            SMG_input = int(input("Is SMG? (integer, 0,1,2,-1) [0]: ") or 0)
                            RG_input = int(input("Is RG? (integer, 0,1,2,-1) [0]: ") or 0)
                            Jet_input = int(input("Is Jet? (integer, 0,1,2,-1) [0]: ") or 0)
                        else:
                            SMG_input = 0; RG_input=0; Jet_input=0;
                        for band in focus_bands:
                            found_string += ' {} {}'.format(detections[band], goodfields[band])
                        found_string += ' {}'.format(SMG_input)
                        found_string += ' {}'.format(RG_input)
                        found_string += ' {}'.format(Jet_input)
                        with open(summary_file, 'a+') as f:
                            f.write("{}\n".format(found_string)) 
                        # plt.close()
                    elif interactive:
                        next_one = int(input("Next one [1/0] [1]: ") or 1)
                        if next_one == 1:
                            if view_mode == 'single':
                                plt.close()
                            elif view_mode == 'multiple':
                                continue
                        else:
                            return 0
                    if (outdir is None) and (view_mode == 'multiple'):
                        is_close = int(input("Close all windows? (1/0) [1]") or 1)
                        if is_close == 0:
                            return 0
                plt.clf()
    except KeyboardInterrupt:
        plt.close()
        return 0

def run_measure_flux_py3(basedir, objs=None, bands=['B6','B7'], focus_band=None,
                   suffix='combine.ms.auto.cont', target_wave=None,
                   resolutions=['0.3arcsec', '0.6arcsec'], central_mask_scale=2.0, 
                   selected_resolution='0.3arcsec', ncol=1, aperture_scale=6.0,
                   summary_file=None, view_mode='multiple', fov_scale=1.5,
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
        # make a summary plot
        nrow = int(1.0 * len(bands) / ncol)
        fig = plt.figure(figsize=(ncol*4.2*len(resolutions),4*nrow)) 
        for obj in objs:
            if obj in objs_finished:
                print("{} already done!".format(obj))
                continue
            print(obj)
            # fig.clf()
            ax = fig.subplots(nrow, len(resolutions)*ncol) 
            obj_focus_band_dir = os.path.join(basedir, focus_band, obj)
            obj_sourcefound = {}
            if obj_match.match(obj):
                print('>>>>> {}'.format(obj))
                # fig.suptitle(obj)
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
                        # sources_found = source_finder(image_fullpath,
                                # methods=['single_aperture', 'peak'],
                                # ax=ax_select, pbcor=True, central_mask_scale=2.0, **kwargs)
                        print("images:", image_fullpath)
                        fitsimage = CalImage(image_file=image_fullpath, pbcor_file=None)
                        sources_found = source_finder(fitsimage, method='sep', plot=False,
                                                      aperture_scale=aperture_scale, 
                                                      detection_threshold=5.0)

                        valid_selection = ((sources_found['radial_distance'] < 0.5*fov_scale*fitsimage.fov)
                            & (sources_found['radial_distance']>central_mask_scale*fitsimage.bmaj*3600))
                        valid_dets = sources_found[valid_selection]
                        sources_found = valid_dets
                        fitsimage.plot(ax=ax_select, show_detections=True, detections=sources_found, 
                                       fontsize=0, aperture_scale=aperture_scale)
                        # show fov:                                                   
                        ellipse_fov = patches.Ellipse((0, 0), width=fov_scale*fitsimage.fov, 
                                                      height=fov_scale*fitsimage.fov,                  
                                                      angle=0, fill=False, facecolor=None, edgecolor='gray',
                                                      alpha=0.8)                                             
                        ax_select.add_patch(ellipse_fov)
                        ax_select.set_title('{}'.format(res))
                        ny, nx = fitsimage.imagesize
                        x_index = (np.arange(0, nx) - nx/2.0) * fitsimage.pixel2deg_ra * 3600 # to arcsec
                        y_index = (np.arange(0, ny) - ny/2.0) * fitsimage.pixel2deg_dec * 3600 # to arcsec
                        ax_select.text(0.5*np.max(x_index), 0.8*np.min(y_index), 
                                "{:.0f}uJy/beam".format(fitsimage.std*1e6), color='white',
                                horizontalalignment='center', verticalalignment='top')
                        if band == focus_band:
                            obj_sourcefound['{}'.format(res)] = sources_found
                        for idx, det in enumerate(sources_found):
                            xdet = (det['x']-fitsimage.imagesize[1]/2.0)*fitsimage.pixel2deg_ra*3600
                            ydet = (det['y']-fitsimage.imagesize[0]/2.0)*fitsimage.pixel2deg_dec*3600
                            ax_select.text(xdet, ydet-1.0, "{}:{:.2f}mJy".format(idx,det['flux']*1e3), 
                                    color='white',horizontalalignment='center', verticalalignment='top', 
                                    fontsize=8, alpha=0.8)
                            # print(idx, xdet, ydet, det['flux'])
                        plt.show()

                # interactive selecting sources
                if summary_file:
                    detections_summary = open(summary_file, 'a+')
                    dets_list = []
                    is_detection = int(input("Dection in band: {} (integer, 0/1) [1]?: ".format(focus_band)) or 1)
                    # is_detection = 1
                    if is_detection > 0:
                        for res in resolutions:
                            detection_idx = str(
                                    input("The index for the true detections of band {} in reselection:{}\nSeperate with comma: ".format(
                                               focus_band, res)))
                            for idx in detection_idx.split(','):
                                if idx == '':
                                    continue
                                # ra, dec, flux, snr = obj_sourcefound['{}'.format(res)][int(idx)]
                                dets_list.append(obj_sourcefound[res][int(idx)])
                        dets_unique = combine_detection_py3(dets_list, units=u.deg)
                        selected_image = "{}_{}_{}.{}.image.fits".format(obj, focus_band, suffix, selected_resolution)
                        selected_image_fullpath = os.path.join(obj_focus_band_dir, selected_image)
                        selected_image_pbcor_fullpath = os.path.join(obj_focus_band_dir, 
                                                        selected_image.replace('image', 'pbcor.image'))
                        fitsimage_select = CalImage(image_file=selected_image_fullpath, 
                                                     pbcor_file=selected_image_pbcor_fullpath)
                        print("Using {} for flux measurements.\n".format(selected_image_fullpath))
                        # sources_flux, sources_flux_snr, sources_radial_distance = flux_measure(
                                # seleted_image_fullpath, coords_unique, target_wave=target_wave,
                                # methods=['adaptive_aperture', 'gaussian', 'peak'], 
                                # calculate_radial_distance=True)
                        sources_flux_aper, sources_flux_err_aper = measure_flux(fitsimage_select, 
                                detections=dets_unique, method='single-aperture', 
                                minimal_aperture_size=fitsimage_select.bmaj*2*3600,
                                aperture_scale=aperture_scale)
                        sources_flux_gaussian, sources_flux_err_gaussian = measure_flux(fitsimage_select, 
                                detections=dets_unique, method='gaussian', aperture_scale=aperture_scale)
                        sources_flux_aper_snr = sources_flux_aper / sources_flux_err_aper
                        sources_flux_gaussian_snr = sources_flux_gaussian / sources_flux_err_gaussian
                        
                        sources_radial_distance = dets_unique['radial_distance']
                        sources_flux_peak = dets_unique['peak_flux']
                        sources_flux_peak_snr = dets_unique['peak_snr']
                        # sources_flux, sources_flux_snr = sources_flux, sources_flux_snr
                        # print("sources_flux", sources_flux)
                        # print("sources_flux_snr", sources_flux_snr)
                        #print(sources_flux)
                        
                        # classify the detections: 0)not sure 1)SMG; 2)radio;
                        n_uniq = len(dets_unique)
                        source_type = []
                        print("Classify the detections: 0)not sure 1)SMG; 2)radio:")
                        for si in range(n_uniq):
                            print("For selected source:")
                            print("     {}: flux={:.2f} snr={:.2f} dist={:.2f}".format(
                                si, sources_flux_peak[si]*1000., sources_flux_peak_snr[si], 
                                sources_radial_distance[si]))
                        for si in range(n_uniq):
                            s_type = int(input("Selected source: {} type:(integer, 0/1/2) [1]?: ".format(si)) or 1)
                            source_type.append(s_type)
                            print('measureed flux:', sources_flux_aper[si]*1e3, sources_flux_err_aper[si]*1e3)
                        print("Totoal {} sources".format(n_uniq))
                        for i in range(n_uniq):
                            detections_summary.write('{} {} {} {:.6f} {:.6f} {} {} {} {} {} {} {}'.format(
                                obj, i, source_type[i], dets_unique[i]['ra'], dets_unique[i]['dec'],
                                sources_flux_aper[i]*1000, sources_flux_aper_snr[i],
                                sources_flux_gaussian[i]*1000, sources_flux_gaussian_snr[i],
                                sources_flux_peak[i]*1000, sources_flux_peak_snr[i],
                                sources_radial_distance[i]))
                            detections_summary.write('\n')
                        detections_summary.close()
                    else:
                        print("No detection is confirmed. Next one...")
                plt.clf()
    except KeyboardInterrupt:
        plt.close()
        print("Failed files:", failed_files)
        return 0
    print("Failed files:", failed_files)

def run_calculate_effarea_py3(imagedir=None, flux=np.linspace(0.01, 10, 500),  objs=None, band=None, 
        suffix='combine.ms.auto.cont', resolution='0.3arcsec', objs_nocenter=None, 
        almacal_catalogue=None, overwrite=True, fov_scale=2.0,
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
            mask_scale = 0.0
            if almacal_catalogue is not None:
                mask_code = almacal_cat[almacal_cat['obj'] == obj]['goodfield_{}'.format(band)].data[0]
                if mask_code < 0:
                    mask_scale= -1.0 * mask_code
            elif objs_nocenter is not None:
                mask_scale = 2.0
            _, objarea = calculate_effectivearea_py3(flux, snr_threshold=snr_threshold, 
                        image_file=image_fullpath, image_pbcor_file=image_pbcorr_fullpath,
                        central_mask_scale=mask_scale, fov_scale=fov_scale, **kwargs)
            effarea[band] += objarea
    if savefile:
        effarea.write(savefile, format='ascii', overwrite=overwrite)

def calculate_simu(snr, simu_name, simu_boosting_file=None, flux_mode='aperture', 
        simu_completeness_file=None, simu_fake_file=None):
    """calculate simulations
    """
    if simu_boosting_file is None:
        simu_boosting_file = simu_name + '_boosting.dat'
    if simu_completeness_file is None:
        simu_completeness_file = simu_name + '_completeness.dat'
    if simu_fake_file is None:
        simu_fake_file = simu_name + '_fake.dat'
    simu_boosting = Table.read(simu_boosting_file, format='ascii')
    simu_flux_boosting = simu_boosting['flux_{}'.format(flux_mode)]/simu_boosting['flux_input']
    simu_completeness = Table.read(simu_completeness_file, format='ascii')
    simu_fake = Table.read(simu_fake_file, format='ascii')

    snr_mid = 0.5*(snr[1:] + snr[:-1])
    n_bins = len(snr)

    flux_boosting = np.zeros((n_bins-1,3))
    comp_array = np.zeros(n_bins-1)
    fake_array = np.zeros(n_bins-1)

    for i in range(n_bins-1):
        simu_boosting_select = (simu_boosting['snr_peak'] > snr[i]) & (simu_boosting['snr_peak'] < snr[i+1])
        simu_comp_select = (simu_completeness['snr_peak'] > snr[i]) & (simu_completeness['snr_peak'] < snr[i+1])
        simu_fake_select = (simu_fake['snr_peak'] > snr[i]) & (simu_fake['snr_peak'] < snr[i+1])
        if np.sum(simu_boosting_select) < 1:
            flux_boosting[i] = np.nan,np.nan, np.nan
            continue
        if np.sum(simu_comp_select) < 1:
            comp_array[i] = np.nan
            continue
        if np.sum(simu_fake_select) < 1:
            fake_array[i] = np.nan
            continue
        #print(flux_boosting[i])
        #print(simu_flux_boosting[simu_boosting_select])
        #print(np.percentile(simu_flux_boosting[simu_boosting_select], 15.8))
        #print(np.percentile(simu_flux_boosting[simu_boosting_select], 84.2))
        flux_boosting[i] = (1.0*np.median(simu_flux_boosting[simu_boosting_select]), np.percentile(simu_flux_boosting[simu_boosting_select], 15.8), 
                                                                               np.percentile(simu_flux_boosting[simu_boosting_select], 84.2))
        comp_array[i] = 1.0*np.sum(simu_completeness[simu_comp_select]['is_recovered'])/np.sum(simu_comp_select)
        fake_array[i] = 1.0*np.sum(simu_fake[simu_fake_select]['is_fake'])/np.sum(simu_fake_select)
    return snr_mid, flux_boosting, comp_array, fake_array

def run_number_counts_p3(flist=None, detections_file=None, effective_area_file=None, band='B6',
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
            simu_file = default_simulation
        else:
            simu_file = os.path.join(simulation_folder, obj, obj+'_simulation')
        if not os.path.isfile(simu_file):
            # raise ValueError('No simulation could be found for {}'.format(obj))
            print('Warning: using the default simulation results!')
            simu_file = default_simulation
        try:
            snr_list, flux_boosting, comp_list, fake_rate_list = calculate_simu(
                np.arange(0.2, 11, 0.2), simu_file)
        except:
            snr_list, flux_boosting, comp_list, fake_rate_list = calculate_simu(
                np.arange(0.2, 11, 0.2), simu_file)
            print("Using default simulations: {}".format(default_simulation))
        cs_comp = interpolate.interp1d(snr_list, comp_list, fill_value='extrapolate')
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
        flux_bootstrap = flux + np.random.randn(len(index_SMG_bootstrap)) * flux_err        
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
    #mode_NN = np.mean(mode_simu[0], axis=0)
    #mode_NN_err = np.std(mode_simu[0], axis=0)
    mode_NN = np.ma.mean(np.ma.masked_invalid(mode_simu[0]),axis=0).data
    mode_NN_err = np.ma.std(np.ma.masked_invalid(mode_simu[0]),axis=0).data
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
        plot_number_counts_py3(flist=flist, NN=mode_NN, NN_err=[NN_err_lower, NN_err_upper], ax=ax, band=band)

def plot_number_counts_py3(datafile=None, flist=None, NN=None, NN_err=None, ax=None, band=None,
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



"""
# the pipeline code

# generate the effective area for all the bands
for band in ['B3','B4','B5','B6','B7','B8']:
    band_good = dets_all['goodfield_{}'.format(band)] != 0
    run_calculate_effarea_py3(
        imagedir='make_all_good_images/{}'.format(band),
        flux=np.logspace(-2, 1.2, 100),
        objs=dets_all['obj'][band_good],
        band=band,
        almacal_catalogue='check_SMGs_sextractor.txt',
        savefile='{0}/{0}_effarea_py3.txt'.format(band))





"""
