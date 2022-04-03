# This is a affiliated file for testing and additional runs for python3 library

# This inital design of this library is to use sep to do source finding

import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
import astropy.units as u
import astropy.constants as const


sys.path.append('/home/jchen/work/projects/code_snippets/python3')
from image_tools import FitsImage, source_finder, measure_flux

def combine_detection_py3(dets, units=u.deg, tolerance=0.3*u.arcsec):
    dets_ulist = [dets[0],]
    for idx in range(1, len(dets)):
        det_candidate = dets[idx]
        for udet in dets_ulist:
            coord_difference = udets['ra','dec']*units - det_candidate['ra','dec']*units
            # print(np.sqrt(coord_difference[0]**2 + coord_difference[1]**2))
            if coord_difference[0]**2 + coord_difference[1]**2 < tolerance**2:
                is_unique = False
        if is_unique:
            dets_ulist.append(det_candidate)
    return vstack(dets_ulist)



def run_check_SMGs_py3(basedir, objs=None, bands=['B6','B7'], suffix='combine.ms.auto.cont', 
                   resolutions=['0.3arcsec', '0.6arcsec'], ncol=1,
                   summary_file=None, view_mode='multiple', focus_bands=None,
                   interactive=True, outdir=None, continue_mode=True, central_mask_radius=2,
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
                        #print('Finding source in:', image_fullpath)
                        if not os.path.isfile(image_fullpath):
                            if band in focus_bands:
                                obj_validbands['{}'.format(band)] = False
                            continue
                        savefile = image_name + '.source_found.txt'
                        fitsimage = FitsImage(image_file=image_fullpath)
                        sources_found = source_finder(fitsimage, method='sep', plot=False,
                                                      aperture_scale=3.0, detection_threshold=5.0)
                        nocenter_dets = sources_found[sources_found['radial_distance'] 
                                                       > central_mask_radius]
                        fitsimage.plot(ax=ax_select, show_detections=True, detections=nocenter_dets)
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
                            for ra, dec, flux, snr in sources_found['ra','dec','peak_flux', 'peak_snr']:
                                obj_summary.write('{} {:.6f} {:.6f} '.format(source_idx, ra, dec))
                                obj_summary.write(" {:.4f} {:.2f} ".format(flux, snr))
                                obj_summary.write('\n')
                                source_idx += 1
                fig.savefig(summary_plot, bbox_inches='tight', dpi=100)
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
                            detection_input=int(input(
                                "Dection in Band:{} (integer, 0,1,2,3,-1) [0]?: ".format(band)) or 0)
                            goodfield_input=int(input(
                                "Usable for Band:{} (integer, 0,1,2,-1,-2) [0]?: ".format(band)) or 0)
                            detections[band] = detection_input
                            goodfields[band] = goodfield_input
                        if len(bands) > 0:
                            SMG_input = int(input("Is SMG? (integer, 0,1,2,-1) [0]: ") or 0)
                            RG_input = int(input("Is RG? (integer, 0,1,2,-1) [0]: ") or 0)
                            Jet_input = int(input("Is Jet? (integer, 0,1,2,-1) [0]: ") or 0)
                        plt.close()
                    for band in focus_bands:
                        found_string += ' {} {}'.format(detections[band], goodfields[band])
                    found_string += ' {}'.format(SMG_input)
                    found_string += ' {}'.format(RG_input)
                    found_string += ' {}'.format(Jet_input)
                    with open(summary_file, 'a+') as f:
                        f.write("{}\n".format(found_string)) 
                else:
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
                plt.close()
    except KeyboardInterrupt:
        return 0

def run_measure_flux_py3(basedir, objs=None, bands=['B6','B7'], focus_band=None,
                   suffix='combine.ms.auto.cont', target_wave=None,
                   resolutions=['0.3arcsec', '0.6arcsec'], central_mask_radius=2.0, 
                   selected_resolution='0.3arcsec', ncol=1, aperture_scale=2.0,
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
                        # sources_found = source_finder(image_fullpath,
                                # methods=['single_aperture', 'peak'],
                                # ax=ax_select, pbcor=True, central_mask_radius=2.0, **kwargs)
                        print("images:", image_fullpath)
                        fitsimage = FitsImage(image_file=image_fullpath, pbcor_file=None)
                        sources_found = source_finder(fitsimage, method='sep', plot=False,
                                                      aperture_scale=aperture_scale, 
                                                      detection_threshold=5.0)
                        nocenter_dets = sources_found[sources_found['radial_distance'] 
                                                       > central_mask_radius]
                        fitsimage.plot(ax=ax_select, show_detections=True, detections=nocenter_dets, 
                                       fontsize=0)
                        ax_select.set_title('{}'.format(res))
                        ny, nx = fitsimage.imagesize
                        x_index = (np.arange(0, nx) - nx/2.0) * fitsimage.pixel2deg_ra * 3600 # to arcsec
                        y_index = (np.arange(0, ny) - ny/2.0) * fitsimage.pixel2deg_dec * 3600 # to arcsec
                        ax_select.text(0.5*np.max(x_index), 0.8*np.min(y_index), 
                                "{:.0f}uJy/beam".format(fitsimage.std*1e6), color='white',
                                horizontalalignment='center', verticalalignment='top')
                        if band == focus_band:
                            obj_sourcefound['{}'.format(res)] = sources_found
                        for idx, det in enumerate(nocenter_dets):
                            xdet = (det['x']-fitsimage.imagesize[1]/2.0)*fitsimage.pixel2deg_ra*3600
                            ydet = (det['y']-fitsimage.imagesize[0]/2.0)*fitsimage.pixel2deg_dec*3600
                            ax_select.text(xdet, ydet-1.0, "{}:{:.2f}mJy".format(idx,det['flux']*1e3), 
                                    color='white',horizontalalignment='center', verticalalignment='top', 
                                    fontsize=8, alpha=0.8)
                            # print(idx, xdet, ydet, det['flux'])

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
                                dets_list.append(sources_found[int(idx)])
                        dets_unique = combine_detection_py3(dets_list)
                        selected_image = "{}_{}_{}.{}.image.fits".format(obj, focus_band, suffix, selected_resolution)
                        selected_image_fullpath = os.path.join(obj_focus_band_dir, selected_image)
                        selected_image_pbcor_fullpath = os.path.join(obj_focus_band_dir, 
                                                        selected_image.replace('image', 'pbcor.image'))
                        fitsimage_select = FitsImage(image_file=selected_image_fullpath, 
                                                     pbcor_file=selected_image_pbcor_fullpath)
                        print("Using {} for flux measurements.\n".format(selected_image_fullpath))
                        # sources_flux, sources_flux_snr, sources_radial_distance = flux_measure(
                                # seleted_image_fullpath, coords_unique, target_wave=target_wave,
                                # methods=['adaptive_aperture', 'gaussian', 'peak'], 
                                # calculate_radial_distance=True)
                        sources_flux_aper, sources_flux_err_aper = measure_flux(fitsimage_select, 
                                detections=dets_unique, method='single_aperture', 
                                aperture_scale=aperture_scale)
                        sources_flux_gaussian, sources_flux_err_gaussian = measure_flux(fitsimage_select, 
                                detections=dets_unique, method='gaussian', aperture_scale=aperture_scale)
                        sources_flux_aper_snr = sources_flux_err_aper / sources_flux_err_aper
                        sources_flux_gaussian_snr = sources_flux_err_gaussian / sources_flux_err_gaussian
                        
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
                        print("Totoal {} sources".format(n_uniq))
                        for i in range(n_uniq):
                            detections_summary.write('{} {} {} {:.6f} {:.6f} {} {} {} {} {} {} {}'.format(
                                obj, i, source_type[i], dets_unique[i]['ra'], dets_unique[i]['dec'],
                                sources_flux_aper[i], sources_flux_aper_snr[i],
                                sources_flux_gaussian[i], sources_flux_gaussian_snr[i],
                                sources_flux_peak[i], sources_flux_peak_snr[i],
                                sources_radial_distance[i]))
                            detections_summary.write('\n')
                    detections_summary.close()
                    plt.close()

                else:
                    next_one = int(input("Next one [1/0] [1]: ") or 1)
                    if next_one == 1:
                        if view_mode == 'single':
                            pass
                        elif view_mode == 'multiple':
                            continue
                    else:
                        return 0
                if not summary_file and (view_mode == 'multiple'):
                    is_close = int(input("Close all windows? (1/0) [1]") or 1)
                    if is_close == 0:
                        return 0
                plt.close()
    except KeyboardInterrupt:
        return 0
    print("Failed files:", failed_files)


