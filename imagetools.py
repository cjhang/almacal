# A collection of tools related with image analysis

# Author: Jianhang Chen
# Email: cjhastro@gmail.com
#

import numpy as np
from scipy import signal, ndimage 
from scipy.interpolate import CubicSpline
from scipy import interpolate
from astropy import units as u
from astropy.coordinates import SkyCoord

from astropy.io import fits
from astropy.stats import sigma_clipped_stats, sigma_clip, SigmaClip
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from photutils import (DAOStarFinder, EllipticalAperture, aperture_photometry, CircularAperture,
        RectangularAperture, find_peaks, Background2D, MedianBackground, SExtractorBackground)
import matplotlib.pyplot as plt
from matplotlib import patches
from astropy.modeling import models, fitting
from astropy.convolution import Gaussian2DKernel, convolve

def source_finder2(fitsimage, outdir='./', savefile=None, model_background=True, 
                  threshold=5.0, debug=False, algorithm='DAOStarFinder', return_image=False,
                  filter_size=None, box_size=None, methods=['aperture','gaussian','peak'],
                  subtract_background=False, figname=None, ax=None, pbfile=True,
                  fov_scale=2.0, mask_threshold=5., second_check=True, baseimage=None,
                  cmap=None, return_flux=True):
    """finding point source in the image

    This is a two stage source finding algorithm. First, DAOStarFinder or find_peak will be used to find
    the sources. Then, additional total and peak SNR checking is used to confirm the detection. The 
    second stage can be configured by 

    mask_threshold: the mask size of known_sources
    central_mask: central mask radius in units of FWHM of final synthesised beam
    """
    with fits.open(fitsimage) as hdu:
        header = hdu[0].header
        wcs = WCS(header)
        data = hdu[0].data
        ny, nx = data.shape[-2:]
        freq = header['CRVAL3']
        lam = (const.c/(freq*u.Hz)).decompose().to(u.um)
        fov = 1.02 * (lam / (12*u.m)).decompose()* 206264.806
        fcenter = [header['CRVAL1'], header['CRVAL2']] # in degrees
        pixel_center = [nx/2., ny/2.]
    if pbfile:
        fitsimage_path = os.path.dirname(fitsimage)
        fitsimage_basename = os.path.basename(fitsimage)
        fitsimage_pbcor = os.path.join(fitsimage_path, fitsimage_basename.replace('pbcor.fits', 
                                       'pb.fits.gz'))
        if not os.path.isfile(fitsimage_pbcor):
            print("No primary beam corrected data found, the final flux can be wrong!")
            pbcor = False
        # load the pbcorred data
        with fits.open(fitsimage_pbcor) as hdu_pbcor:
            pbcor_map = hdu_pbcor[0].data.reshape(ny,nx)
     
    data_masked = np.ma.masked_invalid(data.reshape(ny, nx) * pbcor_map)
    # for DAOStarFinder, in pixel space
    pixel_scale = 1/np.abs(header['CDELT1'])
    fov_pixel = fov / 3600. * pixel_scale
    fwhm = header['BMAJ']*3600*u.arcsec
    fwhm_pixel = header['BMAJ']*pixel_scale
    bmaj, bmin = header['BMAJ']*3600, header['BMIN']*3600 # convert to arcsec
    a, b = header['BMAJ']*pixel_scale, header['BMIN']*pixel_scale
    ratio = header['BMIN'] / header['BMAJ']
    theta = header['BPA']
    beamsize = np.pi*a*b/(4*np.log(2))

    mean, median, std = sigma_clipped_stats(data_masked, sigma=5.0, iters=5)  
    
    if debug:
        print('image shape', ny, nx)
        print("sigma stat:", mean, median, std)
        print('fwhm_pixel:', fwhm_pixel)
        print('a, b, theta', a, b, theta)
        print('beamsize', beamsize)

    if model_background:
        if filter_size is None:
            filter_size = 3#max(int(fwhm_pixel/1.5), 6)
        if box_size is None:
            box_size = int(fwhm_pixel*2)
        if debug:
            print('filter_size', filter_size)
            print('box_size', box_size)
        sigma_clip = SigmaClip(sigma=3.)
        # bkg_estimator = MedianBackground()
        bkg_estimator = SExtractorBackground() 
        bkg = Background2D(data_masked.data, (box_size, box_size), 
                filter_size=(filter_size, filter_size),
                mask=data_masked.mask,
                sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        background = bkg.background
    else:
        background = median
    data_masked_sub = data_masked - background
    
    # find stars
    if algorithm == 'DAOStarFinder': # DAOStarFinder
        daofind = DAOStarFinder(fwhm=fwhm_pixel, threshold=0.8*threshold*std, ratio=ratio, 
                                theta=theta+90, sigma_radius=1.5, sharphi=1.0, sharplo=0.2,)  
        sources_found_candidates = daofind(data_masked_sub)#, mask=known_mask)
        # if debug:
            # print("candidates:", sources_found_candidates)
        if len(sources_found_candidates) < 1:
            print("No point source!")
            # return 0
            sources_found_candidates = None
        else:
            sources_found_x_candidates = sources_found_candidates['xcentroid'] 
            sources_found_y_candidates = sources_found_candidates['ycentroid']
            # sources_found_candidates_peak = sources_found_candidates['peak']
            # peak_select = sources_found_candidates_peak > 0.5*threshold*std #two stage source finding
            # sources_found_candidates = sources_found_candidates[peak_select]
            # sources_found_x_candidates = sources_found_x_candidates[peak_select]
            # sources_found_y_candidates = sources_found_y_candidates[peak_select]
        
    elif algorithm == 'find_peak': # find_peak
        sources_found_candidates = find_peaks(data_masked_sub, threshold=threshold*std, 
                box_size=fwhm_pixel)
        # if debug:
            # print("candidates:", sources_found_candidates)
        if len(sources_found_candidates) < 1:
            print("No point source!")
            # return 0
            sources_found_candidates = None
        else:
            sources_found_x_candidates = sources_found_candidates['x_peak'].data 
            sources_found_y_candidates = sources_found_candidates['y_peak'].data
            #source_found_peak = sources_found['peak_value']

    else:
        raise ValueError("Unsurport algorithm: {}!".format(algorithm))
   
    # numerical aperture correction
    image_gaussian = make_gaussian_image((2.*np.ceil(a), 2.*np.ceil(a)), fwhm=[b,a], 
                offset=[0.5, 0.5], # 0.5 offset comes from central offset of photutils of version 0.5
                theta=theta/180.*np.pi) 
    flux_g = auto_photometry(image_gaussian, bmaj=b, bmin=a, theta=theta/180.*np.pi, 
                             methods='aperture', debug=False, aperture_correction=1.0,
                             beamsize=1.0)
    aperture_correction = 1/flux_g[0]
    if debug:
        print('aperture correction', aperture_correction)

    # flux measurements
    flux_auto = []
    flux_snr_list = []
    sources_found = False
    if sources_found_candidates:
        # print('source_found_peak',source_found_peak)
        sources_found_center_candidates = list(zip(sources_found_x_candidates, 
                                                   sources_found_y_candidates))
        #sources_found_coords_candidates = pixel_to_skycoord(sources_found_x, sources_found_y, wcs)
        sources_found_x = []
        sources_found_y = []
        # sources_found_coords = []


        # aperture photometry based on source finding coordinates
        ## simple aperture flux
        #aper_found = EllipticalAperture(sources_found_center_candidates, 1*a, 1*b, 
        #                                theta=theta+90/180*np.pi)
        #phot_table_found = aperture_photometry(data_masked, aper_found)
        #flux_aper_found = (phot_table_found['aperture_sum'] / beamsize * 1000).tolist() # convert mJy

        # automatically aperture photometry
        seg_radius = 2.0*np.int(a)
        segments = RectangularAperture(sources_found_center_candidates, seg_radius, 
                seg_radius, theta=0)
        try:
            segments_mask = segments.to_mask(method='center')
        except:
            segments_mask = []
        for i,s in enumerate(segments_mask):
            if subtract_background:
                data_cutout = s.cutout(data_masked_sub)
            else:
                data_cutout = s.cutout(data_masked)
            flux_list = auto_photometry(data_cutout, bmaj=b, bmin=a, beamsize=beamsize,
                                        theta=theta/180*np.pi, debug=False, methods=methods,
                                        aperture_correction=aperture_correction)
            # if pbcor:
                # data_pbcor_cutout = s.cutout(data_pbcor_masked)
                # flux_pbcor_list = auto_photometry(data_pbcor_cutout, bmaj=b, bmin=a, 
                        # beamsize=beamsize, theta=theta/180*np.pi, debug=False, methods=methods)
            #print("flux_list", flux_list)
            is_true = True
            if second_check:
                if 'aperture' in methods:
                    if flux_list[methods.index('aperture')] < threshold*std:
                        is_true = False
                if 'gaussian' in methods:
                    if flux_list[methods.index('gaussian')] < threshold*std:
                        is_true = False
                if 'peak' in methods:
                    if flux_list[methods.index('peak')] < threshold*std:
                        is_true = False
                # checking whether within the designed fov
                if ((sources_found_x_candidates[i]-pixel_center[0])**2 
                    + (sources_found_y_candidates[i]-pixel_center[1])**2) >\
                            (fov_scale*fov_pixel*0.5)**2:
                    is_true = False

            if is_true: 
                sources_found_x.append(sources_found_x_candidates[i])
                sources_found_y.append(sources_found_y_candidates[i])
                flux_snr = np.array(flux_list) / std
                if pbfile:
                    pbcor_pixel = pbcor_map[sources_found_y_candidates[i], 
                                            sources_found_x_candidates[i]]
                    flux_auto.append(np.array(flux_list) * 1000. / pbcor_pixel)

                else:
                    flux_auto.append(np.array(flux_list) * 1000.)
                flux_snr_list.append(flux_snr)

        if len(flux_auto)>0:
            sources_found = True
        sources_found_center = list(zip(sources_found_x, sources_found_y))
        sources_found_coords = pixel_to_skycoord(sources_found_x, sources_found_y, wcs)
        # return segments_mask
    if debug:
        if sources_found:
            print("sources_found_center", sources_found_center)
            print("sources_found_coords", sources_found_coords)
            print('auto_photometry:', flux_auto)


    if debug or figname or ax:
        # visualize the results
        if ax is None:
            fig= plt.figure()
            ax = fig.add_subplot(111)
        ny, nx = data_masked.shape[-2:]
        scale = np.abs(header['CDELT1'])*3600
        x_index = (np.arange(0, nx) - nx/2.0) * scale
        y_index = (np.arange(0, ny) - ny/2.0) * scale
        #x_map, y_map = np.meshgrid(x_index, y_index)
        #ax.pcolormesh(x_map, y_map, data_masked)
        extent = [np.min(x_index), np.max(x_index), np.min(y_index), np.max(y_index)]
        ax.imshow(data_masked, origin='lower', extent=extent, interpolation='none', vmax=10.*std, 
                  vmin=-3.0*std, cmap=cmap)
        
        ax.text(0, 0, '+', color='r', fontsize=24, fontweight=100, horizontalalignment='center',
                verticalalignment='center')
        ellipse = patches.Ellipse((0.8*np.min(x_index), 0.8*np.min(y_index)), width=bmin, height=bmaj, 
                                  angle=theta, facecolor='orange', edgecolor=None, alpha=0.8)
        ax.add_patch(ellipse)
        if debug:
            ellipse = patches.Ellipse((0, 0), width=fov_scale*fov, height=fov_scale*fov, 
                                      angle=0, fill=False, facecolor=None, edgecolor='grey', 
                                      alpha=0.8)
            ax.add_patch(ellipse)
        ax.text(0.7*np.max(x_index), 0.9*np.min(y_index), 'std: {:.2f}mJy'.format(std*1000.), 
                color='white', fontsize=10, horizontalalignment='center', verticalalignment='center')
        
        if sources_found:
            for i,e in enumerate(sources_found_center):
                yy = scale*(e[0]-ny/2.)
                xx = scale*(e[1]-ny/2.)
                ellipse = patches.Ellipse((yy, xx), width=3*b*scale, 
                                          height=3*a*scale, angle=theta, facecolor=None, fill=False, 
                                          edgecolor='white', alpha=0.8, linewidth=2)
                ax.add_patch(ellipse)
                ax.text(1.0*yy+1, 1.0*xx, 
                        "[{}]{:.2f}mJy".format(i, flux_auto[i][0]), 
                        color='white', fontsize=10, horizontalalignment='left', 
                        verticalalignment='center')
                        #bbox=dict(boxstyle="round", ec=(1,0.5,0.5), fc=(1,0.8,0.8), alpha=0.3))
                # plt.figure()
                # plt.imshow(data_masked_sub)
                # ellipse_aper = EllipticalAperture(e, 2.*a, 2.*a, theta=0)
                # ap_patch = ellipse_aper.plot(color='white', lw=2)


        if figname:
            ax.set_title(figname)
            fig.savefig(os.path.join(outdir, figname+'.png'), dpi=200)
        if debug:
            plt.show()


    if savefile and sources_found:
        savefile_fullpath = os.path.join(outdir, savefile)
        with open(savefile_fullpath, 'w+') as sfile:
            sfile.write('# ra[arcsec]  dec[arcsec] ')
            for m in methods:
                sfile.write(' '+m+'[mJy] ')
            sfile.write('\n')
            for ra, dec, flux in zip(sources_found_coords.ra, sources_found_coords.dec, flux_auto):
                sfile.write('{:.6f}  {:.6f} '.format(ra.to(u.arcsec).value,
                    dec.to(u.arcsec).value))
                for f_m in flux:
                    sfile.write(" {:.8f} ".format(f_m))
                sfile.write('\n')
    if len(flux_auto)>0:
        if return_flux:
            return list(zip(sources_found_coords.ra.to(u.arcsec).value, 
                            sources_found_coords.dec.to(u.arcsec).value, flux_auto, flux_snr_list))
        else:
            return list(zip(sources_found_coords.ra.to(u.arcsec).value, 
                            sources_found_coords.dec.to(u.arcsec).value))

    else:
        return []
 
def auto_photometry2(image, bmaj=1, bmin=1, theta=0, beamsize=None, debug=False, 
                    methods=['adaptive_aperture','gaussian', 'peak'], aperture_correction=1.068,
                    rms=None):
    """autoumatically measure the flux with different methods
    
    aperture_size = 1.0 # in units of beam size
    """
    ysize, xsize = image.shape
    y_center, x_center = ysize/2., xsize/2.
    flux_list = []
    flux_error_list = []
    if beamsize is None:
        print("Warning: Unit the sum and maximum, without dividing the beam size.")
    sigma2FWHM = 2.35482
    # Aperture Photometry
    if 'aperture' in methods:
        aperture_size=1.0
        aperture = EllipticalAperture((x_center, y_center), aperture_size*bmaj, aperture_size*bmin, theta)
        # extract_area = lambda x: x.area()
        # area_apers = np.array(list(map(extract_area, apertures)))
        
        phot_table = aperture_photometry(image, aperture)
        flux_in_aper = phot_table['aperture_sum'].data[0] * aperture_correction
        # extract_qtable = lambda x: x.data[0]
        # flux_apers = np.array(list(map(extract_qtable, phot_table.columns[3:].values())))

        if debug:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            im = ax.imshow(image, origin='lower', interpolation='none')
            plt.colorbar(im, fraction=0.046, pad=0.04)
            aperture.plot(color='gray', lw=2)
        if beamsize:
            flux_list.append(1.0*flux_in_aper/beamsize)
            flux_error_list.append(rms * np.sqrt(aperture.area())/beamsize)
        else:
            flux_list.append(flux_in_aper)
            flux_error_list.append(rms * np.sqrt(aperture.area()))
        # print('aperture.area;',aperture.area())

    if 'adaptive_aperture' in methods:
        radii = np.arange(0.2*bmaj/2, ysize/2./np.cos(theta), 0.5/np.cos(theta))

        apertures = [EllipticalAperture((x_center, y_center), a, a*bmin/bmaj, theta) for a in radii]
        extract_area = lambda x: x.area()
        area_apers = np.array(list(map(extract_area, apertures)))
        
        phot_table = aperture_photometry(image, apertures)
        extract_qtable = lambda x: x.data[0]
        flux_apers = np.array(list(map(extract_qtable, phot_table.columns[3:].values())))
        snr_apers = flux_apers / (rms*np.sqrt(area_apers))
        print("flux_apers:", 1000*flux_apers/beamsize)
        print("snr_apers:", snr_apers)

        slopes = np.ones_like(snr_apers)
        #slopes[1:] = np.diff(flux_apers)/np.max(flux_apers)/np.diff(area_apers) #Normalized slope
        slopes[1:] = np.diff(snr_apers)/np.max(snr_apers)/np.diff(area_apers) #Normalized slope
        print('slopes', slopes)
        slope_selection = slopes < 1e-8
        if debug:
            print('position:{}/{}'.format(np.min(np.where(slope_selection)[0]), len(radii)))
            # print('slopes', slopes)
        if len(flux_apers[slope_selection])<3:
            flux_apers_stable = flux_apers[-2]
            area_apers_stable = area_apers[-2]
        else:
            flux_apers_stable = flux_apers[slope_selection][0]
            area_apers_stable = area_apers[slope_selection][0]
        slope_selection_index = np.where(slope_selection == 1)[0][0]
        if debug:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            im = ax.imshow(image, origin='lower', interpolation='none')
            plt.colorbar(im, fraction=0.046, pad=0.04)
            for ap in apertures:
                im = ap.plot(color='gray', lw=2)
            im = apertures[slope_selection_index].plot(color='black', lw=4)

        if beamsize:
            flux_apers_stable = 1.0*flux_apers_stable/beamsize
            flux_error_stable = rms*np.sqrt(area_apers_stable)/beamsize
        else:
            flux_apers_stable = flux_apers_stable
            flux_error_stable = rms*np.sqrt(area_apers_stable)
        print('stable:', flux_apers_stable, flux_error_stable)
        flux_list.append(flux_apers_stable)
        flux_error_list.append(flux_error_stable)

    if 'gaussian' in methods:
        yidx, xidx = np.indices((ysize, xsize))
        yrad, xrad = yidx-ysize/2., xidx-xsize/2.

        image_scale = np.max(image)
        image_norm = image / image_scale
        p_init = models.Gaussian2D(amplitude=1, x_stddev=bmaj/sigma2FWHM, y_stddev=bmin/sigma2FWHM, 
                theta=theta, bounds={"x_mean":(-2.,2.), "y_mean":(-2.,2.), 
                                     "x_stddev":(0., xsize), "y_stddev":(0., ysize)})
        fit_p = fitting.LevMarLSQFitter()
        p = fit_p(p_init, xrad, yrad, image_norm, 
                  weights=1/(yrad**2+xrad**2+(0.25*bmaj)**2+(0.25*bmin)**2))
        flux_fitted = 2*image_scale*np.pi*p.x_stddev.value*p.y_stddev.value*p.amplitude.value
        # calculate the area with amplitude < rms
        gaussian_area = np.sum(p(xrad, yrad) > rms/image_scale)
        if debug:
            print("Initial guess:", p_init)
            print("Fitting:", p)
            print("Flux:", flux_fitted)
            fig, ax = plt.subplots(1, 3, figsize=(8, 3.5))
            im0 = ax[0].imshow(image_norm, origin='lower', interpolation='none')
            plt.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)
            gaussian_init = patches.Ellipse((ysize/2., xsize/2.), height=p_init.y_stddev.value, 
                    width=p_init.x_stddev.value, angle=p_init.theta.value/np.pi*180,
                    linewidth=2, facecolor=None, fill=None, edgecolor='orange', alpha=0.8)
            ax[0].add_patch(gaussian_init)
            ax[0].set_title("Data")
            im1 = ax[1].imshow(p(xrad, yrad), origin='lower', interpolation='none')
            plt.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)
            ax[1].set_title("Model")
            im2 = ax[2].imshow(image_norm - p(xrad, yrad), origin='lower', interpolation='none')
            plt.colorbar(im2, ax=ax[2], fraction=0.046, pad=0.04)
            ax[2].set_title("Residual")
        if beamsize:
            flux_list.append(1.0*flux_fitted/beamsize)
            flux_error_list.append(rms*np.sqrt(gaussian_area)/beamsize)
        else:
            flux_list.append(flux_fitted)
            flux_error_list.append(rms*np.sqrt(gaussian_area))
    if 'peak' in methods:
        # the peak value has the same units of one pixel
        # if beamsize:
            # flux_list.append(np.max(image))#)#/beamsize)
        # else:
        flux_list.append(np.max(image))
        flux_error_list.append(rms)
    return flux_list, flux_error_list


def flux_measure2(image, coords_list, methods=['aperture', 'gaussian','peak'], pbcor=True, 
                 model_background=False, jackknif=True, target_wave=None,
                 subtract_background=False, debug=False, ax=None, fov_scale=2.0,
                 calculate_radial_distance=False):
    """measure the flux from the images by given coordinates
    the images can be multiple images on the same field but with different resolutions
    the flux from the coarser resolution will be preferred
    It also offers the validatiy check bettween the sythesized beam and the PSF
    """
    # convert the coordinates list to array
    coords_array = np.array(coords_list)
    known_sources_coords = SkyCoord(ra=coords_array[:,0]*u.arcsec, dec=coords_array[:,1]*u.arcsec)

    # read the images
    with fits.open(image) as hdu:
        header = hdu[0].header
        wcs = WCS(header)
        data = hdu[0].data
        ny, nx = data.shape[-2:]
        deg2pixel = 1/np.abs(header['CDELT1'])
        freq = header['CRVAL3']
        lam = (const.c/(freq*u.Hz)).decompose().to(u.um)
        fov = 1.02 * (lam / (12*u.m)).decompose()* 206264.806
        fov_pixel = fov / 3600. * deg2pixel
        a, b = header['BMAJ']*deg2pixel, header['BMIN']*deg2pixel
        theta = header['BPA']
        bmaj, bmin = header['BMAJ']*3600, header['BMIN']*3600 # convert to arcsec
        beamsize = np.pi*a*b/(4*np.log(2))
        data_masked = np.ma.masked_invalid(data.reshape(ny, nx))
        known_mask = data_masked.mask
        # mask the known data to calculate the sensitivity of the map
        known_sources_pixel = skycoord_to_pixel(known_sources_coords, wcs)
        known_sources_center = list(zip(known_sources_pixel[0], known_sources_pixel[1]))
        for p in known_sources_center:
            aper = CircularAperture(p, 2.0*a)
            aper_mask = aper.to_mask(method='center')
            known_mask = np.bitwise_or(known_mask, aper_mask[0].to_image((ny,nx)).astype(bool))
        data_field = np.ma.array(data_masked, mask=known_mask) 
        # mean, median, std = sigma_clipped_stats(data_field, sigma=5.0, iters=5)  
    if jackknif:
        image_jeckknifed = image.replace('image', 'jackknif_image')
        if not os.path.isfile(image_jeckknifed):
            print("Warning: no jackknifed image found {}, use the origin image instead.".format(
                   image_jeckknifed))
            mean, median, std = calculate_image_sensitivity(image)  
        else:
            mean, median, std = calculate_image_sensitivity(image_jeckknifed)  
    else:
        mean, median, std = calculate_image_sensitivity(image)  
        
    if pbcor:
        image_path = os.path.dirname(image)
        image_basename = os.path.basename(image)
        image_pbcor = os.path.join(image_path, image_basename.replace('image', 'pbcor.image'))
        with fits.open(image_pbcor) as hdu_pbcor:
            data_pbcor = hdu_pbcor[0].data
            data_pbcor_masked = np.ma.masked_invalid(data_pbcor.reshape(ny, nx))
            pbcor_map = (data_masked.filled(0.0) / data_pbcor_masked.filled(np.inf)).reshape(ny, nx)
    
    if model_background:
        if filter_size is None:
            filter_size = 3#max(int(fwhm_pixel/1.5), 6)
        if box_size is None:
            box_size = int(fwhm_pixel*2)
        if debug:
            print('filter_size', filter_size)
            print('box_size', box_size)
        sigma_clip = SigmaClip(sigma=3.)
        # bkg_estimator = MedianBackground()
        bkg_estimator = SExtractorBackground() 
        bkg = Background2D(data_masked.data, (box_size, box_size), 
                filter_size=(filter_size, filter_size),
                mask=data_masked.mask,
                sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        background = bkg.background
    else:
        background = median
    data_masked_sub = data_masked - background

    # numerical aperture correction
    image_gaussian = make_gaussian_image((2.*np.ceil(a), 2.*np.ceil(a)), fwhm=[b,a], 
                offset=[0.5, 0.5], # 0.5 offset comes from central offset of photutils of version 0.5
                theta=theta/180.*np.pi) 
    flux_g = auto_photometry(image_gaussian, bmaj=b, bmin=a, theta=theta/180.*np.pi, 
                             methods='aperture', debug=False, aperture_correction=1.0,
                             beamsize=1.0)
    aperture_correction = 1/flux_g[0]
    
    print("known_sources_coords", known_sources_coords)
    print("known_sources_pixel", known_sources_pixel)
    print("known_sources_center", known_sources_center)
    seg_radius = 2.0*np.int(a)
    segments = RectangularAperture(known_sources_center, 2*seg_radius, 2*seg_radius, theta=0)
    segments_mask = segments.to_mask(method='center')
    
    flux_auto = []
    flux_err_auto = []
    # flux_snr_list = []
    for i,s in enumerate(segments_mask):
        if subtract_background:
            data_cutout = s.cutout(data_masked_sub)
        else:
            data_cutout = s.cutout(data_masked)
        flux_list, flux_err_list = auto_photometry2(data_cutout, bmaj=b, bmin=a, beamsize=beamsize,
                                    theta=theta/180*np.pi, debug=debug, methods=methods,
                                    aperture_correction=aperture_correction, rms=std)
        #flux_err_list.append(flux_snr)
        #flux_snr = flux_list / std
        if target_wave:
            flux_list = flux_interpolate(lam, flux_list, target_wave) 
        if pbcor:
            pbcor_pixel = pbcor_map[np.ceil(known_sources_center[i][0]), 
                                    np.ceil(known_sources_center[i][1])]
            flux_auto.append(np.array(flux_list) * 1000. / pbcor_pixel)
            flux_err_auto.append(np.array(flux_err_list) * 1000. / pbcor_pixel)
        else:
            flux_auto.append(np.array(flux_list) * 1000.)
            flux_err_auto.append(np.array(flux_err_list) * 1000.)
        flux_snr_list = np.array(flux_auto) / np.array(flux_err_auto)
    
    if calculate_radial_distance:
        radial_distance = []
        for coord_pixel in known_sources_center:
            pixel_dist = np.sqrt((coord_pixel[0] - nx/2.0)**2+
                                 (coord_pixel[1] - ny/2.0)**2)
            ang_dist = pixel_dist / deg2pixel * 3600 # in arcsec
            radial_distance.append(ang_dist)
    if debug:
        # visualize the results
        if ax is None:
            fig= plt.figure()
            ax = fig.add_subplot(111)
        ny, nx = data_masked.shape[-2:]
        scale = np.abs(header['CDELT1'])*3600
        x_index = (np.arange(0, nx) - nx/2.0) * scale
        y_index = (np.arange(0, ny) - ny/2.0) * scale
        #x_map, y_map = np.meshgrid(x_index, y_index)
        #ax.pcolormesh(x_map, y_map, data_masked)
        extent = [np.min(x_index), np.max(x_index), np.min(y_index), np.max(y_index)]
        ax.imshow(data_masked, origin='lower', extent=extent, interpolation='none', vmax=10.*std, 
                  vmin=-3.0*std)
        
        ax.text(0, 0, '+', color='r', fontsize=24, fontweight=100, horizontalalignment='center',
                verticalalignment='center')
        ellipse = patches.Ellipse((0.8*np.min(x_index), 0.8*np.min(y_index)), width=bmin, height=bmaj, 
                                  angle=theta, facecolor='orange', edgecolor=None, alpha=0.8)
        ax.add_patch(ellipse)
        if debug:
            ellipse = patches.Ellipse((0, 0), width=fov_scale*fov, height=fov_scale*fov, 
                                      angle=0, fill=False, facecolor=None, edgecolor='grey', 
                                      alpha=0.8)
            ax.add_patch(ellipse)
        ax.text(0.7*np.max(x_index), 0.9*np.min(y_index), 'std: {:.2f}mJy'.format(std*1000.), 
                color='white', fontsize=10, horizontalalignment='center', verticalalignment='center')
        for i,e in enumerate(zip(*known_sources_pixel)):
            yy = scale*(e[0]-ny/2.)
            xx = scale*(e[1]-ny/2.)
            ellipse = patches.Ellipse((yy, xx), width=2*b*scale, 
                                      height=2*a*scale, angle=theta, facecolor=None, fill=False, 
                                      edgecolor='black', alpha=0.8, linewidth=2)
            ax.add_patch(ellipse)
        plt.show()

    print('flux_auto', flux_auto)
    print('flux_snr_list', flux_snr_list)

    if calculate_radial_distance:
        return flux_auto, flux_snr_list, radial_distance
    # return flux_auto, flux_snr_list
    return flux_auto, flux_err_auto
 

def run_flux_measure_almared(fitsimages, summary_file=None):
    if os.path.isdir(fitsimages):
        fitsimages = glob.glob(os.path.join(fitsimages, "*.pbcor.fits"))
    
    if not os.path.isfile(summary_file):
        print('Initializing the output file')
        with open(summary_file, 'w+') as f:
            f.write("obj ra dec flux_aperture flux_err_aperture flux_gaussian flux_err_gaussian flux_peak flux_err_peak")
            f.write('\n')
    name_match = re.compile('uid___[\w]+\.(?P<name>[\w\.]+)_sci.spw')
    for image in fitsimages:
        obj = name_match.search(image).groupdict()['name']
        sf = source_finder2(image, threshold=5.0, return_flux=False)
        targets_coords = []
        for target in sf:
            targets_coords.append([target[0], target[1]])
        sources_flux, sources_flux_err = flux_measure2(image, targets_coords)
            
        if summary_file:
            detections_summary = open(summary_file, 'a+')
            for i in range(len(sf)):
                detections_summary.write('{} {:.6f} {:.6f} {} {} {} {} {} {}'.format(obj,
                                         targets_coords[i][0], targets_coords[i][1],
                                         sources_flux[i][0], sources_flux_err[i][0],
                                         sources_flux[i][1], sources_flux_err[i][1],
                                         sources_flux[i][2], sources_flux_err[i][2]))
                detections_summary.write('\n')
            detections_summary.close()
