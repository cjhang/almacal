# A collection of tools related with image analysis

# Author: Jianhang Chen
# Email: cjhastro@gmail.com
import numpy as np

def auto_photometry2(image, bmaj=1, bmin=1, theta=0, beamsize=None, debug=False, 
                    methods=['aperture','gaussian', 'peak'], aperture_correction=1.068,
                    rms=None):
    """autoumatically measure the flux with different methods
    
    aperture_size = 1.0 # in units of beam size
    """
    ysize, xsize = image.shape
    y_center, x_center = ysize/2., xsize/2.
    flux_list = []
    flux_error = []
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
        else:
            flux_list.append(flux_in_aper)
        # print('aperture.area;',aperture.area())
        flux_error.append(rms * np.sqrt(aperture.area()))
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
        else:
            flux_list.append(flux_fitted)
        flux_error.append(rms*np.sqrt(gaussian_area))
    if 'peak' in methods:
        # the peak value has the same units of one pixel
        # if beamsize:
            # flux_list.append(np.max(image))#)#/beamsize)
        # else:
        flux_list.append(np.max(image))
        flux_error.append(rms)
    if 'aperture_auto' in methods:
        radii = np.arange(0.2*bmaj/2, ysize/2./np.cos(theta), 0.5/np.cos(theta))

        apertures = [EllipticalAperture((x_center, y_center), a, a*bmin/bmaj, theta) for a in radii]
        extract_area = lambda x: x.area()
        area_apers = np.array(list(map(extract_area, apertures)))
        
        phot_table = aperture_photometry(image, apertures)
        extract_qtable = lambda x: x.data[0]
        flux_apers = np.array(list(map(extract_qtable, phot_table.columns[3:].values())))
        snr_apers = flux_apers / (rms*np.sqrt(area_apers))

        slopes = np.ones_like(snr_apers)
        slopes[1:] = np.diff(snr_apers)#/np.diff(area_apers) #Normalized slope
        slope_selection = slopes < 1e-4
        if debug:
            print('position:{}/{}'.format(np.min(np.where(slope_selection)[0]), len(radii)))
            # print('slopes', slopes)
        if len(flux_apers[slope_selection])<3:
            flux_apers_stable = flux_apers[-2]
            area_apers_stable = area_apers[-2]
        else:
            flux_apers_stable = flux_apers[slope_selection][0]
            area_apers_stable = area_apers[slope_selection][0]
        if debug:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            im = ax.imshow(image, origin='lower', interpolation='none')
            plt.colorbar(im, fraction=0.046, pad=0.04)
            for ap in apertures:
                im = ap.plot(color='gray', lw=2)
        if beamsize:
            flux_list.append(1.0*flux_apers_stable/beamsize)
        else:
            flux_list.append(flux_apers_stable)
        flux_error.append(rms*np.sqrt(area_apers_stable))
    return flux_list, flux_error

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
                                    theta=theta/180*np.pi, debug=False, methods=methods,
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
    return flux_auto, flux_snr_list
 
