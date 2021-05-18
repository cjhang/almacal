# functions related with almacal simulation
import numpy as np
from scipy import signal, ndimage 
from scipy.interpolate import CubicSpline
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

def gkern(bmin=1., bmaj=1, theta=0, size=21,):
    """
    creates gaussian kernel with side length l and a sigma of sig
    """
    size = np.max([size, 2*int(bmaj)+1])
    FWHM2sigma = 2.*np.sqrt(2*np.log(2))
    gkernx = signal.gaussian(size, std=bmin/FWHM2sigma).reshape(size, 1)
    gkerny = signal.gaussian(size, std=bmaj/FWHM2sigma).reshape(size, 1)

    kernel = np.outer(gkerny, gkernx)
    kernel_rot = ndimage.rotate(kernel, -1.*theta, reshape=False)
    
    return kernel_rot / np.sum(kernel_rot)

def aspectNC(s):
    return 4.4e5*(s/0.1+1)**(-2.5)

def loguniform(low, high, size=None):
    return 10**(np.random.uniform(low, high, size))

def aspectNC_pdf(radius=None, boundary=[0.01, 10], ):
    """differential number

    Args:
        radius: in arcsec
    """
    xtmp = np.linspace(boundary[0], boundary[1], 10000)
    xtmp_mid = 0.5*(xtmp[1:]+xtmp[:-1])
    ytmp = aspectNC(xtmp_mid)/3600**2*np.pi*radius**2 #* (xtmp[1]-xtmp[0])
    pdf_discrete = ytmp/np.sum(ytmp)
    cs = CubicSpline(xtmp_mid, pdf_discrete)
    return cs

def aspect_sampler(n, niter=100, fluxrange=[0.02, 10], radius=15):
    x = loguniform(np.log10(fluxrange[0]), np.log10(fluxrange[1]), niter) 
    # print(x)
    y = aspectNC(x)/3600**2*np.pi*radius**2 #* (xtmp[1]-xtmp[0])
    y_pdf = y / np.sum(y)
    return np.random.choice(x, n, p=y_pdf)

def make_random_source(direction, direction_units='deg',freq=None, n=None, radius=1, 
        prune=True, prune_threshold=3., debug=False, savefile=None, clname=None,
        fluxrange=[0, 1], fluxunit='mJy', known_sources=None,
        sampler=np.random.uniform, sampler_params={}, budget=None):
    """This function used to add random source around a given direction
        
    This is a proximate solution, which treat the field of view as a plane.
    Then calculate the longtitude and latitue using sin and cos functions.

    Args:
        direction: the pointing direction
        n: the number of points
        radius: the largest radius that the points can be put, in arcsec
                it can be a list whose format is [r_min, r_max] both in arcsec
        prune: remove close pairs, make it easier for pointsource testing
        savefile: save the coordination files
        clname: the component list filename
        budget: The total flux density, conflict with flux
        
    Example:
        direction = 'J2000 13h03m49.215900s -55d40m31.60870s'
   """
    # generate a random radius
    # print(direction)
    skycoord = SkyCoord(direction[5:])
    #if isinstance(direction, str):
    #    try:
    #        skycoord = SkyCoord(direction[5:])
    #    except:
    #        skycoord = SkyCoord(direction)
    #elif isinstance(direction, list):
    #    skycoord = SkyCoord(*direction, units=direction_units)
    #else:
    #    print("{} as a direction is not support!".format(direction))
    if debug:
        print('fluxrange', fluxrange)
        print('radius', radius)
    if isinstance(fluxrange, (int, float)):
        fluxrange = [fluxrange, fluxrange]
    if n:
        theta = 2*np.pi*np.random.uniform(0, 1, n)
        if isinstance(radius, (float, int)):
            rho = radius * np.sqrt(np.random.uniform(0, 1, n))
        elif isinstance(radius, (list, np.ndarray)):
            rho = np.diff(radius)*np.sqrt(np.random.uniform(0, 1, n)) + radius[0]
        flux_sampling = np.array([sampler(**sampler_params) for i in range(n)])
        flux_input = flux_sampling * np.diff(fluxrange) + fluxrange[0]
        if debug:
            print('flux range {} [{}]'.format(fluxrange, fluxunit))
    elif budget:
        theta = []
        rho = []
        flux_input = []
        total_input = 0 # in mJy
        while True:
            theta_one = 2*np.pi*np.random.uniform(0, 1)
            rho_one = radius * np.sqrt(np.random.uniform(0, 1))
            # flux_sampling = sampler(scaler, 1)
            flux_sampling = sampler(**sampler_params)
            if len(flux_sampling) == 1:
                flux_input_one = flux_sampling
            else:
                flux_input_one = np.random.choice(flux_sampling)
            
            total_input += flux_input_one
            if total_input < budget:
                total_input += flux_input_one
                theta.append(theta_one)
                rho.append(rho_one)
                flux_input.append(flux_input_one[0])
            else:
                total_input -= flux_input_one
                break
        theta, rho, flux_input = np.array(theta), np.array(rho), np.array(flux_input)
        if debug:
            print('flux range {} [{}]'.format(fluxrange, fluxunit))
            print('flux input {}, budget: {}'.format(total_input, budget))
            print(flux_input)
    else:
        raise ValueError("You need to specify number of points or the total budget!")
    # plt.hist(flux_input)
    # plt.show()
    # return

    delta_ra = np.array(rho) * np.cos(theta)/(np.cos(skycoord.dec).value)
    delta_dec = np.array(rho) * np.sin(theta)

    if prune:
        n_prune = 0
        selected_ra = [delta_ra[0],] 
        selected_dec = [delta_dec[0],]
        for i in range(1, n):
            coord_i = SkyCoord(ra=delta_ra[i]*u.arcsec, dec=delta_dec[i]*u.arcsec)
            coords_selected = SkyCoord(ra=np.array(selected_ra)*u.arcsec, dec=np.array(selected_dec)*u.arcsec)
            separation_with_select = coord_i.separation(coords_selected)
            if np.sum(separation_with_select<prune_threshold*u.arcsec)>=1:
                if debug:
                    print("pruning one point...")
                n_prune += 1
                continue
            selected_ra.append(delta_ra[i]) 
            selected_dec.append(delta_dec[i])
        delta_ra = np.array(selected_ra)
        delta_dec = np.array(selected_dec)
        print("Pruned {} sources.".format(n_prune))
    ra_random = delta_ra*u.arcsec + skycoord.ra
    dec_random = delta_dec*u.arcsec+ skycoord.dec
   
    if (known_sources is not None) and (len(known_sources) > 0):
        known_sources_coords_ra = []
        known_sources_coords_dec = []
        distance = prune_threshold*u.arcsec
        # print('distance', distance)
        for source in known_sources:
            known_sources_coords_ra.append(source[0])
            known_sources_coords_dec.append(source[1])
        known_sources_coords = SkyCoord(ra=known_sources_coords_ra, dec=known_sources_coords_dec, 
                                        unit='arcsec')

        selected_after_known_ra = []
        selected_after_known_dec = []
        for ra,dec in zip(ra_random, dec_random):
            if (np.sum(np.abs(ra - known_sources_coords.ra) < distance) > 0 and 
                np.sum(np.abs(dec - known_sources_coords.dec) < distance) > 0):
                print('skipping one source')
                continue
            selected_after_known_ra.append(ra)
            selected_after_known_dec.append(dec)
        ra_random = selected_after_known_ra
        dec_random = selected_after_known_dec

    if savefile:
        with open(savefile, 'w+') as sfile:
            sfile.write('# ra[arcsec]  dec[arcsec]  flux[{}]\n'.format(fluxunit))
            for ra, dec, flux in zip(ra_random, dec_random, flux_input):
                sfile.write('{:.5f} {:.6f} {:.4f}\n'.format(ra.value, dec.value, flux))
        # np.savetxt(savefile, [[ra_random], [dec_random], [flux_random]], 
                # header='#ra[arcsec]  dec[arcsec]  flux[{}]'.format(fluxunit))

    if clname:
        # generate component list
        skycoord_list = SkyCoord(ra=ra_random, dec=dec_random, unit='arcsec')
        f = lambda x: ('J2000 '+x.to_string('hmsdms')).encode('utf-8')
        direction_list = list(map(f, skycoord_list))
        os.system('rm -rf {}'.format(clname))
        cl.done()
        for d,f in zip(direction_list, flux_input):
            cl.addcomponent(dir=d, flux=f, fluxunit=fluxunit, 
                            freq=freq, shape='point',index=0)
        cl.rename(clname)
        cl.done()
        return clname
    else:
        return [ra_random, dec_random, flux_input]

def add_random_sources(vis=None, fitsimage=None, mycomplist=None, outdir='./', 
        make_image=True, outname=None, debug=False, uvtaper_scale=None,
        **kwargs):
    """
    The sources will be injected to the original file, so only pass in the copied data!
    Args:
        radius : in arcsec
        budget: can be a single value, or a list, [mean, low_boundary, high_boundary]
        flux: units in mJy

    Notes:
        1. To make the reading and writing more efficiently, the model will be added directly
           into original measurement. It is suggested to make a copy before calling this func
    """
    
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
   
    if vis:
        ft(vis=vis, complist=mycomplist)
        uvsub(vis=vis, reverse=True)
        suffix = '.cont.auto'
        make_cont_img(vis, outdir=outdir, clean=True, niter=1000, 
                      only_fits=True, uvtaper_scale=uvtaper_scale, pblimit=-0.01,
                      fov_scale=2.0, basename=outname, suffix='')
        delmod(vis=vis)
        clearcal(vis=vis)
    if fitsimage:
        # add random sources directly into the image
        hdu = fits.open(fitsimage)
        header = hdu[0].header
        wcs = WCS(header)
        data = hdu[0].data
        ny, nx = data.shape[-2:]
        data_masked = np.ma.masked_invalid(data.reshape(ny, nx))
        # mean, median, std = sigma_clipped_stats(data_masked, sigma=10.0)  

        # for DAOStarFinder, in pixel space
        pixel_scale = 1/np.abs(header['CDELT1'])
        fwhm = header['BMAJ']*3600*u.arcsec
        fwhm_pixel = header['BMAJ']*pixel_scale
        a, b = header['BMAJ']*pixel_scale, header['BMIN']*pixel_scale
        ratio = header['BMIN'] / header['BMAJ']
        theta = header['BPA']
        beamsize = np.pi*a*b/(4*np.log(2))

        # kernel = Gaussian2DKernel(stddev=0.25*(a+b))
        kernel = gkern(bmaj=a, bmin=b, theta=theta)
        # print('mycomplist', mycomplist)
        # hdr = wcs.to_header()
        # hdr['OBSERVER'] = 'Your name'
        # hdr['COMMENT'] = "Here's some comments about this FITS file."
        mycomplist_pixels = []
        blank_image = np.zeros((ny, nx))
        # mycomplist = make_random_source(mydirection, freq=myfreq, n=n, radius=radius, debug=debug, 
                                        # fluxrange=fluxrange, savefile=savefile_fullpath, 
                                        # known_sources=known_sources, budget=budget, **kwargs) 
        for ra, dec, flux in zip(*mycomplist):
            ra_pix, dec_pix = map(np.around, skycoord_to_pixel(SkyCoord(ra, dec, unit='arcsec'), wcs))
            # print(ra_pix, dec_pix, flux)
            mycomplist_pixels.append([int(ra_pix), int(dec_pix), flux])
            blank_image[int(dec_pix), int(ra_pix)] = flux
        fake_image = convolve(blank_image, kernel)/1000*beamsize + data_masked # convert mJy into Jy before added into image
        fake_image = fake_image.filled(np.nan)
        hdu = fits.PrimaryHDU(header=header, data=fake_image)
        hdu.writeto(os.path.join(outdir , outname+'.fits'), overwrite=True)
        return 

def subtract_sources(vis, complist=None, ):
    md = msmdtool()
    if not md.open(vis):
        raise ValueError("Failed to open {}".format(vis))
    phasecenter = md.phasecenter()
    mydirection = phasecenter['refer'] +' '+ SkyCoord(phasecenter['m0']['value'], phasecenter['m1']['value'], unit="rad").to_string('hmsdms')
    myfreq = "{:.2f}GHz".format(np.mean(read_spw(vis)))
    ft(vis=vis, complist=complist)
    uvsub(vis=vis, reverse=False)
    vis_new = vis+'.sourcessub'
    split(vis=vis, datacolumn='corrected', outputvis=vis_new)

def make_gaussian_image(shape, fwhm=None, sigma=None, area=1., offset=(0,0), theta=0):
    """make a gaussian image for testing

    theta: in rad, rotating the gaussian counterclock wise
    """
    image = np.zeros(shape, dtype=float)
    yidx, xidx = np.indices(shape)
    yrad, xrad = yidx-shape[0]/2., xidx-shape[1]/2.
    y = xrad*np.cos(theta) + yrad*np.sin(theta)
    x = yrad*np.cos(theta) - xrad*np.sin(theta)
    if fwhm:
        fwhm2sigma = 2.35482
        sigma = (np.array(fwhm) / fwhm2sigma).tolist()
    if isinstance(sigma, (list, tuple)):
        ysigma, xsigma = sigma
    elif isinstance(sigma, (int, float)):
        ysigma = xsigma = sigma
    if isinstance(offset, (list, tuple)):
        y_offset, x_offset = offset
    elif isinstance(offset, (int, float)):
        y_offset = x_offset = offset
    flux = area * np.exp(-(x-x_offset)**2/2./xsigma**2 - (y-y_offset)**2/2./ysigma**2) / (2*np.pi*xsigma*ysigma)

    return flux

def auto_photometry(image, bmaj=1, bmin=1, theta=0, beamsize=None, debug=False, 
                    methods=['aperture','gaussian', 'peak'], aperture_correction=1.068):
    """automatically measure the flux with different methods
    
    aperture_size = 1.0 # in units of beam size
    """
    
    ysize, xsize = image.shape
    y_center, x_center = ysize/2., xsize/2.
    flux_list = []
    if beamsize is None:
        print("Warning: Unit the sum and maximum, without dividing the beam size.")
    sigma2FWHM = 2.35482
    # Aperture Photometry
    if 'aperture' in methods:
        aperture_size=1.0
        aperture = EllipticalAperture((x_center, y_center), aperture_size*bmaj, aperture_size*bmin, theta)
        extract_area = lambda x: x.area()
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
    if 'gaussian' in methods:
        yidx, xidx = np.indices((ysize, xsize))
        yrad, xrad = yidx-ysize/2., xidx-xsize/2.

        image_norm = image / np.max(image)
        p_init = models.Gaussian2D(amplitude=1, x_stddev=bmaj/sigma2FWHM, y_stddev=bmin/sigma2FWHM, 
                theta=theta, bounds={"x_mean":(-2.,2.), "y_mean":(-2.,2.), 
                                     "x_stddev":(0., xsize), "y_stddev":(0., ysize)})
        fit_p = fitting.LevMarLSQFitter()
        p = fit_p(p_init, xrad, yrad, image_norm, 
                  weights=1/(yrad**2+xrad**2+(0.25*bmaj)**2+(0.25*bmin)**2))
        flux_fitted = 2*np.max(image)*np.pi*p.x_stddev.value*p.y_stddev.value*p.amplitude.value

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
    if 'peak' in methods:
        # the peak value has the same units of one pixel
        # if beamsize:
            # flux_list.append(np.max(image))#)#/beamsize)
        # else:
        flux_list.append(np.max(image))
    if 'aperture_auto' in methods:
        radii = np.arange(0.2*bmaj/2, ysize/2./np.cos(theta), 0.5/np.cos(theta))

        apertures = [EllipticalAperture((x_center, y_center), a, a*bmin/bmaj, theta) for a in radii]
        extract_area = lambda x: x.area()
        area_apers = np.array(list(map(extract_area, apertures)))
        
        phot_table = aperture_photometry(image, apertures)
        extract_qtable = lambda x: x.data[0]
        flux_apers = np.array(list(map(extract_qtable, phot_table.columns[3:].values())))

        slopes = np.ones_like(area_apers)
        slopes[1:] = np.diff(flux_apers)/np.max(flux_apers)/np.diff(area_apers) #Normalized slope
        slope_selection = slopes < 0.001
        if debug:
            print('position:{}/{}'.format(np.min(np.where(slope_selection)[0]), len(radii)))
            # print('slopes', slopes)
        if len(flux_apers[slope_selection])<3:
            flux_apers_stable = flux_apers[-2]
        else:
            flux_apers_stable = flux_apers[slope_selection][0]
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
    return flux_list

def source_finder(fitsimage, outdir='./', sources_file=None, savefile=None, model_background=True, 
                  threshold=5.0, debug=False, algorithm='DAOStarFinder', return_image=False,
                  filter_size=None, box_size=None, methods=['aperture', 'gaussian','peak'],
                  subtract_background=False, known_sources=None, figname=None, ax=None, pbcor=False,
                  fov_scale=2.0, mask_threshold=5., second_check=True, baseimage=None,
                  return_image_cuts=False):
    """finding point source in the image

    This is a two stage source finding algorithm. First, DAOStarFinder or find_peak will be used to find
    the sources. Then, additional total and peak SNR checking is used to confirm the detection. The 
    second stage can be configured by 

    mask_threshold: the mask size of known_sources
    """
    with fits.open(fitsimage) as hdu:
        header = hdu[0].header
        wcs = WCS(header)
        data = hdu[0].data
        ny, nx = data.shape[-2:]
        data_masked = np.ma.masked_invalid(data.reshape(ny, nx))
        mean, median, std = sigma_clipped_stats(data_masked, sigma=10.0)  
        freq = header['CRVAL3']
        lam = (const.c/(freq*u.Hz)).decompose().to(u.um)
        fov = 1.02 * (lam / (12*u.m)).decompose()* 206264.806
        fcenter = [header['CRVAL1'], header['CRVAL2']] # in degrees
        pixel_center = [nx/2., ny/2.]
    if pbcor:
        fitsimage_path = os.path.dirname(fitsimage)
        fitsimage_basename = os.path.basename(fitsimage)
        fitsimage_pbcor = os.path.join(fitsimage_path, fitsimage_basename.replace('image', 
                                       'pbcor.image'))
        if not os.path.isfile(fitsimage_pbcor):
            print("No primary beam corrected data found, the final flux can be wrong!")
            pbcor = False
    if pbcor:
        # load the pbcorred data
        with fits.open(fitsimage_pbcor) as hdu_pbcor:
            data_pbcor = hdu_pbcor[0].data
            data_pbcor_masked = np.ma.masked_invalid(data_pbcor.reshape(ny, nx))
            pbcor_map = (data_masked.filled(0.0) / data_pbcor_masked.filled(np.inf)).reshape(ny, nx)
    
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

    # read known_sources
    known_mask = np.full((ny, nx), 0, dtype=bool)
    if baseimage:
        known_sources = source_finder(baseimage)
    if (known_sources is not None) and (len(known_sources) > 0):
        known_sources_coords_ra = []
        known_sources_coords_dec = []
        distance = mask_threshold*u.arcsec
        # print('distance', distance)
        for source in known_sources:
            known_sources_coords_ra.append(source[0])
            known_sources_coords_dec.append(source[1])
        known_sources_coords = SkyCoord(ra=known_sources_coords_ra, dec=known_sources_coords_dec, 
                                        unit='arcsec')
        known_sources_pixel = skycoord_to_pixel(known_sources_coords, wcs)
        for p in list(zip(*known_sources_pixel)):
            aper = CircularAperture(p, 2.0*a)
            aper_mask = aper.to_mask(method='center')
            known_mask = np.bitwise_or(known_mask, aper_mask[0].to_image((ny,nx)).astype(bool))

        if False:
            plt.figure()
            plt.imshow(known_mask,origin='lower')
            # for p in known_sources_pixel:
                # aperture = CircularAperture(p, r=2.*a)
                # aperture.plot(color='white', lw=2)
            plt.show()

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
    data_masked_sub[known_mask] = std
    
    # find stars
    if algorithm == 'DAOStarFinder': # DAOStarFinder
        daofind = DAOStarFinder(fwhm=fwhm_pixel, threshold=threshold*std, ratio=ratio, 
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
                box_size=fwhm_pixel, mask=known_mask)
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
        segments = RectangularAperture(sources_found_center_candidates, 2.*a, 2.*a, theta=0)
        segments_mask = segments.to_mask(method='center')
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
                if pbcor:
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


    if sources_file:
        sources_input = np.loadtxt(sources_file, skiprows=1)
        if len(sources_input.shape) == 1:
            sources_input_coords = SkyCoord(ra=sources_input[0], dec=sources_input[1], unit='arcsec')
            flux_input = sources_input[-1]
        elif len(sources_input.shape) > 1:
            sources_input_coords = SkyCoord(ra=sources_input[:,0], dec=sources_input[:,1], 
                                            unit='arcsec')
            flux_input = sources_input[:,-1]
        else:
            raise ValueError('Unsupported file: {}'.format(sources_file))
        # flux_input = sources_input[:,-1]
        # sources_input_coords = SkyCoord(ra=sources_input[:,0]*u.arcsec, 
                                        # dec=sources_input[:,1]*u.arcsec)
        sources_input_pixels = skycoord_to_pixel(sources_input_coords, wcs)
        sources_input_center = zip(*sources_input_pixels)
        
        #automatically aperture photometry
        flux_input_auto = []
        segments = RectangularAperture(sources_input_center, 2.*a, 2.*a, theta=0)
        segments_mask = segments.to_mask(method='center')
        for s in segments_mask:
            if subtract_background:
                data_cutout = s.cutout(data_masked_sub)
            else:
                data_cutout = s.cutout(data_masked)
            flux_list = auto_photometry(data_cutout, bmaj=b, bmin=a, beamsize=beamsize,
                                        theta=theta/180*np.pi, debug=False, methods=methods,
                                        aperture_correction=aperture_correction)
            flux_input_auto.append(np.array(flux_list) * 1000) #change to mJy

        if sources_found:
            sources_input_pixels = skycoord_to_pixel(sources_input_coords, wcs)

            idx_found, idx_input, d2d, d3d = sources_input_coords.search_around_sky(
                    sources_found_coords, fwhm)
            # sources_input_found = [flux_input[idx_input], np.array(flux_auto)[idx_found]]
            idx_input_comp = np.array(list(set(range(len(flux_input))) - set(idx_input)), dtype=int)
            idx_found_comp = np.array(list(set(range(len(flux_auto))) - set(idx_found)), dtype=int)
            idxs = [idx_input, idx_found, idx_input_comp, idx_found_comp]
        else:
            # sources_input_found = [np.array([]), np.array([])]
            idxs = [np.array([]), np.array([]), np.array([]), np.array([])]

        if debug:
            print(flux_input)
            print(flux_input_auto)

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
        ax.imshow(data_masked, origin='lower', extent=extent, interpolation='none')
        
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
        
        if sources_file:
            # show the input sources
            for i,e in enumerate(sources_input_center):
                ellipse = patches.Ellipse(e, width=3*b, height=3*a, angle=theta, facecolor=None, 
                                          fill=False, edgecolor='orange', alpha=0.8, linewidth=1)
                ax.add_patch(ellipse)
                if debug:
                    ax.text(e[0], e[1], flux_input[i])

        # show the sources found
        if known_sources:
            for i,e in enumerate(zip(*known_sources_pixel)):
                yy = scale*(e[0]-ny/2.)
                xx = scale*(e[1]-ny/2.)
                ellipse = patches.Ellipse((yy, xx), width=2*b*scale, 
                                          height=2*a*scale, angle=theta, facecolor=None, fill=False, 
                                          edgecolor='black', alpha=0.8, linewidth=2)
                ax.add_patch(ellipse)
        if sources_found:
            for i,e in enumerate(sources_found_center):
                yy = scale*(e[0]-ny/2.)
                xx = scale*(e[1]-ny/2.)
                ellipse = patches.Ellipse((yy, xx), width=3*b*scale, 
                                          height=3*a*scale, angle=theta, facecolor=None, fill=False, 
                                          edgecolor='white', alpha=0.8, linewidth=2)
                ax.add_patch(ellipse)
                ax.text(1.1*yy, 1.0*xx, "({})\n{:.2f}mJy".format(i, flux_auto[i][0]), 
                        color='white', fontsize=12,)
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
    if return_image:
        return data_masked.data, background
    if return_image_cuts:
        data_cutout_list = []
        segments_true = RectangularAperture(sources_found_center, 2.*a, 2.*a, theta=0)
        segments_true_mask = segments_true.to_mask(method='center')
        for i,s in enumerate(segments_true_mask):
            if subtract_background:
                data_cutout = s.cutout(data_masked_sub)
            else:
                data_cutout = s.cutout(data_masked)
            data_cutout_list.append(data_cutout)    
        return data_cutout_list


    if sources_file:
        return np.array(flux_input), np.array(flux_input_auto), np.array(flux_auto), idxs    
    if len(flux_auto)>0:
        return list(zip(sources_found_coords.ra.to(u.arcsec).value, 
                        sources_found_coords.dec.to(u.arcsec).value, flux_auto, flux_snr_list))
    else:
        return []

def gen_sim_images(mode='image', vis=None, imagefile=None, n=20, repeat=1, 
                    snr=(1,20), fov_scale=1.5, outdir='./', basename=None,
                    uvtaper_scale=None, budget=None,
                    debug=False, **kwargs):
    """generate the fake images with man-made sources

    Args:
     mode: image or uv
    """
    # make a copy of the original file
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))

    if mode == 'image':
        if basename is None:
            basename = os.path.basename(imagefile)
        if not imagefile:
            raise ValueError("The image file must be given!")
        with fits.open(imagefile) as hdu:
            data = hdu[0].data
            header = hdu[0].header
        freq = header['CRVAL3']
        lam = (const.c/(freq*u.Hz)).decompose().to(u.um)
        fov = 1.02 * (lam / (12*u.m)).decompose().value * 206264.806
        ny, nx = data.shape[-2:]
        data_masked = np.ma.masked_invalid(data.reshape(ny, nx))
        mean, median, std = sigma_clipped_stats(data_masked, sigma=5.0)  

        # pixel_scale = 1/np.abs(header['CDELT1'])
        # fwhm = header['BMAJ']*3600*u.arcsec
        # fwhm_pixel = header['BMAJ']*pixel_scale
        # a, b = header['BMAJ']*pixel_scale, header['BMIN']*pixel_scale
        # ratio = header['BMIN'] / header['BMAJ']
        # theta = header['BPA']
        # beamsize = np.pi*a*b/(4*np.log(2))
        #gaussian_norm = 2*np.pi*a*b/2.35482**2
        refer = 'J'+str(int(header['EQUINOX']))
        mydirection = refer +' '+ SkyCoord(header['CRVAL1'], header['CRVAL2'], 
                                          unit="deg").to_string('hmsdms')
        myfreq = "{:.2f}GHz".format(header['CRVAL3']/1e9)

        sensitivity = std * 1000 # in mJy/pixel
        fluxrange = np.array(snr) * sensitivity
        known_sources = source_finder(imagefile, fov_scale=fov_scale)
        for i in range(repeat):
            print('run {}'.format(i))
            basename_repeat = basename + '.run{}'.format(i)
            complist_file = os.path.join(outdir, basename_repeat+'.txt')
            mycomplist = make_random_source(mydirection, freq=myfreq, 
                    radius=0.9*0.5*fov_scale*fov, # 0.9 is to compensate the optimal imsize 
                    debug=debug, fluxrange=fluxrange, savefile=complist_file, n=n, 
                    sampler=np.random.uniform, sampler_params={},
                    known_sources=known_sources, budget=budget, **kwargs) 
            add_random_sources(fitsimage=imagefile, mycomplist=mycomplist,
                    outdir=outdir, outname=basename_repeat, debug=debug, **kwargs)
    # adding source in uv has not been test for new scheme
    if mode == 'uv':
        if basename is None:
            basename = os.path.basename(vis)
        if not vis:
            raise ValueError("The visibility must be given!")
        if debug:
            print('basename', basename)
        # vistmp = os.path.join(outdir, basename+'.tmp')
        # split(vis=vis, outputvis=vistmp, datacolumn='data')
        # # read information from vis 
        md = msmdtool()
        if not md.open(vis):
            raise ValueError("Failed to open {}".format(vis))
        phasecenter = md.phasecenter()
        mydirection = phasecenter['refer'] +' '+ SkyCoord(phasecenter['m0']['value'], 
                        phasecenter['m1']['value'], unit="rad").to_string('hmsdms')
        freq_mean = np.mean(read_spw(vis))
        myfreq = "{:.2f}GHz".format(freq_mean)
        if debug:
            print(mydirection)
            print(myfreq)
        tb.open(vis + '/ANTENNA')
        antenna_diameter_list = tb.getcol('DISH_DIAMETER')
        tb.close()
        antenna_diameter = np.max(antenna_diameter_list) * u.m
        wavelength = const.c / (freq_mean * u.GHz) # in um
        fov = (fov_scale * 1.22 * wavelength / antenna_diameter * 206265).decompose().value
        if debug:
            print('fov', fov)
            print('radius', 0.5*fov)
        if imagefile is None:
            imagefile = os.path.join(outdir, basename+'.image.fits')
            if os.path.isfile(imagefile):
                print("Finding the default imagefile: {}".format(imagefile))
            else:
                if debug:
                    print('imagefile', imagefile)
                    print("No image file founded, trying to produce image instead!")
                make_cont_img(vis, outdir=outdir, clean=True, niter=1000, suffix='',
                              only_fits=True, uvtaper_scale=uvtaper_scale, pblimit=-0.01,
                              fov_scale=fov_scale, basename=basename)

        # im_head = imhead(imagefile)
        # im_beam = im_head['restoringbeam']
        # im_incr = im_head['incr']
        im_info = imstat(imagefile)
        # beamsize = np.pi*a*b/(4*np.log(2))
        # beamsize = np.pi/(4*np.log(2))* im_beam['major']['value'] * im_beam['minor']['value'] / (im_incr[0]/np.pi*180*3600)**2
        # gaussian_norm = 2.0*np.pi*im_beam['major']['value'] * im_beam['minor']['value'] / 2.35482**2
        sensitivity = im_info['rms'] * 1000 # in mJy/pixel
        # flux_base = rms[0] #* 2*np.pi # convert the peak flux density into brightness
        if debug:
            print('sensitivity', sensitivity)
            # print('gaussian_norm', gaussian_norm)
            # print('beamsize', beamsize)

        fluxrange = np.array(snr) * sensitivity
        known_sources = source_finder(imagefile, fov_scale=fov_scale)
        clname_fullpath = os.path.join(outdir, basename+'.cl')
        for i in range(repeat):
            print('run {}'.format(i))
            basename_repeat = basename + '.run{}'.format(i)
            complist_file = os.path.join(outdir, basename_repeat+'.txt')
            if debug:
                print('basename_repeat', basename_repeat)
                print('complist_file', complist_file)
            mycomplist = make_random_source(mydirection, freq=myfreq, 
                    radius=[0.45*fov_scale*fov, 0.45*fov*fov_scale], 
                    # radius=0.9*0.5*fov_scale*fov, # 0.9 is to compensate the optimal imsize 
                    debug=debug, fluxrange=fluxrange, savefile=complist_file, n=n, 
                    sampler=np.random.uniform, sampler_params={}, clname=clname_fullpath,
                    known_sources=known_sources, budget=budget, **kwargs) 
            add_random_sources(vis=vis, mycomplist=mycomplist,
                    outdir=outdir, outname=basename_repeat, debug=debug, **kwargs)
            rmtables(clname_fullpath)
    # if budget:
    # #if snr is None:
    # #TODO: better physical motivated
        # # Limitation from extragalactic background
        # EBL = 14 # 14-18Jy/deg2
        # radius = 0.5*fov
        # budget_mean = (np.pi*(radius*u.arcsec)**2 * 20*u.Jy/u.deg**2).to(u.mJy).value
        # print('budget_mean', budget_mean)
        # budget_sigma = 0.5
        # fluxrange = [1.0, 0.0001]
        # for i in range(repeat):
            # # generate the budget for each run
            # budget = budget_sigma * np.random.randn() + budget_mean

            # basename_new = basename+'.run{}'.format(i)
            # if mode == 'uv':
                # add_random_sources(vis=vistmp, n=None, budget=budget, radius=0.5*fov, outdir=outdir,
                        # uvtaper_scale=uvtaper_scale, basename=basename_new, known_file=known_file, 
                        # inverse_image=inverse_image, fluxrange=fluxrange, debug=debug, **kwargs)
            # elif mode == 'image':
                # add_random_sources(fitsimage=fitsimage, n=None, budget=budget, radius=0.5*fov, outdir=outdir,
                        # uvtaper_scale=uvtaper_scale, basename=basename_new, known_file=known_file, 
                        # inverse_image=inverse_image, fluxrange=fluxrange, debug=debug, **kwargs)

def calculate_sim_images(simfolder, vis=None, baseimage=None, n=20, repeat=10, 
        basename=None, savefile=None, fov_scale=1.5, second_check=True,
        threshold=5.0, plot=False, snr_mode='peak', **kwargs):
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

    known_sources = source_finder(baseimage, fov_scale=fov_scale)
    
    flux_match = re.compile('(?P<obj>J\d*[+-]\d*)_(?P<band>B\d+)\w*.snr(?P<snr>\d+.\d+).run(?P<run>\d+)')
    if basename is None:
        basename = os.path.basename(baseimage)
        print('basename', basename)
    # all_fake_images = glob.glob(basename+'*{}.fits'.format(suffix))

    flux_input_list = []
    flux_input_autolist = []
    flux_input_foundlist = []
    flux_found_autolist = []
    snr_input_list = []
    # snr_inputfound_comp = []
    detection_input_array = np.array([[0., 0.]])
    detection_found_array = np.array([[0., 0.]])
    for run in np.arange(repeat):
        n_input = 0
        n_found = 0
        # flux = s*sensitivity
        # flux_input_list.append(flux)
        #print('SNR: {}, RUN:{}'.format(s, run))
        #print('flux:', flux)
        simimage = "{basename}.run{run}.fits".format(basename=basename, run=run)
        simimage_sourcefile = "{basename}.run{run}.txt".format(basename=basename, run=run)
        simimage_fullpath = os.path.join(simfolder, simimage)
        simimage_sourcefile_fullpath = os.path.join(simfolder, simimage_sourcefile)
        # try:
        sf_return = source_finder(simimage_fullpath, sources_file=simimage_sourcefile_fullpath, 
                known_sources=known_sources, threshold=threshold, second_check=second_check, 
                **kwargs)
        # except:
            # continue
        # if sf_return == 0:
            # continue
        flux_input, flux_input_auto, flux_found_auto, idxs = sf_return 
        # print('flux_input', flux_input)
        # print('flux_input_auto', flux_input_auto)
        # print('flux_found_auto', flux_found_auto)
        # print('idxs', idxs)
        
        if len(idxs[0])<1:
            print("Skip run{}".format(run))
            continue
        flux_input_list.append(flux_input)
        flux_input_autolist.append(flux_input_auto)
        flux_found_autolist.append(flux_found_auto)


        # print('peak', flux_input_auto[:,2] / rms_flux)
        # print('integrated', flux_input / rms_flux)
        #calculate the snr
        if snr_mode == 'peak':
            snr_input = flux_input_auto[:,2] / rms_flux
            snr_input_list.append(snr_input)
        elif snr_mode == 'integrated':
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
            detection_input_array = np.vstack([detection_input_array, 
                np.array(zip(snr_input_found, np.ones_like(snr_input_found)))])
        if len(snr_input_failed) > 0:
            detection_input_array = np.vstack([detection_input_array, 
                np.array(zip(snr_input_failed, np.zeros_like(snr_input_failed)))])
        if len(snr_found_input) > 0:
            detection_found_array = np.vstack([detection_found_array, 
                np.array(zip(snr_found_input, np.ones_like(snr_found_input)))])
        if len(snr_found_fake) > 0:
            detection_found_array = np.vstack([detection_found_array, 
                np.array(zip(snr_found_fake, np.zeros_like(snr_found_fake)))])

    # print(snr_input_list)
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
        with open(savefile, 'w+') as fp:
           json.dump(data_saved, fp)
    if plot:
        plot_sim_results(data_saved)
    # return data_saved
    #flux_list, flux_peak_list, flux_found_list, completeness_list

def plot_sim_results(data=None, jsonfile=None, snr = np.arange(1.0, 10, 0.1), plot=True,
                     mode='median'):
    if jsonfile:
        with open(jsonfile) as jf:
            data = json.load(jf)
    snr_input = np.array(data['snr_input'])
    flux_input = np.array(data['flux_input'])
    flux_input_aperture = np.array(data['flux_input_aperture'])
    flux_input_gaussian = np.array(data['flux_input_gaussian'])
    flux_aperture = np.array(data['flux_aperture'])
    flux_gaussian = np.array(data['flux_gaussian'])
    detection_snr = np.array(data['detection_snr'])
    detection_input_array = np.array(data['detection_input_array'])
    detection_found_array = np.array(data['detection_found_array'])

    # print('detection_array', detection_array)
    completeness_list = []
    fake_rate_list = []
    # print(detection_input_array)
    # print(np.array(detection_input_array).shape)
    for i in range(1, len(snr)):
        s = snr[i]
        s_b = snr[i-1]
        # calculate the completeness
        snr_select = np.bitwise_and((detection_input_array[:, 0]<s), (detection_input_array[:, 0]>s_b))
        n_input = np.sum(snr_select)
        n_found = np.sum(np.array(detection_input_array[:,1][snr_select]))
        completeness_list.append(1.0*n_found/n_input)
        
        # calculate the fake detaction rate
        snr_select2 = np.bitwise_and((detection_found_array[:, 0]<s), (detection_found_array[:, 0]>s_b))
        n_found2 = np.sum(snr_select2)
        n_fake = n_found2 - np.sum(np.array(detection_found_array[:,1][snr_select2]))
        fake_rate_list.append(1.0*n_fake/n_found2)

    # calculate flux boosting and completeness
    snr_mid = 0.5*(snr[1:] + snr[:-1])
    aperture_mean = []
    gaussian_mean = []
    completeness_list = []
    fake_rate_list = []
    # total_found_list = []
    # total_fake_list = []
    for i in range(1, len(snr)):
        s = snr[i]
        s_b = snr[i-1]
        # calculate mean flux boosting
        snr_select = np.bitwise_and((snr_input<s), (snr_input>s_b))
        aperture_boosting = flux_input_aperture[snr_select] / flux_input[snr_select]
        gaussian_boosting = flux_input_gaussian[snr_select] / flux_input[snr_select]
        if mode == 'mean':
            aperture_mean.append([np.mean(aperture_boosting), np.std(aperture_boosting)])
            gaussian_mean.append([np.mean(gaussian_boosting), np.std(gaussian_boosting)])
        elif mode == 'median':
            aperture_mean.append([np.median(aperture_boosting), np.std(aperture_boosting)])
            gaussian_mean.append([np.median(gaussian_boosting), np.std(gaussian_boosting)])



        # calculate the completeness
        snr_select1 = np.bitwise_and((detection_input_array[:, 0]<s), (detection_input_array[:, 0]>s_b))
        n_input = np.sum(snr_select1)
        n_found = np.sum(np.array(detection_input_array[:,1][snr_select1]))
        completeness_list.append(1.0*n_found/n_input)
        
        # calculate the fake detaction rate
        snr_select2 = np.bitwise_and((detection_found_array[:, 0]<s), (detection_found_array[:, 0]>s_b))
        n_found2 = np.sum(snr_select2)
        n_fake = n_found2 - np.sum(np.array(detection_found_array[:,1][snr_select2]))
        fake_rate_list.append(1.0*n_fake/n_found2)
        # total_found_list.append(n_found2)
        # total_fake_list.append(n_fake)


    if plot:
        fig = plt.figure(figsize=(12, 3))
        ax = fig.add_subplot(1,3,1)
        ax.set_xlabel('SNR')
        ax.set_ylabel(r'$S_{\rm out}/S_{\rm in}$')
        ax.plot(snr_input, flux_input_aperture/flux_input, 'k.', label='aperture')
        ax.plot(snr_input, flux_input_gaussian/flux_input, 'r.', label='gaussian')
        aperture_mean = np.array(aperture_mean)
        gaussian_mean = np.array(gaussian_mean)
        ax.plot(snr_mid, aperture_mean[:,0], 'ko', label='aperture')
        ax.errorbar(snr_mid, aperture_mean[:,0], yerr=aperture_mean[:,1], color='k', lw=2, capsize=5, elinewidth=2, markeredgewidth=2, alpha=0.8)
        ax.plot(snr_mid, gaussian_mean[:,0], 'ro', label='gaussian')
        ax.errorbar(snr_mid, gaussian_mean[:,0], yerr=gaussian_mean[:,1], color='r', lw=2, capsize=5, elinewidth=2, markeredgewidth=2, alpha=0.8)
        # for i in range(len(flux_input_list)):
            # if len(snr_input_list[i]>0):
                # # print(snr_inputfound_list[i])
                # # print(flux_input_list[i])
                # # print(flux_inputfound_list[i])
                # ax.plot(snr_input_list[i], flux_input_autolist[i][:,0]/flux_input_list[i], 'k.', label='aperture')
                # ax.plot(snr_input_list[i], flux_input_autolist[i][:,1]/flux_input_list[i], 'r.', label='gaussian')
        
        ax = fig.add_subplot(1,3,2)
        # print('snr', snr)
        # print('completeness_list', completeness_list)
        ax.plot(0.5*(snr[1:]+snr[:-1]), completeness_list, 'o')
        # ax.plot(np.array(flux_peak_list)/rms, completeness_list, 'o')
        ax.set_xlabel('SNR')
        ax.set_ylabel(r'Completeness')
        # ax.set_xlim((0., 8))
        ax.set_ylim((-0.1, 1.2))

        ax = fig.add_subplot(1,3,3)
        ax.plot(0.5*(snr[1:]+snr[:-1]), fake_rate_list, 'o')
        # ax.plot(0.5*(snr[1:]+snr[:-1]), total_found_list, 'ro')
        # ax.plot(0.5*(snr[1:]+snr[:-1]), total_fake_list, 'bo')
        ax.set_xlabel('SNR')
        ax.set_ylabel(r'Fake percentage')
        # ax.set_xlim((0., 8))
        # ax.set_ylim((-0.1, 1.2))

        plt.show()
    return snr_mid, [aperture_mean, gaussian_mean], completeness_list, fake_rate_list

def image_sim(image, outdir='./',):
    """calculate the completeness for given fitsfile
    
    image: fits image
    """

    gen_sim_images()
    calculate_sim_images()
    save(fluxbooting)
    save(completeness)
