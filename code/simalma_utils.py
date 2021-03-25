# functions related with almacal simulation
import numpy as np
from scipy import signal, ndimage
from astropy import units as u
from astropy.coordinates import SkyCoord

from astropy.io import fits
from astropy.stats import sigma_clipped_stats, sigma_clip, SigmaClip
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from photutils import (DAOStarFinder, EllipticalAperture, aperture_photometry, 
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

def make_random_source(direction, freq=None, n=1, radius=1, prune=False,
        prune_threshold=3, debug=False, savefile=None, clname=None,
        flux=None, fluxunit='mJy', known_file=None, known_data=None,
        sampler=np.random.power, scaler=1):
    """This function used to add random source around a given direction
        
    This is a proximate solution, which treat the field of view as a plane.
    Then calculate the longtitude and latitue using sin and cos functions.

    Args:
    direction: the pointing direction
    n: the number of points
    radius: the largest radius that the points can be put, in arcsec
    prune: remove close pairs, make it easier for pointsource testing
    savefile: save the coordination files
    clname: the component list filename

    direction: 'J2000 13h03m49.215900s -55d40m31.60870s'
   """
    # generate a random radius
    skycoord = SkyCoord(direction[5:])
    theta = 2*np.pi*np.random.uniform(0, 1, n)
    rho = radius * np.sqrt(np.random.uniform(0, 1, n))
    print('flux', flux)
    if isinstance(flux, (list, tuple, np.ndarray)):
        flux_sampling = sampler(scaler, n)
        flux_input = flux_sampling * np.diff(flux) + flux[0]
    elif isinstance(flux, (int, float)):
        flux_input = np.full(n, flux)
    if debug:
        print('fluxunit', fluxunit)

    print('flux_input',flux_input)
    delta_ra = rho * np.cos(theta)/(np.cos(skycoord.dec).value)
    delta_dec = rho * np.sin(theta)

    if prune:
        selected_ra = [delta_ra[0],] 
        selected_dec = [delta_dec[0],]
        for i in range(1, n):
            coord_i = SkyCoord(ra=delta_ra[i]*u.arcsec, dec=delta_dec[i]*u.arcsec)
            coords_selected = SkyCoord(ra=np.array(selected_ra)*u.arcsec, dec=np.array(selected_dec)*u.arcsec)
            separation_with_select = coord_i.separation(coords_selected)
            if np.sum(separation_with_select<prune_threshold*u.arcsec)>=1:
                print("pruning one point...")
                continue
            selected_ra.append(delta_ra[i]) 
            selected_dec.append(delta_dec[i])
        delta_ra = np.array(selected_ra)
        delta_dec = np.array(selected_dec)
    ra_random = delta_ra*u.arcsec + skycoord.ra
    dec_random = delta_dec*u.arcsec+ skycoord.dec
   
    if known_file:
        try:
            known_data = np.loadtxt(known_file, skiprows=1)
        except:
            print('No known file can be open.')
            known_data = None
    if known_data is not None:
        distance = prune_threshold*u.arcsec
        print('distance', distance)
        selected_after_known_ra = []
        selected_after_known_dec = []
        if len(known_data.shape) == 1:
            known_sources = SkyCoord(ra=known_data[0], dec=known_data[1], unit='arcsec')
        elif len(known_data.shape) > 1:
            known_sources = SkyCoord(ra=known_data[:,0], dec=known_data[:,1], unit='arcsec')
        else:
            raise ValueError('Unsupported file: {}'.format(known_file))
        # print('known_sources ra:', known_sources.ra.value)
        # print('known_sources dec:', known_sources.dec.value)
        for ra,dec in zip(ra_random, dec_random):
            if np.sum(np.abs(ra - known_sources.ra) < distance)>0 and np.sum(np.abs(dec - known_sources.dec) < distance)>0:
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
        skycoord_list = SkyCoord(ra=ra_random, dec=dec_random)
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

def add_random_sources(vis=None, fitsimage=None, n=5, radius=10, outdir='./', make_image=True, 
        basename=None, debug=False, flux=None, known_file=None, uvtaper_scale=None,
        inverse_image=False, **kwargs):
    """
    radius : in arcsec

    Notes:
        1. To make the reading and writing more efficiently, the model will be added directly
           into original measurement. It is suggested to make a copy before calling this func
    """
    
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
    if inverse_image:
        flux = -1.0 * flux
   
    if vis:
        md = msmdtool()
        if not md.open(vis):
            raise ValueError("Failed to open {}".format(vis))
        phasecenter = md.phasecenter()
        mydirection = phasecenter['refer'] +' '+ SkyCoord(phasecenter['m0']['value'], phasecenter['m1']['value'], unit="rad").to_string('hmsdms')
        myfreq = "{:.2f}GHz".format(np.mean(read_spw(vis)))
        if debug:
            print(mydirection)
            print(myfreq)
        clname_fullpath = os.path.join(outdir, basename+'.cl')
        savefile_fullpath = os.path.join(outdir, basename+'.txt')
        # overwrite old files
        rmtables(clname_fullpath)
        os.system('rm -rf {}'.format(savefile_fullpath))
        # generate random sources
        mycomplist = make_random_source(mydirection, freq=myfreq, n=n, radius=radius, debug=debug, flux=flux, 
                                        clname=clname_fullpath, savefile=savefile_fullpath, known_file=known_file, **kwargs)
        ft(vis=vis, complist=mycomplist)
        uvsub(vis=vis, reverse=True)
        
        if make_image:
            suffix = '.cont.auto'
            make_cont_img(vis, outdir=outdir, clean=True, niter=1000, 
                          only_fits=True, uvtaper_scale=uvtaper_scale, pblimit=-0.01,
                          fov_scale=2.0, datacolumn='corrected', usemask='user',
                          basename=basename, suffix=suffix)
            if inverse_image:
                # tb = tbtool()
                # tb.open(os.path.join(outdir, basename+suffix+'.fits'), nomodify=False)
                # tb.put('map', -1.*tb.getcol('map'))
                # tb.flush()
                # tb.close()
                with fits.open(os.path.join(outdir, basename+suffix+'.fits')) as hdu:
                    hdu[0].data = -1.0 * hdu[0].data
                    hdu.writeto(os.path.join(outdir ,basename+suffix+'.fits'), overwrite=True)

            if uvtaper_scale:
                tb = tbtool()
                for taper in uvtaper_scale:
                    tb.open(os.path.join(outdir, basename+suffix+'uvtaper{}.fits'.format(taper)), nomodify=False)
                    tb.put('map', -1.*tb.getcol('map'))
                    tb.flush()
                    tb.close()
        # clean up temperary files
        rmtables(clname_fullpath)
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
        mean, median, std = sigma_clipped_stats(data_masked, sigma=10.0)  

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
        refer = 'J'+str(int(header['EQUINOX']))
        mydirection = refer +' '+ SkyCoord(header['CRVAL1'], header['CRVAL2'], unit="deg").to_string('hmsdms')
        myfreq = "{:.2f}GHz".format(header['CRVAL3']/1e9)
        if debug:
            print(mydirection)
            print(myfreq)
        savefile_fullpath = os.path.join(outdir, basename+'.txt')
        mycomplist = make_random_source(mydirection, freq=myfreq, n=n, radius=radius, debug=debug, flux=flux, 
                                        savefile=savefile_fullpath, known_file=known_file)
        # print('mycomplist', mycomplist)
        mycomplist_pixels = []
        blank_image = np.zeros((ny, nx))
        for ra, dec, flux in zip(*mycomplist):
            ra_pix, dec_pix = map(np.around, skycoord_to_pixel(SkyCoord(ra, dec), wcs))
            # print(ra_pix, dec_pix, flux)
            mycomplist_pixels.append([int(ra_pix), int(dec_pix), flux])
            blank_image[int(dec_pix), int(ra_pix)] = flux
            
        fake_image = convolve(blank_image, kernel)/1000*beamsize + data_masked # convert mJy into Jy before added into image
        hdu[0].data = fake_image.filled(np.nan)
        if inverse_image:
            hdu[0].data = -1.0*hdu[0].data
        hdu.writeto(os.path.join(outdir ,basename+'.fits'), overwrite=True)
        return fake_image 

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

def read_source(sourcefile):
    pass

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

def auto_photometry(image, bmaj=1, bmin=1, theta=0, beamsize=None, debug=False, methods=['aperture','gaussian', 'peak']):
    """automatically measure the flux with different methods
    
    aperture_size = 1.0 # in units of beam size
    """
    
    ysize, xsize = image.shape
    y_center, x_center = ysize/2., xsize/2.
    flux_list = []
    if beamsize is None:
        print("Warning: Unit the sum and maximum, without dividing the beam size.")
    # Aperture Photometry
    if 'aperture' in methods:
        aperture_size=1.0
        aperture_correction = 1.068
        sigma2FWHM = 2.35482
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
            ax.imshow(image, origin='lower', interpolation='none')
            aperture.plot(color='gray', lw=2)
        if beamsize:
            flux_list.append(1.0*flux_in_aper/beamsize)
        else:
            flux_list.append(flux_in_aper)
    if 'gaussian' in methods:
        yidx, xidx = np.indices((ysize, xsize))
        yrad, xrad = yidx-ysize/2., xidx-xsize/2.

        image_norm = image / np.max(image)
        p_init = models.Gaussian2D(amplitude=1, x_stddev=bmaj/sigma2FWHM, y_stddev=bmin/sigma2FWHM, theta=theta, 
                                   bounds={"x_mean":(-2.,2.), "y_mean":(-2.,2.), 
                                           "x_stddev":(0., xsize), "y_stddev":(0., ysize)})
        fit_p = fitting.LevMarLSQFitter()
        p = fit_p(p_init, xrad, yrad, image_norm)
        flux_fitted = 2*np.max(image)*np.pi*p.x_stddev.value*p.y_stddev.value*p.amplitude.value

        if debug:
            print("Initial guess:", p_init)
            print("Fitting:", p)
            print("Flux:", flux_fitted)
            fig, ax = plt.subplots(1, 3, figsize=(8, 3.5))
            ax[0].imshow(image, origin='lower', interpolation='none')
            gaussian_init = patches.Ellipse((ysize/2., xsize/2.), height=p_init.y_stddev.value, 
                    width=p_init.x_stddev.value, angle=p_init.theta.value/np.pi*180,
                    linewidth=2, facecolor=None, fill=None, edgecolor='orange', alpha=0.8)
            ax[0].add_patch(gaussian_init)
            ax[0].set_title("Data")
            ax[1].imshow(p(xrad, yrad), origin='lower', interpolation='none')
            ax[1].set_title("Model")
            ax[2].imshow(image - p(xrad, yrad), origin='lower', interpolation='none')
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
            ax.imshow(image, origin='lower', interpolation='none')
            for ap in apertures:
                ap.plot(color='gray', lw=2)
        if beamsize:
            flux_list.append(1.0*flux_apers_stable/beamsize)
        else:
            flux_list.append(flux_apers_stable)
    return flux_list

def gen_fake_images(vis=None, imagefile=None, known_file=None, n=20, repeat=10, 
                    snr=np.arange(1,20,0.1), fov_scale=1.5, outdir='./', basename=None,
                    uvtaper_scale=None, mode='image', inverse_image=False, 
                    debug=False, **kwargs):
    """generate the fake images with man-made sources

    Args:
     mode: image or uv
    """
    # make a copy of the original file
    if basename is None:
        basename = os.path.basename(vis)
    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
    vistmp = os.path.join(outdir, basename+'.tmp')
    split(vis=vis, outputvis=vistmp, datacolumn='data')
    # read information from vis 
    spw_specrange = read_spw(vistmp)
    freq_mean = np.mean(spw_specrange) # in GHz
    tb.open(vistmp + '/ANTENNA')
    antenna_diameter_list = tb.getcol('DISH_DIAMETER')
    tb.close()
    antenna_diameter = np.max(antenna_diameter_list) * u.m
    wavelength = const.c / (freq_mean * u.GHz) # in um
    fov = (fov_scale * 1.22 * wavelength / antenna_diameter * 206265).decompose()
    if debug:
        print('fov', fov)
        print('radius', 0.5*fov)
    # read information from image files
    imagefile_base, imagefile_extension = os.path.splitext(imagefile)
    print(imagefile_base, imagefile_extension)
    if False:#imagefile_extension == 'fits':
        fitsimage = imagefile
        hdu = fits.open(fitsimage)
        header = hdu[0].header
        data = np.ma.masked_invalid(hdu[0].data)
        pixel_scale = 1/np.abs(header['CDELT1'])
        fwhm = header['BMAJ']*3600*u.arcsec
        fwhm_pixel = header['BMAJ']*pixel_scale
        a, b = header['BMAJ']*pixel_scale, header['BMIN']*pixel_scale
        ratio = header['BMIN'] / header['BMAJ']
        theta = header['BPA']
        beamsize = np.pi*a*b/(4*np.log(2))
        gaussian_norm = 2*np.pi*a*b/2.35482**2
        sensitivity = np.ma.std(data) * 1000 # in mJy/pixel
        # flux_base = rms
    else:
        im_head = imhead(imagefile)
        im_beam = im_head['restoringbeam']
        im_incr = im_head['incr']
        im_info = imstat(imagefile)
        # beamsize = np.pi*a*b/(4*np.log(2))
        beamsize = np.pi/(4*np.log(2))* im_beam['major']['value'] * im_beam['minor']['value'] / (im_incr[0]/np.pi*180*3600)**2
        gaussian_norm = 2.0*np.pi*im_beam['major']['value'] * im_beam['minor']['value'] / 2.35482**2
        sensitivity = im_info['rms'] * 1000 # in mJy/pixel
        # flux_base = rms[0] #* 2*np.pi # convert the peak flux density into brightness
        if mode == 'image':
            fitsimage = imagefile+'.fits'
            exportfits(imagename=imagefile, fitsimage=imagefile+'.fits', overwrite=True)
    # print('sensitivity', sensitivity)
    # print('gaussian_norm', gaussian_norm)
    # print('beamsize', beamsize)

    # if snr_mode == 'peak':
    flux_base = sensitivity #* gaussian_norm
    if len(snr)>3:
        # the old method, which will be deprecated in the near future
        for s in snr:
            # if s in snr_old:
                # continue
            print(">>>>>>>\n>> snr={}\n>>>>>>>>".format(s))
            for i in range(repeat):
                basename_new = basename+'.snr{}.run{}'.format(s, i)
                if mode == 'uv':
                    add_random_sources(vis=vistmp, n=20, radius=0.5*fov, outdir=outdir,
                            uvtaper_scale=uvtaper_scale, basename=basename_new, 
                            flux=s*flux_base, known_file=known_file, 
                            inverse_image=inverse_image)
                elif mode == 'image':
                    add_random_sources(fitsimage=fitsimage, n=20, radius=0.5*fov, 
                            outdir=outdir,uvtaper_scale=uvtaper_scale, 
                            basename=basename_new, flux=s*flux_base, known_file=known_file, 
                            inverse_image=inverse_image)
        # clear temperary files
        rmtables(vistmp)
    elif len(snr)==2:
        flux = np.array(snr)*flux_base
        for i in range(repeat):
            basename_new = basename+'.run{}'.format(i)
            if mode == 'uv':
                add_random_sources(vis=vistmp, n=20, radius=0.5*fov, outdir=outdir,uvtaper_scale=uvtaper_scale, 
                                   basename=basename_new, flux=flux, known_file=known_file, 
                                   inverse_image=inverse_image, **kwargs)
            elif mode == 'image':
                add_random_sources(fitsimage=fitsimage, n=20, radius=0.5*fov, outdir=outdir,uvtaper_scale=uvtaper_scale, 
                                   basename=basename_new, flux=flux, known_file=known_file, inverse_image=inverse_image,
                                   **kwargs)
        # clear temperary files
        rmtables(vistmp)

def source_finder(fitsimage, sources_file=None, savefile=None, model_background=True, 
                  threshold=5.0, debug=False, algorithm='find_peak', return_image=False,
                  filter_size=None, box_size=None, methods=['aperture', 'gaussian','peak'],
                  subtract_background=False, known_file=None):
    """finding point source in the image
    """

    hdu = fits.open(fitsimage)
    header = hdu[0].header
    wcs = WCS(header)
    data = hdu[0].data
    ny, nx = data.shape[-2:]
    data_masked = np.ma.masked_invalid(data.reshape(ny, nx))
    mean, median, std = sigma_clipped_stats(data_masked, sigma=10.0)  

    # for DAOStarFinder, in pixel space
    pixel_scale = 1/np.abs(header['CDELT1'])
    fwhm = header['BMAJ']*3600*u.arcsec
    fwhm_pixel = header['BMAJ']*pixel_scale
    a, b = header['BMAJ']*pixel_scale, header['BMIN']*pixel_scale
    ratio = header['BMIN'] / header['BMAJ']
    theta = header['BPA']
    beamsize = np.pi*a*b/(4*np.log(2))

    # read known_sources
    known_mask = np.full((ny, nx), 0, dtype=bool)
    if known_file:
        try:
            known_data = np.loadtxt(known_file, skiprows=1)
        except:
            print('cannot open {}'.format(known_file))
            known_data = None
        if known_data is not None:
            if len(known_data.shape) == 1:
                known_sources = SkyCoord(ra=known_data[0], dec=known_data[1], unit='arcsec')
                known_sources_pixel = skycoord_to_pixel(known_sources, wcs)
                known_aper = RectangularAperture(known_sources_pixel, 2.*a, 2.*a, theta=0)
            elif len(known_data.shape) > 1:
                known_sources = SkyCoord(ra=known_data[:,0], dec=known_data[:,1], unit='arcsec')
                known_sources_pixel = skycoord_to_pixel(known_sources, wcs)
                known_aper = RectangularAperture(list(zip(*known_sources_pixel)), 2.*a, 2.*a, theta=0)
            # known_sources = SkyCoord(ra=known_data[:,0], dec=known_data[:,1], unit='arcsec')
            # known_sources_pixel = skycoord_to_pixel(known_sources, wcs)
            # known_aper = RectangularAperture(list(zip(*known_sources_pixel)), 2.*a, 2.*a, theta=0)
            known_aper_mask = known_aper.to_mask(method='center')
            for m in known_aper_mask:
                known_mask = np.bitwise_or(known_mask, m.to_image((ny,nx)).astype(bool))
            if False:
                plt.figure()
                plt.imshow(known_mask,origin='lower')
                plt.show()
    
    if debug:
        print('image shape', ny, nx)
        print("sigma stat:", mean, median, std)
        print('fwhm_pixel:', fwhm_pixel)
        print('a, b', a, b)
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
        daofind = DAOStarFinder(fwhm=fwhm_pixel, threshold=threshold*std*0.5, ratio=ratio, 
                                theta=theta+90, sigma_radius=1, sharphi=0.7, sharplo=0.2,)  
        data_masked[known_mask] = std
        sources_found = daofind(data_masked_sub)#, mask=known_mask)
        if debug:
            print(sources_found)
        if len(sources_found) < 1:
            print("No point source!")
            # return 0
            sources_found = None
        else:
            sources_found_x, sources_found_y = sources_found['xcentroid'], sources_found['ycentroid']
            source_found_peak = sources_found['peak']
            peak_select = source_found_peak > threshold #two stage source finding
            source_found_x = sources_found_x[peak_select]
            source_found_y = sources_found_y[peak_select]
            source_found_peak = source_found_peak[peak_select]
        
    elif algorithm == 'find_peak': # find_peak
        sources_found = find_peaks(data_masked_sub, threshold=threshold*std, 
                box_size=fwhm_pixel, mask=known_mask)
        if debug:
            print(sources_found)
        if len(sources_found) < 1:
            print("No point source!")
            # return 0
            sources_found = None
        else:
            sources_found_x = sources_found['x_peak'].data 
            sources_found_y = sources_found['y_peak'].data
            source_found_peak = sources_found['peak_value']

    else:
        raise ValueError("Unsurport algorithm: {}!".format(algorithm))
   

    flux_auto = []
    if sources_found is not None:
        # print('source_found_peak',source_found_peak)
        sources_found_center = list(zip(sources_found_x, sources_found_y))
        sources_found_coords = pixel_to_skycoord(sources_found_x, sources_found_y, wcs)

        # aperture photometry based on source finding coordinates
        ## simple aperture flux
        aper_found = EllipticalAperture(sources_found_center, 1*a, 1*b, theta=theta+90/180*np.pi)
        phot_table_found = aperture_photometry(data_masked, aper_found)
        flux_aper_found = (phot_table_found['aperture_sum'] / beamsize * 1000).tolist() # convert mJy

        # automatically aperture photometry
        segments = RectangularAperture(sources_found_center, 2.*a, 2.*a, theta=0)
        segments_mask = segments.to_mask(method='center')
        for s in segments_mask:
            if subtract_background:
                data_cutout = s.cutout(data_masked_sub)
            else:
                data_cutout = s.cutout(data_masked)
            flux_list = auto_photometry(data_cutout, bmaj=b, bmin=a, beamsize=beamsize,
                                        theta=theta/180*np.pi, debug=False, methods=methods)
            flux_auto.append(np.array(flux_list) * 1000) #from Jy to mJy
        # return segments_mask
    if debug:
        if sources_found:
            print("sources_found_center", sources_found_center)
            print("aper_found.positions", aper_found.positions)
            print('flux in aperture', flux_aper_found)
            print('auto_photometry:', flux_auto)


    if sources_file:
        sources_input = np.loadtxt(sources_file, skiprows=1)
        if len(sources_input.shape) == 1:
            sources_input_coords = SkyCoord(ra=sources_input[0], dec=sources_input[1], unit='arcsec')
            flux_input = sources_input[-1]
        elif len(sources_input.shape) > 1:
            sources_input_coords = SkyCoord(ra=sources_input[:,0], dec=sources_input[:,1], unit='arcsec')
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
                                        theta=theta/180*np.pi, debug=False, methods=methods)
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
            flux_auto = []
            idxs = [np.array([]), np.array([]), np.array([]), np.array([])]

        if debug:
            print(flux_input)
            print(flux_input_auto)

    if debug:
        # visualize the results
        fig= plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(data_masked, origin='lower')
        
        if sources_file:
            # show the input sources
            for i,e in enumerate(sources_input_center):
                ellipse = patches.Ellipse(e, width=3*b, height=3*a, angle=theta, facecolor=None, fill=False, edgecolor='orange', alpha=0.8, linewidth=1)
                ax.add_patch(ellipse)
                if debug:
                    ax.text(e[0], e[1], flux_input[i])

        # show the sources found
        if sources_found:
            for i,e in enumerate(sources_found_center):
                ellipse = patches.Ellipse(e, width=3*b, height=3*a, angle=theta, facecolor=None, fill=False, edgecolor='red', alpha=0.4, linewidth=4)
                ax.add_patch(ellipse)
                if debug:
                    ax.text(e[0], e[1], flux_aper_found[i])

        plt.show()

    if savefile and sources_found:
        with open(savefile, 'w+') as sfile:
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

    if sources_file:
        return np.array(flux_input), np.array(flux_input_auto), np.array(flux_auto), idxs    
    return flux_auto


