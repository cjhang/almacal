# functions related with almacal simulation
import numpy as np
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

def make_random_source(direction, freq=None, n=1, radius=1, prune=False,
        prune_threshold=2, debug=False, savefile=None, clname=None,
        flux=None, fluxunit='mJy', known_file=None):
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
    if isinstance(flux, (list, tuple)):
        flux_input = np.random.uniform(flux[0], flux[1], n)
    elif isinstance(flux, (int, float)):
        flux_input = np.full(n, flux)
    if debug:
        print('fluxunit', fluxunit)

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
        distance = 2*u.arcsec
        print('distance', distance)
        selected_after_known_ra = []
        selected_after_known_dec = []
        known_data = np.loadtxt(known_file, skiprows=1)
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

def add_random_sources(vis, n=5, radius=10, outdir='./', make_image=True, 
        basename=None, debug=False, flux=None, known_file=None, 
        uvtaper_scale=[1.0, 2.0]):
    """
    radius : in arcsec
    """
    
    # clear the working place
    if basename is None:
        basename = os.path.basename(vis) + '.tmp'
    vis_testfile = os.path.join(outdir, basename)
    rmtables(vis_testfile)

    if not os.path.isdir(outdir):
        os.system('mkdir -p {}'.format(outdir))
   
    # make a copy of the original file
    split(vis=vis, outputvis=vis_testfile, datacolumn='data') 
    
    md = msmdtool()
    if not md.open(vis_testfile):
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
    mycomplist = make_random_source(mydirection, freq=myfreq, n=n, radius=radius, debug=debug, prune=True, flux=flux, 
                                    clname=clname_fullpath, savefile=savefile_fullpath, known_file=known_file)
    ft(vis=vis_testfile, complist=mycomplist)
    uvsub(vis=vis_testfile, reverse=True)
    
    vis_testfile_new = vis_testfile+'.new'
    rmtables(vis_testfile_new)
    split(vis=vis_testfile, datacolumn='corrected', outputvis=vis_testfile_new)
    if make_image:
        make_cont_img(vis_testfile_new, outdir=outdir, clean=True, niter=1000, 
                      only_fits=True, uvtaper_scale=uvtaper_scale, pblimit=-0.01,
                      fov_scale=2.0,
                      basename=basename)

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

def make_gaussian_image(shape, sigma, area=1., offset=(0,0), theta=0):
    """make a gaussian image for testing

    theta: in rad, rotating the gaussian counterclock wise
    """
    image = np.zeros(shape, dtype=float)
    yidx, xidx = np.indices(shape)
    yrad, xrad = yidx-shape[0]/2., xidx-shape[1]/2.
    y = xrad*np.cos(theta) + yrad*np.sin(theta)
    x = yrad*np.cos(theta) - xrad*np.sin(theta)
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

def auto_photometry(image, bmaj=1, bmin=1, theta=0, debug=False, methods=['aperture','gaussian']):
    """automatically measure the flux with different methods
    """
    
    ysize, xsize = image.shape
    y_center, x_center = ysize/2., xsize/2.
    
    flux_list = []
    # Aperture Photometry
    if 'aperture' in methods:
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
        flux_list.append(flux_apers_stable)
    if 'gaussian' in methods:
        yidx, xidx = np.indices((ysize, xsize))
        yrad, xrad = yidx-ysize/2., xidx-xsize/2.

        p_init = models.Gaussian2D(amplitude=1, x_stddev=1.*bmaj, y_stddev=1.*bmin, theta=theta)
        fit_p = fitting.LevMarLSQFitter()
        p = fit_p(p_init, xrad, yrad, image)
        flux_fitted = 2*np.pi*p.x_stddev.value*p.y_stddev.value*p.amplitude.value

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
        flux_list.append(flux_fitted)

    return flux_list

def source_finder(fitsimage, sources_file=None, savefile=None, model_background=True, 
                  threshold=5.0, debug=False, algorithm='find_peak', return_image=False,
                  filter_size=None, box_size=None, methods=['aperture', 'gaussian'],
                  known_file=None):
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
        known_data = np.loadtxt(known_file, skiprows=1)
        if len(known_data.shape) == 1:
            known_sources = SkyCoord(ra=known_data[0], dec=known_data[1], unit='arcsec')
            known_sources_pixel = skycoord_to_pixel(known_sources, wcs)
            known_aper = RectangularAperture(known_sources_pixel, 2.*a, 2.*a, theta=0)
        elif len(known_data.shape) > 1:
            known_sources = SkyCoord(ra=known_data[:,0], dec=known_data[:,1], unit='arcsec')
            known_sources_pixel = skycoord_to_pixel(known_sources, wcs)
            known_aper = RectangularAperture(list(zip(*known_sources_pixel)), 2.*a, 2.*a, theta=0)
        else:
            raise ValueError('Unsupported file: {}'.format(known_file))
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
    
    # find stars
    if algorithm == 'DAOStarFinder': # DAOStarFinder
        daofind = DAOStarFinder(fwhm=fwhm_pixel, threshold=threshold*std, ratio=ratio, 
                                theta=theta+90, sigma_radius=1, sharphi=0.7, sharplo=0.2,
                                mask=known_mask)  
        sources_found = daofind(data_masked - background)
        source_found_peak = sources_found['peak']
        if debug:
            print(sources_found)
        if len(sources_found) < 1:
            print("No point source!")
            return 0
        sources_found_x, sources_found_y = sources_found['xcentroid'], sources_found['ycentroid']
        
    elif algorithm == 'find_peak': # find_peak
        sources_found = find_peaks(data_masked - background, threshold=threshold*std, box_size=fwhm_pixel,
                                   mask=known_mask)
        if debug:
            print(sources_found)
        if len(sources_found) < 1:
            print("No point source!")
            return 0
        sources_found_x, sources_found_y = sources_found['x_peak'].data, sources_found['y_peak'].data
        source_found_peak = sources_found['peak_value']

    else:
        raise ValueError("Unsurport algorithm: {}!".format(algorithm))
    sources_found_center = list(zip(sources_found_x, sources_found_y))
    sources_found_coords = pixel_to_skycoord(sources_found_x, sources_found_y, wcs)

    # aperture photometry based on source finding coordinates
    ## simple aperture flux
    if True:
        aper_found = EllipticalAperture(sources_found_center, 1*a, 1*b, theta=theta+90/180*np.pi)
        phot_table_found = aperture_photometry(data_masked, aper_found)
        flux_aper_found = (phot_table_found['aperture_sum'] / beamsize * 1000).tolist() # convert mJy
    # automatically aperture photometry
    if True:
        flux_auto = []
        segments = RectangularAperture(sources_found_center, 2.*a, 2.*a, theta=0)
        segments_mask = segments.to_mask(method='center')
        for s in segments_mask:
            flux_list = auto_photometry(s.cutout(data_masked), bmaj=b, bmin=a, 
                                        theta=theta/180*np.pi, debug=debug, methods=methods)
            flux_auto.append(np.array(flux_list) / beamsize * 1000) #from Jy to mJy
        # return segments_mask
    if debug:
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

        idx_found, idx_input, d2d, d3d = sources_input_coords.search_around_sky(
                sources_found_coords, fwhm)

        # print(idx_found)
        # print(idx_input)
    
        # print(flux_input[idx_input])
        # print(np.array(flux_auto)[idx_found])

        sources_input_found = [flux_input[idx_input], np.array(flux_auto)[idx_found], np.array(source_found_peak[idx_found])]
        # aperture photometry based on input coordinates
        # for s in source_input_pixels:
            # pos_list.append([s['xcentroid'], s['ycentroid']])
        sources_input_center = zip(*sources_input_pixels)
        aper_input = EllipticalAperture(sources_input_center, 3*b, 3*a, 
                                  theta=theta/180*np.pi)
        phot_table_input = aperture_photometry(data_masked, aper_input)
        flux_aper_input = (phot_table_input['aperture_sum'] / beamsize).tolist() # convert mJy
        # automatically aperture photometry
        if True:
            flux_input_auto = []
            segments = RectangularAperture(sources_input_center, 2.*a, 2.*a, theta=0)
            segments_mask = segments.to_mask(method='center')
            for s in segments_mask:
                flux_list = auto_photometry(s.cutout(data_masked), bmaj=b, bmin=a, 
                                            theta=theta/180*np.pi, debug=debug, methods=methods)
                flux_input_auto.append(np.array(flux_list) / beamsize) #already mJy
 
        if debug:
            print(flux_input)
            print(flux_aper_input)
            print(flux_auto)

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
        for i,e in enumerate(sources_found_center):
            ellipse = patches.Ellipse(e, width=3*b, height=3*a, angle=theta, facecolor=None, fill=False, edgecolor='red', alpha=0.4, linewidth=4)
            ax.add_patch(ellipse)
            if debug:
                ax.text(e[0], e[1], flux_aper_found[i])

        plt.show()

    if savefile:
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
        return flux_input, flux_auto, sources_input_found
    return flux_auto


