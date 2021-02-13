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

def make_random_source(direction, freq=None, n=1, radius=1, prune=False,
        prune_threshold=0.5, debug=False, savefile=None, clname=None,
        flux_range=[0, 1], fluxunit='mJy'):
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
    flux_random = np.random.uniform(flux_range[0], flux_range[1], n)
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
                continue
            selected_ra.append(delta_ra[i]) 
            selected_dec.append(delta_dec[i])
        delta_ra = np.array(selected_ra)
        delta_dec = np.array(selected_dec)
    ra_random = delta_ra*u.arcsec + skycoord.ra
    dec_random = delta_dec*u.arcsec+ skycoord.dec
    

    if savefile:
        with open(savefile, 'w+') as sfile:
            sfile.write('# ra[arcsec]  dec[arcsec]  flux[{}]\n'.format(fluxunit))
            for ra, dec, flux in zip(ra_random, dec_random, flux_random):
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
        for d,f in zip(direction_list, flux_random):
            cl.addcomponent(dir=d, flux=f, fluxunit=fluxunit, 
                            freq=freq, shape='point')
        cl.rename(clname)
        cl.done()
        return clname
    else:
        return [ra_random, dec_random, flux_random]

    # skycoord_list = SkyCoord(ra=ra_random, dec=dec_random)
    # f = lambda x: ('J2000 '+x.to_string('hmsdms')).encode('utf-8')
    # return list(map(f, skycoord_list))
   
    # The old method
    # direction_list = []
    # while len(direction_list) < n:
        # # random value for ra and dec
        # # two_delta = (np.random.random_sample(2) - 0.5) * 2 * radius        
        # delta1 = (np.random.random_sample() - 0.5) * 2 #* radius        
        # delta2 = (np.random.random_sample() - 0.5) * 2 #* radius        
        # two_delta = np.array([delta1, delta2]) * radius
        # # print(two_delta[0], two_delta[1])
        # if (two_delta[0]**2 + two_delta[1]**2) >= radius**2:
            # continue
        # ra_random = two_delta[0] + skycoord.ra
        # dec_random = two_delta[1] + skycoord.dec
        # skycoord_tmp = SkyCoord(ra=ra_random, dec=dec_random)
        # direction_list.append('J2000 ' + skycoord_tmp.to_string('hmsdms'))
   
    # return direction_list

def add_raondom_source(vis, n=5, radius=10, outdir='./', make_image=False, 
        basename=None, debug=False):
    """
    radius : in arcsec
    """
    
    # clear the working place
    vis_testfile = os.path.join(outdir, basename+'.ms')
    rmtables(vis_testfile)
   
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
    mycomplist = make_random_source(mydirection, freq=myfreq, n=n, radius=radius, debug=debug, prune=True, clname=clname_fullpath, savefile=savefile_fullpath)
    ft(vis=vis_testfile, complist=mycomplist)
    uvsub(vis=vis_testfile, reverse=True)
    
    vis_testfile_new = vis_testfile+'.new'
    rmtables(vis_testfile_new)
    split(vis=vis_testfile, datacolumn='corrected', outputvis=vis_testfile_new)
    if make_image:
        make_cont_img(vis_testfile_new, outdir=outdir, clean=True, 
                      only_fits=True, basename=basename)

def make_gaussian_image(shape, sigma, offset=(0,0), theta=0):
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
    flux = np.exp(-(x-x_offset)**2/2./xsigma**2 - (y-y_offset)**2/2./ysigma**2) / (2*np.pi*xsigma*ysigma)

    return flux

def auto_photometry(image, a=1, b2a=1, theta=0, debug=False, ):
    """automatically measure the flux with different methods
    """
    
    ysize, xsize = image.shape
    y_center, x_center = ysize/2., xsize/2.
    
    # Aperture Photometry
    if method == 'aperture':
        radii = np.arange(1., ysize/2./np.cos(theta), 1/np.cos(theta))

        apertures = [EllipticalAperture((x_center, y_center), 1.*a, 1.*a*b2a, theta) for a in radii]
        extract_area = lambda x: x.area()
        area_apers = np.array(list(map(extract_area, apertures)))
        
        phot_table = aperture_photometry(image, apertures)
        extract_qtable = lambda x: x.data[0]
        flux_apers = np.array(list(map(extract_qtable, phot_table.columns[3:].values())))

        slopes = np.ones_like(area_apers)
        slopes[1:] = np.diff(flux_apers) / np.diff(area_apers)
        slope_selection = slopes < 0.001
        flux_apers_stable = flux_apers[slope_selection][0]
        print(flux_apers_stable)
        if debug:
            print('radii', radii)
            print(area_apers)
            print(flux_apers)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.imshow(image, origin='lower')
            for ap in apertures:
                ap.plot(color='gray', lw=2)
        return flux_apers_stable
    elif method == 'gaussian':
        pass
        return 0



def source_finder(fitsimage, sources_file=None, savefile=None, model_background=True, 
                  threshold=5.0, debug=False, algorithm='find_peak', return_image=False,
                  filter_size=None, box_size=None):
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
    beamsize = np.pi*a*b/4*np.log(2)
    
    if debug:
        print('image shape', ny, nx)
        print("sigma stat:", mean, median, std)
        print('fwhm_pixel:', fwhm_pixel)
        print('a, b', a, b)

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
                                theta=theta+90, sigma_radius=1, sharphi=0.7, sharplo=0.2)  
        sources_found = daofind(data_masked - background)
        if debug:
            print(sources_found)
        if len(sources_found) < 1:
            print("No point source!")
            return 0
        sources_found_x, sources_found_y = sources_found['xcentroid'], sources_found['ycentroid']
        
    elif algorithm == 'find_peak': # find_peak
        sources_found = find_peaks(data_masked - background, threshold=threshold*std)
        if debug:
            print(sources_found)
        if len(sources_found) < 1:
            print("No point source!")
            return 0
        sources_found_x, sources_found_y = sources_found['x_peak'].data, sources_found['y_peak'].data

    else:
        raise ValueError("Unsurport algorithm: {}!".format(algorithm))
    sources_found_center = list(zip(sources_found_x, sources_found_y))
    sources_found_coords = pixel_to_skycoord(sources_found_x, sources_found_y, wcs)
    # aperture photometry based on source finding coordinates
    aper_found = EllipticalAperture(sources_found_center, 1*a, 1*b, theta=theta/180*np.pi)
    phot_table_found = aperture_photometry(data_masked, aper_found)
    flux_aper_found = (phot_table_found['aperture_sum'] * beamsize).tolist() # convert mJy
    if True:
        segments = RectangularAperture(sources_found_center, 2.*a, 2.*a, theta=0)
        print(segments)
        segments_mask = segments.to_mask(method='center')
        print(segments_mask)
        return segments_mask
    if debug:
        print("sources_found_center", sources_found_center)
        print("aper_found.positions", aper_found.positions)
        print('flux in aperture', flux_aper_found)


    if sources_file:
        sources_input = np.loadtxt(sources_file, skiprows=1)
        flux_input = sources_input[:,-1]
        sources_input_coords = SkyCoord(ra=sources_input[:,0]*u.arcsec, 
                                        dec=sources_input[:,1]*u.arcsec)

        sources_input_pixels = skycoord_to_pixel(sources_input_coords, wcs)

        idx_found, idx_input, d2d, d3d = sources_input_coords.search_around_sky(
                sources_found_coords, fwhm)
    
        # aperture photometry based on input coordinates
        # for s in source_input_pixels:
            # pos_list.append([s['xcentroid'], s['ycentroid']])
        aperture_input_centers = zip(*sources_input_pixels)
        aper_input = EllipticalAperture(zip(*sources_input_pixels), 3*a, 3*b, 
                                  theta=theta/180*np.pi)
        phot_table_input = aperture_photometry(data_masked, aper_input)
        flux_aper_input = (phot_table_input['aperture_sum'] * beamsize).tolist() # convert mJy

        if debug:
            print(flux_input)
            print(flux_aper_input)

    if debug:
        # visualize the results
        fig= plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(data_masked, origin='lower')
        
        if sources_file:
            # show the input sources
            for i,e in enumerate(aperture_input_centers):
                ellipse = patches.Ellipse(e, width=3*b, height=3*a, angle=theta, facecolor=None, fill=False, edgecolor='orange', alpha=0.8, linewidth=1)
                ax.add_patch(ellipse)
                if debug:
                    ax.text(e[0], e[1], flux_input[i])

        # show the sources found
        for i,e in enumerate(sources_found_center):
            print(e)
            ellipse = patches.Ellipse(e, width=3*b, height=3*a, angle=theta, facecolor=None, fill=False, edgecolor='red', alpha=0.4, linewidth=4)
            ax.add_patch(ellipse)
            if debug:
                ax.text(e[0], e[1], flux_aper_found[i])

        plt.show()

    if savefile:
        with open(savefile, 'w+') as sfile:
            sfile.write('# ra[arcsec]  dec[arcsec]  flux[mJy]\n')
            for ra, dec, flux in zip(sources_found_coords.ra, sources_found_coords.dec, flux_aper_found):
                sfile.write('{:.5f} {:.6f} {:.4f}\n'.format(ra.to(u.arcsec).value,
                    dec.to(u.arcsec).value, flux))
    if return_image:
        return data_masked.data, background

    return sources_found_coords

