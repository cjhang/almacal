
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord


def make_random_source(direction, n=1, radius=1*u.arcsec):
    """This function used to add random source around a given direction
        
    This is a proximate solution, which treat the field of view as a plane.
    Then calculate the longtitude and latitue using sin and cos functions.

    direction: 'J2000 13h03m49.215900s -55d40m31.60870s'
   """
    # generate a random radius
    skycoord = SkyCoord(direction[5:])
    
    direction_list = []
    
    theta = 2*np.pi*np.random.uniform(0, 1, n)
    rho = radius * np.sqrt(np.random.uniform(0, 1, n))

    ra_random = rho * np.cos(theta) + skycoord.ra
    dec_random = rho * np.sin(theta) + skycoord.dec

    # f = lambda x: ('J2000 '+x.to_string('hmsdms')).encode('utf-8')
    # return list(map(f, skycoord_list))
   
    # The old method
    while len(direction_list) < n:
        # random value for ra and dec
        # two_delta = (np.random.random_sample(2) - 0.5) * 2 * radius        
        delta1 = (np.random.random_sample() - 0.5) * 2 #* radius        
        delta2 = (np.random.random_sample() - 0.5) * 2 #* radius        
        two_delta = np.array([delta1, delta2]) * radius
        # print(two_delta[0], two_delta[1])
        if (two_delta[0]**2 + two_delta[1]**2) >= radius**2:
            continue
        ra_random = two_delta[0] + skycoord.ra
        dec_random = two_delta[1] + skycoord.dec
        skycoord_tmp = SkyCoord(ra=ra_random, dec=dec_random)
        direction_list.append('J2000 ' + skycoord_tmp.to_string('hmsdms'))
   
    return direction_list


def add_raondom_source(vis, n=5, tmpdir='./tmp', make_image=False):
    """
    """
    md = msmdtool()
    if not md.open(vis):
        raise ValueError("Failed to open {}".format(vis))
    phasecenter = md.phasecenter()
    mydirection = phasecenter['refer'] +' '+ SkyCoord(phasecenter['m0']['value'], phasecenter['m1']['value'], unit="rad").to_string('hmsdms')
    myfreq = "{:.2f}GHz".format(np.mean(read_spw(vis)))
    print(mydirection)
    print(myfreq)

    os.system('rm -rf random_points.cl')
    cl.done()
    for d in make_random_source(mydirection, n=n, radius=15*u.arcsec):
        print(d, d.encode('utf-8'))
        cl.addcomponent(dir=d.encode('utf-8'), flux=5, fluxunit='mJy', 
                        freq=myfreq, shape='point')
    cl.rename('random_points.cl')
    cl.done()

    ft(vis=vis, complist='random_points.cl')
    uvsub(vis=vis, reverse=True)
    
    vis_new = os.path.join(tmpdir, vis+'.new')
    rmtables(vis_new)
    split(vis=vis, datacolumn='corrected', outputvis=vis_new)
    if make_image:
        make_cont_img(vis_new, outdir=tmpdir, clean=True)

