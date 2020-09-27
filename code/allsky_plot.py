import matplotlib.pyplot as plt
import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u

almacal_info_file = '/tmp/almacal.info.txt'

data = np.loadtxt(almacal_info_file, skiprows=1, delimiter=' ',
        dtype={'names': ('obj', 'B3','B4','B5','B6','B7','B8','B9'), 
               'formats': ('S10', 'f4', 'f4','f4','f4','f4','f4','f4')})


band = 5
# for band in [3,4,5,6,7,8,9,10]:
if True:
    
    x_coords = []
    y_coords = []
    depth = []
    
    for row in data:
        obj = (row['obj']).decode('ascii')
        # print('obj', obj)
        # print('{}h{}m00s'.format(obj[1:3], obj[3:5]),
                 # '{}d{}m00s'.format(obj[5:8], obj[8:10]))
        coord = SkyCoord('{}h{}m00s'.format(obj[1:3], obj[3:5]), 
                '{}d{}m00s'.format(obj[5:8], obj[8:10]), frame='icrs')
        ra = coord.ra.wrap_at(180*u.degree)
        dec = coord.dec
        x_coords.append(ra.radian)
        y_coords.append(dec.radian)
        depth.append(row['B{}'.format(band)]/60)

    #print('depth', depth)
    # print(x_coords)
    obs_valid = np.array(depth)>0.1
    xx = np.array(x_coords)[obs_valid]
    yy = np.array(y_coords)[obs_valid]
    depth_valid = np.array(depth)[obs_valid]



    #make figure
    fig = plt.figure(figsize=(12,4))
    fig.suptitle('Band {}'.format(band), fontsize=16)
    ax = fig.add_subplot(121, projection='hammer' ) # or mollweide, aitoff, or lambert. [example](http://matplotlib.org/examples/pylab_examples/geo_demo.html) 
    im = ax.scatter(xx, yy, marker='o', s=5, c=depth_valid, cmap='cool', alpha=0.9)
    cbar = plt.colorbar(im)
    cbar.set_label('Minutes')
    ax.grid(True)
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']) #use if you want to change to a better version of RA. 
    ax.set_xlabel('R.A.')
    ax.set_ylabel(r'$\delta$')

    ax = fig.add_subplot(122)
    im = plt.hist(np.log10(depth_valid), bins=20)
    ax.set_xlabel(r'$\log (t_{\rm obs}[Minutes])$')
    ax.set_ylabel('Counts')

    plt.show()
    fig.savefig('almacal_stat/Band{}_stat.pdf'.format(band))
