import os
import matplotlib.pyplot as plt
import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u

plt.style.use('dark_background')

almacal_info_file = os.path.join(os.path.expanduser('~'), 
        'Documents/projects/almacal/data/almacal_timeOnSource.txt')

data = np.loadtxt(almacal_info_file, skiprows=1,
        dtype={'names': ('obj', 'B3','B4','B5','B6','B7','B8','B9','B10'), 
               'formats': ('S10', 'f4', 'f4','f4','f4','f4','f4','f4','f4')})

bands = ['B6',]
for band in bands:
# if True:
    
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
        depth.append(row['{}'.format(band)])

    #print('depth', depth)
    # print(x_coords)
    obs_valid = np.array(depth)>1e-6
    xx = np.array(x_coords)[obs_valid]
    yy = np.array(y_coords)[obs_valid]
    depth_valid = np.log10(np.array(depth)[obs_valid])



    #make figure
    fig = plt.figure(figsize=(12,4))
    fig.suptitle('Band {}'.format(band), fontsize=16)
    ax = fig.add_subplot(121, projection='hammer' ) # or mollweide, aitoff, or lambert. [example](http://matplotlib.org/examples/pylab_examples/geo_demo.html) 
    im = ax.scatter(xx, yy, marker='o', s=5, c=depth_valid, cmap='cool', alpha=0.9)
    cbar = plt.colorbar(im)
    cbar.set_label(r'$\log(t_{\rm obs} {\rm [Minutes]})$')
    ax.grid(True)
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']) #use if you want to change to a better version of RA. 
    ax.set_xlabel('R.A.')
    ax.set_ylabel(r'$\delta$')

    ax = fig.add_subplot(122)
    im = plt.hist(depth_valid, bins=20)
    plt.xticks([0, 1, 2, 3], [r'$1$', r'$10$', r'$10^2$', r'$10^3$'])
    ax.set_xlabel(r'$\log (t_{\rm obs}{\rm [Minutes]})$')
    ax.set_ylabel('Counts')

    plt.show()
    # fig.savefig('/tmp/Band{}_stat.pdf'.format(band))
