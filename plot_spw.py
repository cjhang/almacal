# this program try to plot all the spw distributions in almacal
import os
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table

# read the time information of all the calibrators
almacal_info_file = '/Users/jchen/Documents/projects/almacal/data/almacal_timeOnSource.txt'
all_obs = Table.read(almacal_info_file, format='ascii')


# read the spw information in each calibrator
base_dir = '/Users/jchen/Public/almacal_test/almacal_stat'

plot_bands = ['B5', 'B6', 'B7', 'B8']

if True: # plot the band information
    band_in_plot = plot_bands
    fig = plt.figure(figsize=(3*len(plot_bands),5))
    # fig.suptitle(os.path.basename(objfolder))
    ax = fig.add_subplot(111)

    #ALMA band information, in GHz
    band_list = {'B3':[84, 116], 'B4':[125, 163], 'B5':[163, 211], 
            'B6':[211, 275], 'B7':[275, 373], 'B8':[385, 500], \
            'B9':[602, 720], 'B10':[787, 950]}

    band_min = 1000 # arbitraty number beyond alma true band range
    band_max = 10
    for band in band_in_plot:
        if band_min > band_list[band][0]:
            band_min = np.min(band_list[band])
        if band_max < band_list[band][1]:
            band_max = np.max(band_list[band])
    
        ax.broken_barh([(band_list[band][0], np.diff(band_list[band])[0]),], 
                       (0, 1), facecolor='lightblue', edgecolor='grey', \
                       linewidth=1, alpha=0.2)
        ax.text(np.mean(band_list[band]), 1.1, "Band"+band[1:], 
                horizontalalignment='center', verticalalignment='center')
    ax.set_xlim(band_min-10, band_max+10)
    ax.set_ylim(-0.2, 1.2)
    ax.set_xlabel('Frequency [GHz]')
    ax.set_ylabel(r'$t_{\rm on\,source}$ fraction')
    ax.tick_params(axis='y', labelcolor='w', top='off', bottom='on', left='off', right='off', labelsize=2)


for obj in os.listdir(base_dir)[:]:
    # if obj != 'J0522-3627':
        # continue
    try:
      objfits = fits.open(base_dir+'/'+obj+'/{}.fits'.format(obj))
    except:
      print("Error reading fits file of ", obj)

    for band in plot_bands:
        # total_time = 0
        # print("Band: ", band)
        band_header = objfits[band].header
        if band_header['total'] < 1e-6:
          continue
        band_data = objfits[band].data
        # print(band_data)

        band_cal_time = band_header['total']
        band_allcals_time = np.sum(all_obs[band])
        h = 0
        time_slot = None
            
        ax.text(np.mean(band_list[band]), -0.1, 
                r"$t_{{\rm total}}$ = {:.2f} h".format(band_allcals_time/60), \
                horizontalalignment='center', verticalalignment='center')
        
        for data in band_data:
            if time_slot is None:
                time_slot = data[-1]
                # total_time = total_time + time_slot
            dh = time_slot/band_cal_time
            # alpha = np.log(band_cal_time) / np.log(band_allcals_time)
            #alpha = 3 * band_cal_time / band_allcals_time
            alpha = 0.003
            ax.broken_barh([(data[0], np.diff(data[:-1])[0]),], 
                           (h, dh), facecolors='salmon', edgecolors='none', \
                           alpha=alpha)
            if data[-1] != time_slot:
                h = h + dh
                # ax.hlines(y=h, xmin=band_list[band][0], xmax=band_list[band][1], 
                        # color='r', linestyle='-', alpha=0.1, linewidth=1)
                time_slot = data[-1]
                # total_time = total_time + time_slot
        # print('total_time:', total_time)
        # print('band_cal_time', band_cal_time)
    objfits.close()

plt.show()
fig.savefig('spw_all.png', bbox_inches='tight', dpi=400)
# if figname:
    # fig.savefig(figname, bbox_inches='tight', dpi=400)
    # plt.close(fig)

