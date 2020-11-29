# this program try to plot all the spw distributions in almacal
import os
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy import constants as const

# read the time information of all the calibrators
almacal_info_file = '/Users/jchen/Documents/projects/almacal/data/almacal_timeOnSource.txt'
all_obs = Table.read(almacal_info_file, format='ascii')


# read the spw information in each calibrator
base_dir = '/Users/jchen/Public/almacal_test/almacal_stat'

plot_bands = ['B5', 'B6', 'B7', 'B8']

if True: # plot the band information
    band_in_plot = plot_bands
    # fig = plt.figure(figsize=(3*len(plot_bands),5))
    # fig.suptitle(os.path.basename(objfolder))
    # ax = fig.add_subplot(111)

    #ALMA band information, in GHz
    band_list = {'B3':[84, 116], 'B4':[125, 163], 'B5':[163, 211], 
            'B6':[211, 275], 'B7':[275, 373], 'B8':[385, 500], \
            'B9':[602, 720], 'B10':[787, 950]}
B6_freq_range = np.linspace(band_list['B6'][0], band_list['B6'][1], np.diff(band_list['B6'])+1)
B7_freq_range = np.linspace(band_list['B7'][0], band_list['B7'][1], np.diff(band_list['B7'])+1)
cal_counts_B6 = np.zeros_like(B6_freq_range)
cal_counts_B7 = np.zeros_like(B7_freq_range)
cal_counts_B6_time = np.zeros_like(B6_freq_range)
cal_counts_B7_time = np.zeros_like(B7_freq_range)
for obj in os.listdir(base_dir)[:]:
    try:
        with open(base_dir+'/'+obj+'/{}.json'.format(obj)) as f:
            fj = json.load(f)
        spw_discrete = discrete_spw(fj)
        spw_discrete_time = discrete_spw(fj, return_time=True)
    except:
        print("Error in ", obj)
        continue
    cal_counts_B6 = cal_counts_B6 + np.array(spw_discrete['B6'])
    cal_counts_B7 = cal_counts_B7 + np.array(spw_discrete['B7'])
    cal_counts_B6_time = cal_counts_B6_time + np.array(spw_discrete_time['B6'])
    cal_counts_B7_time = cal_counts_B7_time + np.array(spw_discrete_time['B7'])


if True:
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(211)
    ax.set_title('Band 6')
    ax.step(B6_freq_range, cal_counts_B6, where='mid', color='blue', label='Calibrator Numbers')
    ax.set_xlabel('Frequency')
    ax.tick_params(axis='y', labelcolor='blue')
    ax.set_ylabel('Numbers of Calibrator', color='blue')
    B6_fwhm = 1.22*((const.c/(250*u.GHz)).to(u.mm)/(12*u.m)).decompose() * 206265
    ax2 = ax.twinx()
    # ax2.step(B6_freq_range, cal_counts_B6 * 1.5 * B6_fwhm, linewidth=0)
    ax2.step(B6_freq_range, cal_counts_B6_time, where='mid', color='red', label='Total time')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_ylabel('Total Exposure Time [Min]', color='red')
    #plt.legend()
    # ax2.set_ylabel(r'Total Mapping Area [arcsec$^2$]')
    # factor = 1.5 #* B6_fwhm
    # secax = ax.secondary_yaxis('right', functions=(forward, inverse))
    # secax.set_ylabel('radians')
    for center_freq in [228.0311, 230.06710, 244.31904, 246.35503]:
        ax.broken_barh([(center_freq-1.875/2, 1.875),], 
                           (0, 400), facecolors='salmon', edgecolors='none', \
                           alpha=0.2)
    ax.set_ylim(0, 400)

    ax = fig.add_subplot(212)
    ax.set_title('Band 7')
    ax.step(B7_freq_range, cal_counts_B7, where='mid', color='blue', label='Calibrator Numbers')
    ax.set_xlabel('Frequency')
    ax.tick_params(axis='y', labelcolor='blue')
    ax.set_ylabel('Numbers of Calibrator', color='blue')
    ax2 = ax.twinx()
    ax2.step(B7_freq_range, cal_counts_B7_time, where='mid', color='red', label='Total time')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_ylabel('Total Exposure Time [Min]', color='red')
    # B7_fwhm = 1.22*((const.c/(344*u.GHz)).to(u.mm)/(12*u.m)).decompose() * 206265
    # ax2.step(B7_freq_range, cal_counts_B7 * 1.5 * B6_fwhm, linewidth=0)
    # ax2.set_ylabel(r'Total Mapping Area [arcsec$^2$]')
    for center_freq in [342.55565, 344.59165, 354.77161, 356.80760]:
        ax.broken_barh([(center_freq-1.875/2, 1.875),], 
                           (0, 400), facecolors='salmon', edgecolors='none', \
                           alpha=0.2)
    # plt.legend()
    ax.set_ylim(0, 320)
     
    plt.show()
    fig.savefig('spw_statistics.pdf', bbox_inches='tight')
