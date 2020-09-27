# This file include the several help functions to read the basic information from the casa asdm file

import os
import re
import numpy as np
import matplotlib.pyplot as plt 

def read_spw(vis):
    """
    """
    tb = tbtool()
    tb.open(vis + '/SPECTRAL_WINDOW')
    col_names = tb.getvarcol('NAME')
    col_freq = tb.getvarcol('CHAN_FREQ')

    #print(col_names)

    spw_specrange = {}

    for key in col_names.keys():
        freq_max = np.max(col_freq[key]) / 1e9
        freq_min = np.min(col_freq[key]) / 1e9
        spw_specrange[key] = [freq_min, freq_max]

    return spw_specrange

def spw_stat(objfolder, plot=False):
    p_obs = re.compile('uid___')

    base_dir = objfolder

    spw_list = {'B3':[], 'B4':[], 'B5':[], 'B6':[],
            'B7':[], 'B8':[], 'B9':[], 'B10':[]}

    for obs in os.listdir(base_dir +'/'):
        if p_obs.match(obs):
            obs_filename = base_dir +'/'+ '/'+obs
            try:
                band_match = re.compile('_(?P<band>B\d{1,2})')
                if band_match.search(obs):
                    band = band_match.search(obs).groupdict()['band']
                    # print("Band: ", band)
                else:
                    print("Error in band match.")
                spw_specrange = read_spw(obs_filename)
                spw_list[band].append(spw_specrange.values())
            except:
                print("Error: in", obs_filename)
    if plot:
        band_in_plot = ['B5', 'B6', 'B7', 'B8']
        fig = plt.figure(figsize=(12,5))
        ax = fig.add_subplot(111)

        #ALMA band information, in GHz
        band_list = {'B3':[84, 116], 'B4':[125, 163], 'B5':[163, 211], 'B6':[211, 275],
                'B7':[275, 373], 'B8':[385, 500], 'B9':[602, 720], 'B10':[787, 950]}

        band_min = 1000
        band_max = 10
        for band in band_in_plot:
            if band_min > band_list[band][0]:
                band_min = np.min(band_list[band])
            if band_max < band_list[band][1]:
                band_max = np.max(band_list[band])
        
            ax.broken_barh([(band_list[band][0], np.diff(band_list[band])[0]),], 
                           (0, 1), facecolor='lightblue', edgecolor='grey', \
                           linewidth=1, alpha=0.5)
            ax.text(np.mean(band_list[band]), 1.1, "Band"+band[1:], 
                    horizontalalignment='center', verticalalignment='center')
        ax.set_xlim(band_min-10, band_max+10)
        ax.set_ylim(-0.2, 1.2)
        ax.set_xlabel('Frequency [GHz]')
        ax.tick_params(axis='y', labelcolor='w', top='off', bottom='on', left='off', right='off', labelsize=2)

        for band in band_in_plot:
            n_obs = len(spw_list[band])
            if n_obs < 1:
                continue
            h = 0
            dh = 1./n_obs
            for obs in spw_list[band]:
                for spw in obs:
                    ax.broken_barh([(spw[0], np.diff(spw)[0]),], 
                                   (h, dh), facecolors='salmon', edgecolors='none', \
                                   alpha=0.5)
                    ax.hlines(y=h, xmin=band_list[band][0], xmax=band_list[band][1], 
                            color='r', linestyle='-', alpha=0.1, linewidth=0.2)
                h = h + dh

        plt.show()





    return spw_list


if __name__ == '__main__':
    pass
