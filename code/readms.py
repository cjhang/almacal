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

def spw_stat(obsfolder, plot=False):
    p_obj = re.compile('J\d+[+-]\d')
    p_obs = re.compile('uid___')

    base_dir = obsfolder

    spw_list = {'B3':[], 'B4':[], 'B5':[], 'B6':[],
            'B7':[], 'B8':[], 'B9':[], 'B10':[]}

    for obj in os.listdir(base_dir):
        if p_obj.match(obj):
            for obs in os.listdir(base_dir +'/'+ obj):
                if p_obs.match(obs):
                    obs_filename = base_dir +'/'+ obj+'/'+obs
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
        fig = plt.figure()
        ax = fig.add_subplot(111)

        #ALMA band information, in GHz
        band_list = {'B3':[84, 116], 'B4':[125, 163], 'B5':[163, 211], 'B6':[211, 275],
                'B7':[275, 373], 'B8':[385, 500], 'B9':[602, 720], 'B10':[787, 950]}

        band_min = 1000
        band_max = 10
        for band in ['B6', 'B7']:
            if band_min > band_list[band][0]:
                band_min = np.min(band_list[band])
            if band_max < band_list[band][1]:
                band_max = np.max(band_list[band])
        
            ax.broken_barh([(band_list[band][0], np.diff(band_list[band])[0]),], 
                           (0, 1), facecolors='grey', edgecolors=None)
        ax.set_xlim(band_min-10, band_max+10)
        ax.set_ylim(-0.2, 1.2)

        for band in ['B6', 'B7']:
            n_obs = len(spw_list[band])
            h = 0
            dh = 1./n_obs
            print("number of obs:", n_obs)
            print('dh = ', dh)
            for obs in spw_list[band]:
                print(obs)
                for spw in obs:
                    print(spw)
                    ax.broken_barh([(spw[0], np.diff(spw)[0]),], 
                                   (h, dh), facecolors='red', edgecolors=None)
                h = h + dh

        plt.show()





    return spw_list


if __name__ == '__main__':
    pass
