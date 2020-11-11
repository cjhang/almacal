# This file include the several help functions to read the basic information from the casa asdm file

import os
import re
import json
import numpy as np
import matplotlib.pyplot as plt 
import analysisUtils as au

plt.ioff()

def read_spw(vis):
    """
    """
    tb = tbtool()
    tb.open(vis + '/SPECTRAL_WINDOW')
    col_names = tb.getvarcol('NAME')
    col_freq = tb.getvarcol('CHAN_FREQ')
    tb.close()

    #print(col_names)

    spw_specrange = {}

    for key in col_names.keys():
        freq_max = np.max(col_freq[key]) / 1e9
        freq_min = np.min(col_freq[key]) / 1e9
        spw_specrange[key] = [freq_min, freq_max]

    return spw_specrange.values()

def read_refdir(vis):
    """read the reference direction and return the standard direction string
    """
    tb = tbtool()
    tb.open(vis+'/FIELD')
    reference_dir = tb.getcol('REFERENCE_DIR').flatten()
    tb.close()
    
    rad2deg = 180./np.pi
    direction = "J2000 " + SkyCoord(reference_dir[0]*rad2deg, refval[1]*rad2deg, 
                                      unit="deg").to_string('hmsdms')
    return direction

def spw_stat(objfolder, plot=False, plotbands=['B5', 'B6', 'B7', 'B8'], 
             figname=None, showfig=False, filename=None, savedata=False):
    """make the statistics about one calibrator

    """

    base_dir = objfolder
    obj = os.path.basename(objfolder)

    p_obs = re.compile('uid___')
    spw_list = {'B3':{'name':[], 'time':[], 'freq':[]}, 
                'B4':{'name':[], 'time':[], 'freq':[]}, 
                'B5':{'name':[], 'time':[], 'freq':[]}, 
                'B6':{'name':[], 'time':[], 'freq':[]}, 
                'B7':{'name':[], 'time':[], 'freq':[]}, 
                'B8':{'name':[], 'time':[], 'freq':[]}, 
                'B9':{'name':[], 'time':[], 'freq':[]}, 
                'B10':{'name':[], 'time':[], 'freq':[]},} 

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
                time_on_source = au.timeOnSource(obs_filename, verbose=False, debug=False)
                time_minutes = time_on_source[0]['minutes_on_source']
                if time_minutes < 1e-6:
                    print('No valid on source time!')
                    continue
                spw_list[band]['time'].append(time_minutes)
                spw_list[band]['name'].append(obs)
                spw_specrange = read_spw(obs_filename)
                spw_list[band]['freq'].append(list(spw_specrange))
            except:
                print("Error: in", obs_filename)
    if plot:
        band_in_plot = plotbands
        fig = plt.figure(figsize=(3*len(band_in_plot),5))
        fig.suptitle(os.path.basename(objfolder))
        ax = fig.add_subplot(111)

        #ALMA band information, in GHz
        band_list = {'B3':[84, 116], 'B4':[125, 163], 'B5':[163, 211], 
                'B6':[211, 275], 'B7':[275, 373], 'B8':[385, 500], \
                'B9':[602, 720], 'B10':[787, 950]}

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
        ax.set_ylabel(r'$t_{\rm on\,source}$ fraction')
        ax.tick_params(axis='y', labelcolor='w', top='off', bottom='on', left='off', right='off', labelsize=2)

        for band in band_in_plot:
            h = 0
            band_total_time = np.sum(spw_list[band]['time'])
            ax.text(np.mean(band_list[band]), -0.1, 
                    r"$t_{{\rm total}}$ = {:.2f} min".format(band_total_time), \
                    horizontalalignment='center', verticalalignment='center')
            n_obs = len(spw_list[band]['time'])
            if n_obs < 1:
                continue
            for i in range(n_obs):
                dh = spw_list[band]['time'][i]/band_total_time
                for spw in spw_list[band]['freq'][i]:
                    ax.broken_barh([(spw[0], np.diff(spw)[0]),], 
                                   (h, dh), facecolors='salmon', edgecolors='none', \
                                   alpha=0.5)
                    ax.hlines(y=h, xmin=band_list[band][0], xmax=band_list[band][1], 
                            color='r', linestyle='-', alpha=0.1, linewidth=1/n_obs)
                h = h + dh
        if showfig:
            plt.show()
        if figname:
            fig.savefig(figname, bbox_inches='tight')
            plt.close(fig)

    if savedata:
        with open(filename, 'w') as fp:
            json.dump(spw_list, fp)
    return spw_list

def discrete_spw(spw_list, return_time=False):
    """this function used for statistics of the histogram of spw distribution
    """
    band_list = {'B3':[84, 116], 'B4':[125, 163], 'B5':[163, 211], 
            'B6':[211, 275], 'B7':[275, 373], 'B8':[385, 500], \
            'B9':[602, 720], 'B10':[787, 950]}
    count_list = {}
    count_time_list = {'B3':0, 'B4':0, 'B5':0, 'B6':0, 'B7':0, 'B8':0, 'B9':0, 
                       'B10':0}
    for band in spw_list.keys():
        freq_low, freq_high = band_list[band]
        freq_range = np.linspace(freq_low, freq_high, np.diff(band_list[band])+1, dtype=int)
        freq_counts = np.zeros_like(freq_range, dtype=int)
        freq_cumulate_time = np.zeros_like(freq_range, dtype=float)
        for idx, obs in enumerate(spw_list[band]['freq']):
            for spw in obs:
                covered_freq = (freq_range >= spw[0]) & (freq_range <= spw[1])
                # freq_counts[covered_freq] = freq_counts[covered_freq] + 1
                freq_counts[covered_freq] = 1
                freq_cumulate_time[covered_freq] = freq_cumulate_time[covered_freq] + spw_list[band]['time'][idx]
        count_list[band] = freq_counts.tolist()
        count_time_list[band] = freq_cumulate_time
    if return_time:
        return count_time_list
    return count_list

