# This file include the several help functions to read the basic information from the casa asdm file

import os
import re
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

    #print(col_names)

    spw_specrange = {}

    for key in col_names.keys():
        freq_max = np.max(col_freq[key]) / 1e9
        freq_min = np.min(col_freq[key]) / 1e9
        spw_specrange[key] = [freq_min, freq_max]

    return spw_specrange

def spw_stat(objfolder, plot=False, plotbands=['B5', 'B6', 'B7', 'B8'], 
             figname=None, showfig=False, filename=None, savedata=False):
    # make the statistics about one calibrator

    base_dir = objfolder
    obj = os.path.basename(objfolder)

    p_obs = re.compile('uid___')
    spw_list = {'B3':{'time':[], 'freq':[]}, 
                'B4':{'time':[], 'freq':[]}, 
                'B5':{'time':[], 'freq':[]}, 
                'B6':{'time':[], 'freq':[]}, 
                'B7':{'time':[], 'freq':[]}, 
                'B8':{'time':[], 'freq':[]}, 
                'B9':{'time':[], 'freq':[]}, 
                'B10':{'time':[], 'freq':[]},} 

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
                spw_specrange = read_spw(obs_filename)
                spw_list[band]['freq'].append(list(spw_specrange.values()))
            except:
                print("Error: in", obs_filename)
    if plot:
        band_in_plot = plotbands
        fig = plt.figure(figsize=(12,5))
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

    if savedata: #save the spw_list as fits file
        table_list = []
        from astropy.io import fits
        from astropy.table import Table
        hdr = fits.Header()
        hdr['OBJ'] = obj
        hdr['COMMENT'] = "SPW and exposure time information"
        table_list.append(fits.PrimaryHDU(header=hdr))
        for band in spw_list.keys(): 
            hdr_band = fits.Header()
            hdr_band['BAND'] = band
            hdr_band['TOTAL'] = np.sum(spw_list[band]['time'])
            spw_array = spw_list[band]['freq']
            time_array = spw_list[band]['time']

            for i,obs in enumerate(spw_array):
                obs_time = time_array[i]
                for spw_range in obs:
                    spw_range.append(obs_time)

            if len(spw_array) > 0:
                time_freq_table = np.concatenate(spw_array)
                table_list.append(fits.BinTableHDU(name=band, data=Table(time_freq_table), 
                                  header=hdr_band))
            else:
                table_list.append(fits.BinTableHDU(name=band, 
                                  header=hdr_band))
        

        hdus = fits.HDUList(table_list)
        if filename:
            hdus.writeto(filename, overwrite=True)


    return spw_list


if __name__ == '__main__':
    pass
