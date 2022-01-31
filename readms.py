# This file include the several help functions to read the basic information from the casa asdm file
# Copywrite @ Jianhang Chen
# Email: cjhastro@gmail.com

# History:
# 2019.09.27: first release, v0.1


import os
import re
import json
import numpy as np
import matplotlib.pyplot as plt 
import analysisUtils as au

def check_intent(vis, unique=True):
    tb = tbtool()
    tb.open(vis)
    intent_idx = np.unique(tb.getcol('STATE_ID'))
    tb.close()
    tb.open(vis+'/STATE')
    intents = tb.getcol('OBS_MODE')
    intents_list = []
    for idx in intent_idx:
        intents_list.append(intents[idx])
    # print(intents_list)

    intents_valid = []
    for item in intents_list:
        # print('item', item)
        if "BANDPASS" in item:
            intents_valid.append('BANDPASS')
        if "PHASE" in item:
            intents_valid.append('PHASE')
        if "FLUX" in item:
            intents_valid.append('FLUX')
    if unique:
        intents_valid = np.unique(intents_valid).tolist()
    return intents_valid

def read_spw(vis):
    """read the spectral windows
    """
    if isinstance(vis, str):
        vis = [vis, ]
    
    if not isinstance(vis, list):
        raise ValueError("read_spw: Unsupported measurements files!")
    
    spw_specrange = {}

    for v in vis:
        tb.open(v + '/SPECTRAL_WINDOW')
        col_names = tb.getvarcol('NAME')
        col_freq = tb.getvarcol('CHAN_FREQ')
        tb.close()

        for key in col_names.keys():
            freq_max = np.max(col_freq[key]) / 1e9
            freq_min = np.min(col_freq[key]) / 1e9
            spw_specrange[key] = [freq_min, freq_max]

    return spw_specrange.values()

def read_refdir(vis, return_coord=False):
    """read the reference direction and return the standard direction string
    """
    tb = tbtool()
    tb.open(vis+'/FIELD')
    reference_dir = tb.getcol('REFERENCE_DIR').flatten()
    tb.close()
    
    rad2deg = 180./np.pi
    ref_coord = SkyCoord(reference_dir[0]*rad2deg, reference_dir[1]*rad2deg, 
                                      unit="deg")
    direction = "J2000 " + ref_coord.to_string('hmsdms')
    if return_coord:
        return ref_coord

    return direction

def spw_stat(vis=None, jsonfile=None, plot=False, savefile=None, 
        bands=['B3','B4','B5', 'B6', 'B7', 'B8','B9','B10'], figname=None, showfig=False,  
        time_select=False, start_time='2010-01-01T00:00:00', end_time='2050-01-01T00:00:00', 
        z=0, lines=None, lines_names=None, exclude_aca=True, debug=False):
    """make the statistics about one calibrator

    Args:
        objfolder (str): the folder contains all the visibility
        vis: a single visibility or a list
        jsonfile: the json file to be read
        plot (bool): whehther plot the results interactively (default is False)
            - plotbands (list): the band to be plotted (default is ['B5', 'B6', 
              'B7', 'B8'])
            - figname (bool): the figname of the saved figure, if None, no figures 
              will be saved (default is None)
        savedata (bool): whether save the statistics into file, default to save as
            json file
            - filename (str): the filename of the saved data
        
        lines (list): the spectral lines to be added to the plots, in units of GHz,
                      like: [115.3, 230.5, 345.8]
        lines_names (list): the names of the lines, like: ['CO1-0, CO2-1, CO3-2']
        z (float): the redshift of the source
    

    """
    spw_list = {}
    for band in bands:
        spw_list[band] = {'name':[], 'time':[], 'freq':[]}
    filelist = []

    if vis is not None:
        if debug:
            print('vis', vis)
        if isinstance(vis, str):
            filelist = [vis,]
        elif isinstance(vis, list):
            filelist = vis
    elif jsonfile:
        if debug:
            print('jsonfile', jsonfile)
        with open(jsonfile, 'r') as f:
            spw_list = json.load(f)
    else:
        raise ValueError("No valid files have been given!")

    band_match = re.compile('_(?P<band>B\d{1,2})')
    for obs in filelist:
        try:
            if band_match.search(obs):
                band = band_match.search(obs).groupdict()['band']
                if debug:
                    print("Band: ", band)
            else:
                if debug:
                    print("Error in band match.")
                continue
            if band not in bands:
                continue

            if exclude_aca:
                try:
                    tb.open(obs + '/ANTENNA')
                    antenna_diameter = np.mean(tb.getcol('DISH_DIAMETER'))
                    tb.close()
                except:
                    continue
                if antenna_diameter < 12.0:
                    if debug:
                        print("Excuding data from {}".format(antenna_diameter))
                    continue

            if time_select:
                start_time = Time(start_time)
                end_time = Time(end_time)
                try:
                    tb.open(obs)
                    obs_time = Time(tb.getcol('TIME').max()/24/3600, format='mjd')
                    tb.close()
                except:
                    if debug:
                        print("Error in opening the visibility!")
                    continue
                if debug:
                    print('> obs_time', obs_time.iso)
                    print('> start_time', start_time.iso)
                    print('> end_time', end_time.iso)
                if obs_time < start_time or obs_time > end_time:
                    if debug:
                        print(">> Skip by wrong observation time: {}".format(obs))
                    continue

            time_on_source = au.timeOnSource(obs, verbose=False, debug=False)
            time_minutes = time_on_source[0]['minutes_on_source']
            if debug:
                print('time_on_source', time_on_source)
            if time_minutes < 1e-6:
                print('No valid on source time!')
                continue
            spw_list[band]['time'].append(time_minutes)
            spw_list[band]['name'].append(os.path.basename(obs))
            spw_specrange = read_spw(obs)
            spw_list[band]['freq'].append(list(spw_specrange))
        except:
            print("Error: in", obs)
    if plot:
        fig = plt.figure(figsize=(3*len(bands),5))
        # fig.suptitle(os.path.basename(objfolder))
        ax = fig.add_subplot(111)

        #ALMA band information, in GHz
        band_list = {'B3':[84, 116], 'B4':[125, 163], 'B5':[163, 211], 
                'B6':[211, 275], 'B7':[275, 373], 'B8':[385, 500], \
                'B9':[602, 720], 'B10':[787, 950]}

        band_min = 1000
        band_max = 10
        for band in bands:
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
        ax.tick_params(axis='x', which='major', labelsize=8)
        ax.tick_params(axis='x', which='minor', labelsize=6)
        ax.tick_params(axis='y', labelcolor='w', top='off', bottom='on', left='off', right='off', labelsize=2)

        # plot the spectral window
        for band in bands:
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
        
        # plot the spectral lines
        if lines:
            if isinstance(lines, (float,int)):
                lines = [lines,]
            if lines_names:
                if isinstance(lines_names, (str)):
                    lines = [lines_names,]
            for idx, line in enumerate(lines):
                line_obs = line / (1.0 + z)
                ax.vlines(line_obs, 0, 1, alpha=0.6, color='b', linewidth=1)
                if lines_names:
                    ax.text(line_obs, 0.8, lines_names[idx], fontsize=12, alpha=0.6, horizontalalignment='center')

        if showfig:
            plt.show()
        if figname:
            fig.savefig(figname, bbox_inches='tight')
            plt.close(fig)

    if savefile is not None:
        with open(savefile, 'w') as fp:
            json.dump(spw_list, fp)
    return spw_list

def discrete_spw(spw_list, bands=None, freq_range=None, width=1, return_time=False):
    """this function used for statistics of the histogram of spw distribution

    The default width is in unit of GHz
    """
    band_list = {'B3':[84, 116], 'B4':[125, 163], 'B5':[163, 211], 
            'B6':[211, 275], 'B7':[275, 373], 'B8':[385, 500], \
            'B9':[602, 720], 'B10':[787, 950]}
    if bands is None:
        bands = band_list.keys()
    stats_bands = {}

    for band in bands:
        if freq_range is None:
            freq_low, freq_high = band_list[band]
            freq_array = np.arange(freq_low, freq_high+width, width, dtype=int)
        else:
            freq_array = np.arange(freq_range[0], freq_range[-1]+width, width, dtype=int)
        freq_counts = np.zeros_like(freq_array, dtype=int)
        freq_cumulate_time = np.zeros_like(freq_array, dtype=float)
        for idx, obs in enumerate(spw_list[band]['freq']):
            for spw in obs:
                covered_freq = (freq_array >= spw[0]) & (freq_array <= spw[1])
                # freq_counts[covered_freq] = freq_counts[covered_freq] + 1
                freq_counts[covered_freq] += 1
                freq_cumulate_time[covered_freq] = (freq_cumulate_time[covered_freq] 
                                                    + spw_list[band]['time'][idx])
        stats_bands[band] = np.array([freq_array, freq_counts, freq_cumulate_time])
    return stats_bands

def readms(vis, bmaj=None, bmin=None):
    """The uniform function to get all the neccessary information for a given visibility
    """
    ms_info = {}
    ms_info['flux_units'] = 'mJy'
    ms_info['sensitivity'] = 1000 * calculate_sensitivity(vis) # convert into mJy
    beamsize = np.pi/(4*np.log(2))*bmaj*bmin 
    ms_info['peak'] = ms_info['sensitivity'] /1000 * beamsize / (2*np.pi*bmaj*bmin)
    return ms_info

def readimage(image):
    """The uniform function to read images
    """
    im_info = {}
    im_head = imhead(image)
    im_beam = im_head['restoringbeam']
    im_incr = im_head['incr']
    im_stat = imstat(image)
    # beamsize = np.pi*a*b/(4*np.log(2))
    beamsize = np.pi/(4*np.log(2))* im_beam['major']['value'] * im_beam['minor']['value'] / (im_incr[0]/np.pi*180*3600)**2
    im_info['beamsize'] = beamsize
    rms = im_stat['rms'] * 1000 # in mJy/beam
    im_info['flux_units'] = 'mJy'
    im_info['sensitivity'] = rms * beamsize
    im_info['rms'] = rms
    im_info['bmaj'] = im_beam['major']['value'] / (im_incr[0]/np.pi*180*3600) 
    im_info['bmin'] = im_beam['minor']['value'] / (im_incr[0]/np.pi*180*3600) 
    im_info['peak'] = im_stat['max']
    return im_info

def read_num(vis):
    all_num = []
    for item in vis:
        if isinstance(item, list):
            all_num.append(len(item))
        else:
            all_num.append(1)
    return all_num

def read_onSourceTime(vis):
    all_onSourceTime = []
    for item in vis:
        if isinstance(item, list):
            item_results = []
            for obs in item:
                try:
                    time_on_source = au.timeOnSource(obs, verbose=False, debug=False)
                    time_minutes = time_on_source[0]['minutes_on_source']
                except:
                    time_minutes = 0.0
                item_results.append(time_minutes)
            all_onSourceTime.append(item_results)
        else:
            try:
                time_on_source = au.timeOnSource(item, verbose=False, debug=False)
                time_minutes = time_on_source[0]['minutes_on_source']
            except:
                time_minutes = 0.0
            all_onSourceTime.append(time_minutes)

    return all_onSourceTime

def read_flux(vis):
    flux_list = []
    for item in vis:
        if isinstance(item, list):
            item_list = []
            for obs in item:
                obs_flux = search_flux(obs)
                item_list.append(obs_flux)
            flux_list.append(item_list)
        else:
            flux_list.append(search_flux(item))
    return flux_list

def read_radec(vis):
    if isinstance(vis, list):
        radec_list = []
        for item in vis:
            if isinstance(item, list):
                item_list = []
                for obs in item:
                    obs_refdir = read_refdir(obs, return_coord=True)
                    item_list.append((obs_refdir.ra.value, obs_refdir.dec.value))
                radec_list.append(item_list)
            else:
                item_refdir = read_refdir(item, return_coord=True)
                radec_list.append(item_refdir.ra.value, item_refdir.dec.value)
    else:
        vis_refdir = read_refdir(vis, return_coord=True)
        radec_list = [vis_refdir.ra.value, vis_refdir.dec.value]
    return radec_list

def read_azel(vis):
    if isinstance(vis, list):
        alel_list = []
        for item in vis:
            if isinstance(item, list):
                item_list = []
                for obs in item:
                    obs_alel = au.computeAzElForMS(obs)
                    item_list.append(obs_alel)
                alel_list.append(item_list)
            else:
                alel_list.append(au.computeAzElForMS(item))
        return alel_list
    else:
        return au.computeAzElForMS(vis)

def read_intents(vis):
    if isinstance(vis, list):
        intents_list = []
        for item in vis:
            if isinstance(item, list):
                item_list = []
                for obs in item:
                    obs_intent = check_intent(obs)
                    item_list.append(obs_intent)
                intents_list.append(item_list)
            else:
                intents_list.append(check_intent(item))
        return intents_list
    else:
        return check_intent(vis)

