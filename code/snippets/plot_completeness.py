
import json
import numpy as np
import matplotlib.pyplot as plt

J1303_5540_B7_daostarfinder = '../../data/simulation/J1303-5540_B7_threshold5_daostarfinder.json'
J1303_5540_B7_peakfinder = '../../../data/simulation/J1303-5540_B7_threshold5_peakfind.json'
J1933_6942_B7_daostarfinder = '../../../data/simulation/J1933-6942_B7_threshold5_daostarfinder.json'
J1933_6942_B7_peakfinder = '../../../data/simulation/J1933-6942_B7_threshold5_findpeak.json'
J1933_6942_B6_daostarfinder = '../../../data/simulation/J1933-6942_B7_threshold5_daostarfinder.json'
J1933_6942_B6_peakfinder = '../../../data/simulation/J1933-6942_B7_threshold5_findpeak.json'

def plot_completeness(jsonfile, snr=np.arange(1.0, 10, 0.1), ax=None):
    """
    """
    
    with open(jsonfile) as jf:
        data = json.load(jf)

    snr_input = np.array(data['snr_input'])
    flux_input = np.array(data['flux_input'])
    flux_input_aperture = np.array(data['flux_input_aperture'])
    flux_input_gaussian = np.array(data['flux_input_gaussian'])
    flux_aperture = np.array(data['flux_aperture'])
    flux_gaussian = np.array(data['flux_gaussian'])
    detection_snr = np.array(data['detection_snr'])
    detection_input_array = np.array(data['detection_input_array'])
    detection_found_array = np.array(data['detection_found_array'])

    print(snr_input)
    # print('detection_array', detection_array)
    completeness_list = []
    fake_rate_list = []
    for i in range(1, len(snr)):
        s = snr[i]
        s_b = snr[i-1]
        # calculate the completeness
        snr_select = np.bitwise_and((detection_input_array[:, 0]<s), (detection_input_array[:, 0]>s_b))
        n_input = np.sum(snr_select)
        n_found = np.sum(np.array(detection_input_array[:,1][snr_select]))
        completeness_list.append(1.0*n_found/n_input)
        
        # calculate the fake detaction rate
        snr_select2 = np.bitwise_and((detection_found_array[:, 0]<s), (detection_found_array[:, 0]>s_b))
        n_found2 = np.sum(snr_select2)
        n_fake = n_found2 - np.sum(np.array(detection_found_array[:,1][snr_select2]))
        fake_rate_list.append(1.0*n_fake/n_found2)

    print(completeness_list)
    if ax is not None:
        fig, ax = plt.subplots(1,3, figsize=(12, 3))
    
    ax0 = ax[0]
    ax0.set_xlabel('SNR')
    ax0.set_ylabel(r'$S_{\rm out}/S_{\rm in}$')
    ax0.plot(snr_input, flux_input_aperture/flux_input, 'k.', label='aperture')
    ax0.plot(snr_input, flux_input_gaussian/flux_input, 'r.', label='gaussian')
    
    ax1 = ax[1]
    ax1.plot(0.5*(snr[1:]+snr[:-1]), completeness_list, 'o')
    ax1.set_xlabel('SNR')
    ax1.set_ylabel(r'Completeness')
    ax1.set_ylim((-0.1, 1.2))

    ax2 = ax[2]
    ax2.plot(0.5*(snr[1:]+snr[:-1]), fake_rate_list, 'o')
    ax2.set_xlabel('SNR')
    ax2.set_ylabel(r'Fake percentage')
    ax2.set_xlim((0., 8))
    ax2.set_ylim((-0.1, 1.2))

    # plt.show()

if __name__ == '__main__':
    fig, ax = plt.subplots(1,3, figsize=(12,3))
    plot_completeness(J1303_5540_B7_daostarfinder,ax=ax)
    plt.show()
