
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

J1303_5540_B7_daostarfinder = '../../data/simulation/J1303-5540_B7_threshold5_daostarfinder.json'
J1303_5540_B7_peakfinder = '../../data/simulation/J1303-5540_B7_threshold5_peakfinder.json'
J1933_6942_B7_daostarfinder = '../../data/simulation/J1933-6942_B7_threshold5_daostarfinder.json'
J1933_6942_B7_peakfinder = '../../data/simulation/J1933-6942_B7_threshold5_peakfinder.json'
J1933_6942_B6_daostarfinder = '../../data/simulation/J1933-6942_B6_threshold5_daostarfinder.json'
J1933_6942_B6_peakfinder = '../../data/simulation/J1933-6942_B6_threshold5_peakfinder.json'

def plot_fluxboosting(jsonfile, ax=None, snr=np.arange(4, 36, 2.), label=None):
    """
    """
    with open(jsonfile) as jf:
        data = json.load(jf)

    snr_input = np.array(data['snr_input'])
    flux_input = np.array(data['flux_input'])
    flux_input_aperture = np.array(data['flux_input_aperture'])
    flux_input_gaussian = np.array(data['flux_input_gaussian'])

    snr_mid = 0.5*(snr[1:] + snr[:-1])
    aperture_mean = []
    gaussian_mean = []

    for i in range(1, len(snr)):
        s = snr[i]
        s_b = snr[i-1]
        snr_select = np.bitwise_and((snr_input<s), (snr_input>s_b))
        aperture_boosting = flux_input_aperture[snr_select] / flux_input[snr_select]
        gaussian_boosting = flux_input_gaussian[snr_select] / flux_input[snr_select]
        aperture_mean.append([np.median(aperture_boosting), np.std(aperture_boosting)])
        gaussian_mean.append([np.median(gaussian_boosting), np.std(gaussian_boosting)])

    aperture_mean = np.array(aperture_mean)
    gaussian_mean = np.array(gaussian_mean)
    print(snr_mid.shape)
    print(aperture_mean.shape)
    print(gaussian_mean.shape)

    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(4, 4))
    
    ax.set_xlabel('SNR')
    ax.set_ylabel(r'$S_{\rm out}/S_{\rm in}$')
    ax.plot(snr_input, flux_input_aperture/flux_input, 'k.', markersize='1', alpha=0.5)
    ax.plot(snr_input, flux_input_gaussian/flux_input, 'r.', markersize='1', alpha=0.5)
    
    ax.plot(snr_mid, aperture_mean[:,0], 'ko', label=label+' aperture')
    ax.errorbar(snr_mid, aperture_mean[:,0], yerr=aperture_mean[:,1], color='k', lw=2, capsize=5, elinewidth=2, markeredgewidth=2, alpha=0.8)
    ax.plot(snr_mid, gaussian_mean[:,0], 'ro', label=label+' gaussian')
    ax.errorbar(snr_mid, gaussian_mean[:,0], yerr=gaussian_mean[:,1], color='r', lw=2, capsize=5, elinewidth=2, markeredgewidth=2, alpha=0.8)

def plot_completeness(jsonfile, snr=np.arange(1.0, 10, 0.3), ax=None, label=None, color=None):
    """completeness
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

    completeness_list = np.ma.masked_invalid(completeness_list)
    snr_mid = 0.5*(snr[1:]+snr[:-1])
    cs = CubicSpline(snr_mid, completeness_list.filled(0.))


    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(4, 4))
   
    ax.plot(snr_mid, completeness_list, 'o', label=label, color=color)
    ax.plot(snr_mid, cs(snr_mid), '--', color=color)
    ax.set_xlabel('SNR')
    ax.set_ylabel(r'Completeness')
    ax.set_ylim((-0.1, 1.2))
    ax.legend()

def plot_fakerate(jsonfile, snr=np.arange(1.0, 10, 0.3), ax=None, label=None):
    """plot the rate of fake sources
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

    completeness_list = np.ma.masked_invalid(completeness_list)
    snr_mid = 0.5*(snr[1:]+snr[:-1])
    cs = CubicSpline(snr_mid, completeness_list.filled(0.))

    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(4, 4))
    ax.plot(0.5*(snr[1:]+snr[:-1]), fake_rate_list, 'o', label=label)
    ax.set_xlabel('SNR')
    ax.set_ylabel(r'Fake percentage')
    ax.set_xlim((0., 8))
    ax.set_ylim((-0.1, 1.2))
    ax.legend()



if __name__ == '__main__':
    # plot the flux boosting
    fig, ax = plt.subplots(1,2, figsize=(12,4))
    plot_fluxboosting(J1303_5540_B7_daostarfinder,ax=ax[0], label='J1303-5540 daostarfinder B7')
    ax[0].legend()
    plot_fluxboosting(J1933_6942_B7_daostarfinder,ax=ax[1], label='J1933-6942 daostarfinder B7')
    ax[1].legend()
    # plt.show()
    # plt.show()

    # plot the completeness
    fig, ax = plt.subplots(1,2, figsize=(12,4))
    ax[0].set_title('Sextractor')
    plot_completeness(J1303_5540_B7_daostarfinder,ax=ax[0], label='J1303-5540 B7')
    plot_completeness(J1933_6942_B7_daostarfinder,ax=ax[0], label='J1933-6942 B7')
    plot_completeness(J1933_6942_B6_daostarfinder,ax=ax[0], label='J1933-6942 B6')
    
    ax[1].set_title('Peakfinder')
    plot_completeness(J1303_5540_B7_peakfinder,ax=ax[1], snr=np.arange(1.0, 10, 0.1), label='J1303-5540 B7')
    plot_completeness(J1933_6942_B7_peakfinder,ax=ax[1], snr=np.arange(1.0, 10, 0.1), label='J1933-6942 B7')
    plot_completeness(J1933_6942_B6_peakfinder,ax=ax[1], snr=np.arange(1.0, 10, 0.1), label='J1933-6942 B6')


    # plot the fake rate
    fig, ax = plt.subplots(1,1, figsize=(6,4))
    plot_fakerate(J1303_5540_B7_daostarfinder,ax=ax, label='J1303-5540 B7')
    # plot_fakerate(J1933_6942_B7_daostarfinder,ax=ax, label='J1933-6942 B7')
    plot_fakerate(J1933_6942_B6_daostarfinder,ax=ax, label='J1933-6942 B7')
    
    
    plt.show()
