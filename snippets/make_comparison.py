# This is the file trying to plot the difference between the data up to now with that of 
# the Iran Oteo in 2016
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

# timeOnSource_all = '../data/make_all_goodimages.txt'
# data_all = np.loadtxt(timeOnSource_all, skiprows=0,
        # dtype={'names': ('obj', 'B3','B4','B5','B6','B7','B8','B9','B10'), 
               # 'formats': ('S10', 'f4', 'f4','f4','f4','f4','f4','f4','f4')})
timeOnSource_all = '../data/all_obstime.txt'
data_all = Table.read(timeOnSource_all, format='ascii')
timeOnSource_main = '../data/main_obstime.txt'
data_main = Table.read(timeOnSource_main, format='ascii')
timeOnSource_oteo2016 = '../data/almacal.oteo2016.info.txt'
data_oteo2016 = np.loadtxt(timeOnSource_oteo2016, skiprows=1, delimiter=' ',
        dtype={'names': ('obj', 'B3','B4','B5','B6','B7','B8','B9'), 
               'formats': ('S10', 'f4', 'f4','f4','f4','f4','f4','f4')})
# final selected data
data_select = {}
timeOnSource_select_B7 = '../data/B7_obstime.txt'
data_select_B7 = Table.read(timeOnSource_select_B7, format='ascii')
data_select['B7'] = data_select_B7
timeOnSource_select_B6 = '../data/B6_obstime.txt'
data_select_B6 = Table.read(timeOnSource_select_B6, format='ascii')
data_select['B6'] = data_select_B6


# printing the fraction of data lose
print("For Band 6:")
print("{:.4f}".format(np.sum(data_select_B6['B6'])/np.sum(data_main['B6'])))
print("For Band 7:")
print("{:.4f}".format(np.sum(data_select_B7['B7'])/np.sum(data_main['B7'])))


# data = data_oteo2016
# save_fig_name = 'data_oteo2016'
# data = data_new
# save_fig_name = 'data_new'

# bands = ['B6', 'B7']

# for the data used in oteo2016
if False:
    band_depth_new = {'B6':[], 'B7':[]}
    band_depth_good = {'B6':[], 'B7':[]}
    band_depth_oteo2016 = {'B6':[], 'B7':[]}
    
    # for new data
    for row in data_new:
        for band in bands:
            obj = (row['obj']).decode('ascii')
            band_depth_new[band].append(row[band]/60) # in minutes
    for band, depth in band_depth_new.items():
        obs_valid = np.array(depth)>0.1
        band_depth_new[band] = np.array(depth)[obs_valid]

    # for good data
    for row in data_good:
        for band in bands:
            obj = (row['obj']).decode('ascii')
            band_depth_good[band].append(row[band]) # in minutes
    for band, depth in band_depth_good.items():
        obs_valid = np.array(depth)>0.1
        band_depth_good[band] = np.array(depth)[obs_valid]

    # for oteo2016
    for row in data_oteo2016:
        for band in bands:
            obj = (row['obj']).decode('ascii')
            band_depth_oteo2016[band].append(row[band]/60) # in minutes
    for band, depth in band_depth_oteo2016.items():
        obs_valid = np.array(depth)>0.1
        band_depth_oteo2016[band] = np.array(depth)[obs_valid]

for band in ['B6','B7']:
    fig = plt.figure(figsize=(6,5))
    plt.style.use('default')
    # fig.suptitle('Comparison', fontsize=16)

    # # plot the oteo2016
    # ax = fig.add_subplot(221)
    # ax.set_title('Oteo et al. (2016)')
    # ax.text(0.1, 9, 'Band 6')
    # im = plt.hist(np.log10(band_depth_oteo2016['B6']), bins=20)
    # ax.set_xlabel(r'$\log (t_{\rm obs}[Minutes])$')
    # ax.set_ylabel('Number')
    # # for 50 minutes, 30 antenna, at 243GHz
    # ax.axvline(x=1.7, color='b', linestyle='-.', alpha=0.5)
    # ax.text(1.7, 7.5, r'$25\,\mu$Jy', color='b', alpha=0.5)
    # ax.set_xlim(0, 1.8)
    # ax.set_ylim(0, 9.5)

    # ax = fig.add_subplot(223)
    # ax.text(0.1, 5.5, 'Band 7')
    # im = plt.hist(np.log10(band_depth_oteo2016['B7']), bins=20)
    # ax.set_xlabel(r'$\log (t_{\rm obs}[Minutes])$')
    # ax.set_ylabel('Number')
    # # for 50 minutes, 30 antenna, at 345GHz
    # ax.axvline(x=1.7, color='r', linestyle='-.', alpha=0.5)
    # ax.text(1.7, 5, r'$39\,\mu$Jy', color='r', alpha=0.5)
    # ax.set_xlim(0, 1.8)
    # ax.set_ylim(0, 6)



    # plot the new data
    ax = fig.add_subplot(111)
    # ax.set_xscale('log')
    ax.set_title(band)
    # ax.text(3.2, 87, 'Band 6')
    data_all_valid = data_all[band][data_all[band] > 0.1]
    data_main_valid = data_main[band][data_main[band] > 0.1]
    data_oteo2016_valid = data_oteo2016[band][data_oteo2016[band] > 0.1] / 60.
    data_select_valid = data_select[band][band][data_select[band][band] > 0.1]
    # im = plt.hist(np.log10(data_all_valid), bins=20, label='All (main+ACA)', alpha=0.5)
    im = plt.hist(np.log10(data_main_valid), bins=20, label='Main', alpha=0.7)
    im = plt.hist(np.log10(data_oteo2016_valid), bins=20, label='Oteo et al. (2016)', alpha=0.9)
    im = plt.hist(np.log10(data_select_valid), bins=20, label='Selected', alpha=0.5)
    ax.set_xlabel(r'$\log (t_{\rm obs}[{\rm Minutes}])$')
    ax.set_ylabel('Number')
    # for 1000 minutes, 43 antenna, at 243GHz
    ax.axvline(x=3, color='b', linestyle='-.', alpha=0.5)
    ax.text(3.1, 60, r'$4\,\mu$Jy', color='b', alpha=0.5)
    # for 50 minutes, 30 antenna, at 243GHz
    ax.axvline(x=1.7, color='b', linestyle='-.', alpha=0.5)
    ax.text(1.7, 60, r'$25\,\mu$Jy', color='b', alpha=0.5)
    ax.set_xlim(0, 3.7)
    ax.set_ylim(0, 95)
    # locs, labels = plt.xticks()
    plt.xticks([0, 1, 2, 3], [r'$1$', r'$10$', r'$10^2$', r'$10^3$'])
    # print(locs, labels)
    ax.legend()
    plt.show()
    #fig.savefig('../results/comparison.pdf')

    # ax = fig.add_subplot(122)
    # # ax.text(3.2, 55, 'Band 7')
    # ax.set_title('Band 7')
    # im = plt.hist(np.log10(band_depth_new['B7']), bins=20, alpha=0.5)
    # im = plt.hist(np.log10(band_depth_good['B7']), bins=20, alpha=0.7)
    # im = plt.hist(np.log10(band_depth_oteo2016['B7']), bins=20, alpha=0.9)
    # plt.xticks([0, 1, 2, 3], [r'$1$', r'$10$', r'$10^2$', r'$10^3$'])
    # ax.set_xlabel(r'$\log (t_{\rm obs}[{\rm Minutes}])$')
    # ax.set_ylabel('Number')
    # # for 1000 minutes, 43 antenna, at 345 GHz
    # ax.axvline(x=3, color='r', linestyle='-.', alpha=0.5)
    # ax.text(3.1, 50, r'$6\,\mu$Jy', color='r', alpha=0.5)
    # # for 50 minutes, 30 antenna, at 345GHz
    # ax.axvline(x=1.7, color='r', linestyle='-.', alpha=0.5)
    # ax.text(1.7, 50, r'$39\,\mu$Jy', color='r', alpha=0.5)
    # ax.set_xlim(0, 3.7)
    # ax.set_ylim(0, 60)

