# This is the file trying to plot the difference between the data up to now with that of 
# the Iran Oteo in 2016
import numpy as np
import matplotlib.pyplot as plt

almacal_info_file_new = '../data/almacal.new.info.txt'
data_new = np.loadtxt(almacal_info_file_new, skiprows=1, delimiter=' ',
        dtype={'names': ('obj', 'B3','B4','B5','B6','B7','B8','B9'), 
               'formats': ('S10', 'f4', 'f4','f4','f4','f4','f4','f4')})

almacal_info_file_good = '../data/make_all_goodimages.txt'
data_good = np.loadtxt(almacal_info_file_good, skiprows=0,
        dtype={'names': ('obj', 'B3','B4','B5','B6','B7','B8','B9','B10'), 
               'formats': ('S10', 'f4', 'f4','f4','f4','f4','f4','f4','f4')})

almacal_info_file_oteo2016 = '../data/almacal.oteo2016.info.txt'
data_oteo2016 = np.loadtxt(almacal_info_file_oteo2016, skiprows=1, delimiter=' ',
        dtype={'names': ('obj', 'B3','B4','B5','B6','B7','B8','B9'), 
               'formats': ('S10', 'f4', 'f4','f4','f4','f4','f4','f4')})

data = data_oteo2016
save_fig_name = 'data_oteo2016'
data = data_new
save_fig_name = 'data_new'

bands = ['B6', 'B7']

# for the data used in oteo2016
if True:
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

if True:

    fig = plt.figure(figsize=(11,5))
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
    ax = fig.add_subplot(121)
    # ax.set_xscale('log')
    ax.set_title('Band 6')
    # ax.text(3.2, 87, 'Band 6')
    im = plt.hist(np.log10(band_depth_new['B6']), bins=20, label='Whole dataset before 2020', alpha=0.5)
    im = plt.hist(np.log10(band_depth_good['B6']), bins=20, label='Usable dataset before 2020', alpha=0.7)
    im = plt.hist(np.log10(band_depth_oteo2016['B6']), bins=20, label='Oteo et al. (2016)', alpha=0.9)
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

    ax = fig.add_subplot(122)
    # ax.text(3.2, 55, 'Band 7')
    ax.set_title('Band 7')
    im = plt.hist(np.log10(band_depth_new['B7']), bins=20, alpha=0.5)
    im = plt.hist(np.log10(band_depth_good['B7']), bins=20, alpha=0.7)
    im = plt.hist(np.log10(band_depth_oteo2016['B7']), bins=20, alpha=0.9)
    plt.xticks([0, 1, 2, 3], [r'$1$', r'$10$', r'$10^2$', r'$10^3$'])
    ax.set_xlabel(r'$\log (t_{\rm obs}[{\rm Minutes}])$')
    ax.set_ylabel('Number')
    # for 1000 minutes, 43 antenna, at 345 GHz
    ax.axvline(x=3, color='r', linestyle='-.', alpha=0.5)
    ax.text(3.1, 50, r'$6\,\mu$Jy', color='r', alpha=0.5)
    # for 50 minutes, 30 antenna, at 345GHz
    ax.axvline(x=1.7, color='r', linestyle='-.', alpha=0.5)
    ax.text(1.7, 50, r'$39\,\mu$Jy', color='r', alpha=0.5)
    ax.set_xlim(0, 3.7)
    ax.set_ylim(0, 60)

    plt.show()
    #fig.savefig('../results/comparison.pdf')
