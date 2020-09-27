# This is the file trying to plot the difference between the data up to now with that of 
# the Iran Oteo in 2016
import numpy as np
import matplotlib.pyplot as plt

almacal_info_file_new = '../data/almacal.new.info.txt'
data_new = np.loadtxt(almacal_info_file_new, skiprows=1, delimiter=' ',
        dtype={'names': ('obj', 'B3','B4','B5','B6','B7','B8','B9'), 
               'formats': ('S10', 'f4', 'f4','f4','f4','f4','f4','f4')})

almacal_info_file_oteo2016 = '../data/almacal.oteo2016.info.txt'
data_oteo2016 = np.loadtxt(almacal_info_file_oteo2016, skiprows=1, delimiter=' ',
        dtype={'names': ('obj', 'B3','B4','B5','B6','B7','B8','B9'), 
               'formats': ('S10', 'f4', 'f4','f4','f4','f4','f4','f4')})

data = data_oteo2016
save_fig_name = 'data_oteo2016'

bands = ['B6', 'B7']

# for the data used in oteo2016
if True:
    band_depth = {'B6':[], 'B7':[]}
    for row in data:
        for band in bands:
            obj = (row['obj']).decode('ascii')
            band_depth[band].append(row[band]/60) # in minutes
    for band, depth in band_depth.items():
        obs_valid = np.array(depth)>0.1
        band_depth[band] = np.array(depth)[obs_valid]

if True:

    fig = plt.figure(figsize=(12,4))
    fig.suptitle('Data Oteo2016', fontsize=16)

    ax = fig.add_subplot(121)
    ax.set_title('Band 6')
    im = plt.hist(np.log10(band_depth['B6']), bins=20)
    ax.set_xlabel(r'$\log (t_{\rm obs}[Minutes])$')
    ax.set_ylabel('Counts')

    ax = fig.add_subplot(122)
    ax.set_title('Band 7')
    im = plt.hist(np.log10(band_depth['B7']), bins=20)
    ax.set_xlabel(r'$\log (t_{\rm obs}[Minutes])$')
    ax.set_ylabel('Counts')
    # for 3000 minutes
    ax.axvline(x=3.5, color='r', linestyle='-.', alpha=0.5)
    ax.text(3.5, 5, r'$3.5\mu$Jy', color='r', alpha=0.5)

    plt.show()
    fig.savefig('results/{}.pdf'.format(save_fig_name))
