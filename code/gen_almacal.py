# It is a new program to read the exporture time of whole almacal observation

from astropy import table
from astropy.io import fits
import analysisUtils as au


# base_dir = '/science-ALMACAL/data'
# outputfile = '/home/jchen/Documents/almacal.new.info.txt'

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



base_dir = '/Users/jchen/Public/almacal_test/data'
outputfile = '/tmp/almacal.new.info.txt'

p_obj = re.compile('J\d+[+-]\d')
p_obs = re.compile('uid___')
bands = ['B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10']


bad_obs = '/Users/jchen/Desktop/projects/almacal/data/broken_obs.txt'

with open(bad_obs) as f:
    all_bad_obs = f.readlines()
for i in range(len(all_bad_obs)):
    all_bad_obs[i] = all_bad_obs[i].strip()


hdr = fits.Header()
hdr['OBSERVER'] = 'ALMACAL'
hdr['COMMENT'] = "Statistic information of the whole ALMACAL sample"
primary_hdu = fits.PrimaryHDU(header=hdr)
#data = np.diag([1,2,3,4])
#image_hdu = fits.ImageHDU(data, name="image")
c1 = fits.Column(name='a', array=np.array([1, 2]), format='K') # double int
c2 = fits.Column(name='b', array=np.array([4, 5]), format='K') # doubel int
c3 = fits.Column(name='c', array=np.array([7., 8.]), format='D') # double float
bintable_hdu = fits.BinTableHDU.from_columns([c1, c2, c3], name='bintable')
hdus = fits.HDUList([primary_hdu, bintable_hdu])
hdus.writeto('data/multitable.fits', overwrite=True)



# with open(outputfile, 'w+') as f_out:
    # # f_out.write('obj B3 B4 B5 B6 B7 B8 B9 B10\n')
    # f_out.write('{:<12s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}\n'.format('obj', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10'))


band_match = re.compile('_(?P<band>B\d{1,2})$')
obj_match = re.compile('J\d{4}[-+]\d{4}')

for obj in os.listdir(base_dir):
    obj_exptime = {'B3':0, 'B4':0, 'B5':0, 'B6':0, 
                'B7':0, 'B8':0,  'B9':0, 'B10':0}
    hdr = fits.Header()
    hdr['OBJ'] = obj
    hdr['COMMENT'] = "SPW and exposure time information"
    if not obj_match.match(obj):
        print('Error load obj:', obj)
        continue
    print("OBJ:", obj)
        
    for obs in os.listdir(base_dir +'/'+ obj):
        if obs in bad_obs:
            continue
        print("OBS:", obs)
        if p_obs.match(obs):
            if band_match.search(obs):
                band = band_match.search(obs).groupdict()['band']
            else:
                print('Error in match the band from', listobs)
                continue
            obs_filename = base_dir +'/'+ obj+'/'+obs
            try:
                time_on_source = au.timeOnSource(obs_filename, verbose=False, debug=False)
                time_minutes = time_on_source[0]['minutes_on_source']
                spw_specrange = read_spw(obs_filename)

            except:
                print("Error in read on source time:", obs_filename)
                continue
            if band in obj_exptime.keys():
                obj_exptime[band] = obj_exptime[band] + time_minutes


    column_list = []
    for band in bands: 
        hdr[band] = obj_exptime[band]
        column_list.append(fits.Column(name=band, array=**, format='D') # double int

    bintable_hdu = fits.BinTableHDU.from_columns(column_list, name='bintable')
    hdus = fits.HDUList([primary_hdu, bintable_hdu])

    with open(outputfile, 'a+') as f_out:
        f_out.write('{:<12s} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f}\n'.format(
            obj, 
            obj_exptime['B3'],\
            obj_exptime['B4'],\
            obj_exptime['B5'],\
            obj_exptime['B6'],\
            obj_exptime['B7'],\
            obj_exptime['B8'],\
            obj_exptime['B9'],\
            obj_exptime['B10'],\
            ))
 
