# It is a new program to read the exporture time of whole almacal observation

import analysisUtils as au


# base_dir = '/science-ALMACAL/data'
# outputfile = '/home/jchen/Documents/almacal.new.info.txt'

base_dir = '/Users/jchen/Public/almacal_test/data'
outputfile = '/tmp/almacal.new.info.txt'

p_obj = re.compile('J\d+[+-]\d')
p_obs = re.compile('uid___')

bad_obs = '/Users/jchen/Desktop/projects/almacal/data/broken_obs.txt'

with open(bad_obs) as f:
    all_bad_obs = f.readlines()
for i in range(len(all_bad_obs)):
    all_bad_obs[i] = all_bad_obs[i].strip()

with open(outputfile, 'w+') as f_out:
    # f_out.write('obj B3 B4 B5 B6 B7 B8 B9 B10\n')
    f_out.write('{:<12s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}\n'.format('obj', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10'))


band_match = re.compile('_(?P<band>B\d{1,2})$')
obj_match = re.compile('J\d{4}[-+]\d{4}')

for obj in os.listdir(base_dir):
    obj_exptime = {'B3':0, 'B4':0, 'B5':0, 'B6':0, 
                'B7':0, 'B8':0,  'B9':0, 'B10':0}
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

            except:
                print("Error: in", obs_filename)
            if band in obj_exptime.keys():
                obj_exptime[band] = obj_exptime[band] + time_minutes
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
 
