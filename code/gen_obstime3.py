# It is a new program to read the exporture time of whole almacal observation

import analysisUtils as au


# base_dir = '/science-ALMACAL/data'
# output_dir = '/home/jchen/Documents/almacal_stat'

base_dir = '/Users/jchen/Public/almacal_test/data'
output_dir = '/tmp'

p_obj = re.compile('J\d+[+-]\d')
p_obs = re.compile('uid___')

bad_obs = '/Users/jchen/Desktop/projects/almacal/data/broken_obs.txt'

with open(bad_obs) as f:
    all_bad_obs = f.readlines()
for i in range(len(all_bad_obs)):
    all_bad_obs[i] = all_bad_obs[i].strip()

band_match = re.compile('_(?P<band>B\d{1,2})$')
obj_match = re.compile('J\d{4}[-+]\d{4}')

for obj in os.listdir(base_dir):
    obj_exptime = {'B3':0, 'B4':0, 'B5':0, 'B6':0, 
                'B7':0, 'B8':0,  'B9':0, 'B10':0}
    if not obj_match.match(obj):
        print('Error load obj:', obj)
        continue
    print("OBJ:", obj)
        
    obj_dirname = base_dir +'/'+ obj
    obj_output_dir = output_dir + '/' + obj
    os.system('mkdir {}'.format(obj_output_dir))
    spw_stat(obj_dirname, plot=True, showfig=False, \
            figname=obj_output_dir+'/'+obj+'.pdf', \
            plotbands=['B5','B6','B7','B8'], 
            savedata=True, filename=obj_output_dir+'/'+ obj+'.fits')

