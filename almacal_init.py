import os
import sys

if 'ALMACAL_NUMBERCOUNTS_HOME' in os.environ.keys():
    root_path = os.environ['ALMACAL_NUMBERCOUNTS_HOME']
else:
    root_path = os.path.join(os.path.expanduser('~'), 'projects/almacal/number_counts')
print("Project path: {}".format(root_path))

sys.path.append(root_path+'/code')

execfile(root_path+'/code/almacal.py')
execfile(root_path+'/code/imaging_utils.py')
execfile(root_path+'/code/simalma_utils.py')
execfile(root_path+'/code/readms.py')
execfile(root_path+'/code/almacal_run.py')
