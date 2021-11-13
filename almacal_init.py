import os
import sys

# ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
if 'ALMACAL_NUMBERCOUNTS_HOME' in os.environ.keys():
    ROOT_DIR = os.environ['ALMACAL_NUMBERCOUNTS_HOME']
else:
    ROOT_DIR = os.path.join(os.path.expanduser('~'), 'work/projects/almacal/number_counts')
ROOT_DIR_CODE = os.path.join(ROOT_DIR, 'code')
print("Project path: {}".format(ROOT_DIR))
sys.path.append(ROOT_DIR)

execfile(os.path.join(ROOT_DIR_CODE, 'almacal.py'))
execfile(os.path.join(ROOT_DIR_CODE, 'imaging_utils.py'))
execfile(os.path.join(ROOT_DIR_CODE, 'simalma_utils.py'))
execfile(os.path.join(ROOT_DIR_CODE, 'readms.py'))
execfile(os.path.join(ROOT_DIR_CODE, 'almacal_run.py'))
