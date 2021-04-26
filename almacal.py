# The main class file for ALMACAL project

import os
import sys
import re


class ALMACAL(object):
    """The Project class
    """
    def __init__(self, basedir=None):
        self.basedir = basedir
        self.objlist = self.get_objlist()

    def get_objlist(self):
        p_obj = re.compile('J\d+[+-]\d')
        objlist = []
        for item in os.listdir(self.basedir):
            if p_obj.match(item):
                objlist.append(os.path.join(self.basedir, item))
        return objlist
        

class Calibrator(object):
    """The Calibrator class
    """

    def __init__(self, basedir=None):
        self.basedir = basedir
        self.obslist = self.get_obslist()
        
    def get_obslist(self):
        p_obs = re.compile('uid___')
        obslist = []
        for item in os.listdir(self.basedir):
            if p_obs.match(item):
                obslist.append(os.path.join(self.basedir, item))
        return obslist

    def get_obs(self, idx):
        return SingleObs(self.obslist[idx])
        

class SingleObs(object):
    """The single observation
    """
    def __init__(self, basedir=None):
        self.basedir = basedir
        self.name = os.path.basename(basedir)

