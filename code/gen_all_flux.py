# this is a quick program to gather all the flux information from the individual files.

import re
import os

base_dir = '/science-ALMACAL/temp/FLUXVAL'
outfile = './allcal.fluxval'

os.system('rm -rf {}'.format(outfile))
# obsname_match = re.compile('(?P<obsname>uid___\w*\.ms(\.split\.cal)?\.(?P<objname>J\d*[+-]+\d*)_(?P<band>B\d+))')
obsname_match = re.compile('(?P<obsname>uid___\w*\.ms(\.split\.cal)?\.(?P<objname>[\s\w+-]+)_(?P<band>B\d+))')
flux_match = re.compile('flux:\s+(?P<flux>\-?\d+(.\d+)?([Ee][+-]?\d+)?)\s+Jy')

#print("total:", len(os.listdir(base_dir)))
for obs in os.listdir(base_dir):
    #print(obs)
    try:
        fname = obsname_match.search(obs).groupdict()
    except:
        print("Error in match file: {}\n".format(obs))
        continue
    try:
        with open(os.path.join(base_dir, obs)) as f:
            flux_string = f.readline()
    except:
        print("Error in read file: {}".format(obs))
        continue
    #print(flux_string)
    objname = fname['objname']
    obsname = fname['obsname']
    band = fname['band']
    try:
        fflux = flux_match.search(flux_string).groupdict()
    except:
        print("Error in read the flux from: {}".format(obs))
        continue
    flux = fflux['flux']

    with open(outfile, 'a+') as fout:
        fout.write("{} {} {} {}\n".format(objname, obsname, band, flux))
    #print(obsname, band, flux)
