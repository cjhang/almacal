#!/usr/bin/env python

# The function to read the listobs output into structured data

# Author: Jianhang Chen
# Email: cjhastro@gmail.com
# History:

import os
import re
from datetime import datetime, timedelta
import glob
import numpy as np


versions = '0.0.1'

def date_format(date):
    # format the date include abbriviation
    month_dict = {'Jan':'01', 'Feb':'02', 'Mar':'03', 'Apr':'04', 'May':'05', 
                  'Jun':'06', 'Jul':'07', 'Aug':'08', 'Sep':'09', 'Oct':'10',
                  'Nov':'11', 'Dec':'12'}
    for m in month_dict.keys():
        if m in date:
            date_new = date.replace(m, month_dict[m])

    return date_new
    

def read_obstime(listobs):
    # read the obstime information
    with open(listobs) as f:
        all_lines = f.readlines() 
    n_start = 0
    p_select = re.compile('Date\s+Timerange')

    for ind in range(len(all_lines)):
        p_matched = p_select.search(all_lines[ind])
        if p_matched:
            n_start = ind
            break
    match_validtime = re.compile('^\s+(\d+-\w+-\d{4}/){0,1}\d{2}:\d{2}:\d{2}.\d+')
    match_info = re.compile('(?P<Date>\d+-\w+-\d{4})\/')
    match_timerange = re.compile('(?P<starttime>\d{2}:\d{2}:\d{2}.\d+)\s-\s(?P<endtime>\d{2}:\d{2}:\d{2}.\d+)')
    match_scanintent = re.compile('\[(?P<ScanIntent>[A-Z_,#]+)\]')
    match_end = re.compile('^[A-Z]\w+')
    #obstime = {'scanintent':[], 'duration':0, 'freq':None}
    obstime = {'scanintent':[], 'duration':0, 'freq':None, 'startdate':None}
    for line in all_lines[n_start+1:]:
        if not match_validtime.match(line):
            continue
        if match_end.match(line):
            break
        if match_info.search(line):
            date = date_format(match_info.search(line).groupdict()['Date'])
            if obstime['startdate'] is None:
                obstime['startdate'] = date
            # print(date)
        
        # calculate the exposure time
        timerange = match_timerange.search(line).groupdict()
        starttime = timerange['starttime']
        endtime = timerange['endtime']
        starttime_iso = datetime.strptime(date+' '+starttime, '%d-%m-%Y %H:%M:%S.%f')
        endtime_iso = datetime.strptime(date+' '+endtime, '%d-%m-%Y %H:%M:%S.%f')
        if endtime_iso > starttime_iso:
            delta_time =  endtime_iso - starttime_iso
        else:
            delta_time = endtime_iso + timedelta(hours=24) - starttime_iso
        # print("diff time:", delta_time.total_seconds())
        obstime['duration'] = obstime['duration'] + delta_time.total_seconds()

        # find the role of the calibrator
        scanintent = match_scanintent.search(line).groupdict()['ScanIntent']
        obstime['scanintent'].append(scanintent.split(','))
    obstime['scanintent'] = np.unique(np.concatenate(obstime['scanintent'])).flatten().tolist()

    return obstime


def read_listobs(logfile):
    """
    """
    with open(logfile) as logfile:
        all_lines = logfile.readlines()

        # read the name and band in the filename
        name = logfile[-13:-3]
        band = logfile[-2:]

        # calculate the exposure time of each obs
        obstime = read_obstime(all_lines)

    return band, 


if __name__ == '__main__':
    listobs_folder = '/Users/jchen/Public/almacal_test/listobs'
    outputfile = '/tmp/almacal.new.info2.txt'
    print("Hello")
    with open(outputfile, 'w+') as f_out:
        f_out.write('obj B3 B4 B5 B6 B7 B8 B9 B10\n')
    # obj_scanintent = {'B3':[], 'B4':[], 'B5':[], 'B6':[], 
                      # 'B7':[], 'B8':[], 'B9':[], 'B10':[]}

    for obj in os.listdir(listobs_folder):
        obj_exptime = {'B3':0, 'B4':0, 'B5':0, 'B6':0, 
                    'B7':0, 'B8':0,  'B9':0, 'B10':0}
        obj_match = re.compile('J\d{4}[-+]\d{4}')
        if not obj_match.match(obj):
            print('Error load obj:', obj)
            continue
        for listobs in glob.glob(listobs_folder + '/{}/*.listobs.txt'.format(obj)):
            #print(listobs)
            listobs_basename = os.path.basename(listobs)
            band_match = re.compile('(?P<band>B\d{1,2})\.listobs\.txt')
            if band_match.search(listobs_basename):
                band = band_match.search(listobs_basename).groupdict()['band']
            else:
                print('Error in match the band from', listobs)
                continue
            if band in obj_exptime.keys():
                # with open(listobs) as f:
                    # all_lines = f.readlines() 
                obstime = read_obstime(listobs)
                # print('obstime:', obstime)
                
                #if obstime['startdate'] > '01-07-2015':
                #    continue

                obj_exptime[band] = obj_exptime[band] + obstime['duration']
                # obj_scanintent[band].append(obstime['scanintent'])
            
            # for band in obj_scanintent.keys():
                # if len(obj_scanintent[band])>0:
                    # obj_scanintent[band] = ','.join(np.unique(np.concatenate(obj_scanintent[band])))
            #out_string = '{} '.format(obj)
            #for key in obj_exptime:
            #    out_string += ' {} {}'.format(obj_exptime[key])
        with open(outputfile, 'a+') as f_out:
            f_out.write('{} {} {} {} {} {} {} {} {}\n'.format(
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
        


