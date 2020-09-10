#!/discover/nobackup/jframe/anaconda3/bin/python
# This is a script to generate parameters for 
# Running noah-mp with the PLUMBER-2 data.

import pandas as pd
import numpy as np
import shutil, errno
import os
import os.path
from os import path
import sys
from io import StringIO

proj_dir = '/discover/nobackup/jframe/noah-da-gp/'
parms_dir = proj_dir + 'setup/parms/'
plum_parms_dir = parms_dir + 'plumber-2-parms/'
extract_dir = parms_dir + 'extract_parms/'

# Read in list of sites
plum = pd.read_csv(proj_dir+'setup/plumber-2-sites.csv')
plum.set_index('site', inplace=True)
  
for s in list(plum.index.values):
    # screen report
    print('Making parameter files for site {}'.format(s))
    stryr, endyr, source, lat, lon = plum.loc[s,:]
    yrs = str(stryr)+'-'+str(endyr)
    print(stryr, endyr, source, lat, lon)

    # Copy the namelist templace and replace with site specific lat/lon/start/end
    startdate = str(stryr)+"01010000"
    enddate = str(endyr+1)+"01010000"
    startdate = "199801010900"
    enddate = "199901020000"
    namelist = "namelist_"+s
    cmd = "cp namelist_template "+namelist
    os.system(cmd)
    cmd = "sed -i 's/_startdate_/"+startdate+"/g' "+namelist 
    os.system(cmd)
    cmd = "sed -i 's/_enddate_/"+enddate+"/g' "+namelist
    os.system(cmd)
    cmd = "sed -i 's/_lat_/"+str(lat)+"/g' "+namelist 
    os.system(cmd)
    cmd = "sed -i 's/_lon_/"+str(lon)+"/g' "+namelist
    os.system(cmd)
    
    # Now run the extraction script
    cmd = "./driver.exe "+namelist
    os.system(cmd)
    cmd = "mv parms.out ../plumber-2-parms/parms_"+s+".txt"
    os.system(cmd)
    cmd = "mv cal_parms.out ../plumber-2-parms/cal_parms_"+s+".txt"
    os.system(cmd)
    cmd = "mv time_parms.txt ../plumber-2-parms/time_parms_"+s+".txt"
    os.system(cmd)
    cmd = "mv OUTPUT.txt ../plumber-2-parms/OUTPUT_"+s+".txt"
    os.system(cmd)
    cmd = "mv state_save.txt ../plumber-2-parms/state_save_"+s+".txt"
    os.system(cmd)

    # Lastly remove the site namelist file
    cmd = "mv "+namelist+" ../../../z_trash/"
    os.system(cmd)

# And just for good measure, call the extraction with the original namelist
cmd = "./driver.exe namelist"
os.system(cmd)

