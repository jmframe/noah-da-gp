#!/discover/nobackup/jframe/anaconda3/bin/python
# This is a script to generate directories for 
# Running noah-mp with the PLUMBER-2 data.
# Later on will make specific setup files for DA and GP
# For no just making sure we can run the model with the data

import pandas as pd
import numpy as np
import shutil, errno
import os
import os.path
from os import path
import sys
from io import StringIO

proj_dir = '/discover/nobackup/jframe/noah-da-gp/'

# Since this program will override the initialization files, 
# be careful running it. Check to make sure you meant to run.
print('this setup will override existing run files!')
continue_setup = str(input('Are you sure?'))
if not continue_setup in ['yes', 'Yes', 'Y', 'y', 'continue', 'Continue', 'affirmative']:
    print('Exiting program based on your answer.')
    sys.exit()

print('is this the correct flux directory: /discover/nobackup/data/plumber-2-flux-txt/ ?')
print('is this the correct met directory: /discover/nobackup/data/plumber-2-met-txt/ ?')
is_flux_dir = str(input('/discover/nobackup/data/plumber-2-flux-txt/?'))
if not is_flux_dir in ['yes', 'Yes', 'Y', 'y', 'continue', 'Continue', 'affirmative']:
    data_dir = str(input('type in the data directory?'))
    sys.exit()
else:
    flx_dir = '/discover/nobackup/jframe/data/plumber-2-flux-txt/'
    met_dir = '/discover/nobackup/jframe/data/plumber-2-met-txt/'

# Read in list of sites
plum = pd.read_csv(proj_dir+'setup/plumber-2-sites.csv')
plum.set_index('site', inplace=True)
  
# --- Set Up Runtime Directories ------------------------------------
# delete and reinit test dirs
if not path.exists('runs/'):
    cmd = 'mkdir runs/'

for s in list(plum.index.values):
    # screen report
    print('Setting up run directory for site {}'.format(s))
    stryr, endyr, source, lat, lon = plum.loc[s,:]
    yrs = str(stryr)+'-'+str(endyr)
    print(stryr, endyr, source, lat, lon)

    # working directory
    wdir = proj_dir+'runs/'+s

    # delete runtime directory for the site
    cmd = 'rm -rf ' + wdir
    os.system(cmd)
    cmd = 'mkdir ' + wdir
    os.system(cmd)
    cmd = 'mkdir ' + wdir + '/reports'
    os.system(cmd)

    # executables
    cmd = 'ln -s ../../noah-mp/noah_mp.exe ' + wdir + '/noah_mp.exe'  
    os.system(cmd)
    cmd = 'ln -s ../../setup/periodic_cleanup.sh ' + wdir + '/periodic_cleanup.sh'  
    os.system(cmd)

    cmd = 'cp setup/cal_parms.tpl ' + wdir
    os.system(cmd) 

    # generic noah-mp input files
    cmd = 'cp setup/da_flag_0.txt ' + wdir + '/da_flag.txt'
    os.system(cmd)
    cmd = 'cp setup/init.txt ' + wdir
    os.system(cmd)
    cmd = 'cp setup/tbot.txt ' + wdir
    os.system(cmd)

    # Large data files. Forcing and observation.
    tup = (s,yrs,source,'Met.txt')
    met_file_name = '_'.join(tup)
    print(met_file_name)
    cmd = 'ln -s ../../../data/plumber-2/met-txt/' + met_file_name + ' '+wdir+ '/forcing.txt'
    os.system(cmd)
    tup = (s,yrs,source,'Flux.txt')
    flux_file_name = '_'.join(tup)
    cmd = 'ln -s ../../../data/plumber-2/flux-txt/' + flux_file_name + ' '+wdir+ '/obs.txt'
    os.system(cmd)

#    # site-specific noah-mp input files
#    ignore_these_parms = [3,4,5,23,25,26,27,30,34,45,46,52,53,54]
#    P=[]
#    with open(data_dir+'parms/extract_parms/site_data/parms_' + S + '.txt', 'r') as f:
#        for i, x in enumerate(f):
#            if i not in ignore_these_parms:
#                P.append(x)
#    with open(wdir+'parms.txt', mode='wt', encoding='utf-8') as f:
#          f.writelines("%s" % p for p in P)
#    cmd = 'cp '+data_dir+'parms/extract_parms/site_data/cal_parms_' + S + '.txt ' + wdir + '/cal_parms.txt'  
#    os.system(cmd)
#    cmd = 'cp '+data_dir+'parms/extract_parms/site_data/time_parms_' + S + '.txt ' + wdir + '/time_parms.txt'
#    os.system(cmd)
#    offset = pals_sites[int(sites[s,0])-1,2]
    
    # TEMPORARY, just getting sonething that will run. NEED TO EXTRACT REAL VALUES!!!
    offset=0
    cmd = 'cp setup/parms/plumber-2-parms/time_parms_'+s+'.txt ' + wdir + '/time_parms.txt'
    os.system(cmd)
    cmd = 'cp setup/parms/plumber-2-parms/cal_parms_'+s+'.txt ' + wdir + '/cal_parms.txt'
    os.system(cmd)
    cmd = 'cp setup/parms/plumber-2-parms/parms_'+s+'.txt ' + wdir + '/parms.txt'
    os.system(cmd)
    cmd = 'cp setup/soil_init.txt ' + wdir + '/soil_init.txt'
    os.system(cmd)
    cmd = 'cp setup/plant_init.txt ' + wdir + '/plant_init.txt'
    os.system(cmd)

    # get site information
    #forcing = np.genfromtxt(wdir + '/forcing.txt',skip_header=1)
    forcing = pd.read_csv(wdir + '/forcing.txt')
    forcing = np.array(forcing)

    # site-specific noah-mp input files - number of timesteps
    Nt = forcing.shape[0]
    fname = wdir + '/num_times.txt'
    with np.printoptions(precision=0, suppress=True):
        with open(fname, 'w') as F:
            #################################################
            ######## Might wanst to calibrate to a shorter period of time
            ######## If so, just change this to the number of time steps.
            F.write(str(Nt))

    # site-specific noah-mp input files - lat/lon
    fname = wdir + '/lat_lon.txt'
    with np.printoptions(precision=7, suppress=True):
        with open(fname, 'w') as F:
            F.write('%f\n%f' % (lat, lon)) 


    # site-specific noah-mp input files - startdate
    if offset < 0:
        startdate = 200012311200 + (12+offset)*100  # offset here is negative to the west
    else:
        startdate = 200101011200 - (12-offset)*100

    startdate = str(int(startdate));
    fname = wdir + '/startdate.txt'
    with open(fname,'w') as F:
        F.write(startdate)
# --- End Script ---------------------------------------------------
