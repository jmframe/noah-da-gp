#!/discover/nobackup/jframe/anaconda3/bin/python
# This is a script to generate directories for 
# Running noah-mp with the PLUMBER-2 data.
# Later on will make specific setup files for DA and GP
# For no just making sure we can run the model with the data

import pandas as pd
import numpy as np
import shutil, errno
import os
from os.path import dirname as up
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
site_list = pd.read_csv(proj_dir+'setup/plumber-2-sites.csv')
print(site_list)
exit()

# --- Set Up Runtime Directories ------------------------------------
# delete and reinit test dirs
cmd = 'rm -rf ' + exp_type + '/cal_dirs'
os.system(cmd)
cmd = 'mkdir ' + exp_type + '/cal_dirs'
os.system(cmd)
cmd = 'mkdir ' + exp_type + '/cal_dirs/reports'
os.system(cmd)
cmd = 'cp setup_dir/cleanCalFiles.py ' + exp_type + '/cal _dirs/'

for s in range(0, Ns):

    # get sites and year, need in string
    S = str(int(sites[s,0]))
    Y = str(int(sites[s,1]))
    siteyear = S + '_' + Y

    # screen report
    print('Setting up run directory for site = %d, year = %d, ' % (sites[s,0],  sites[s,1]))
    print('siteyear: ', siteyear)

    # working directory
    wdir = exp_type + '/cal_dirs/run_' + siteyear

    # delete runtime directory for the site
    cmd = 'rm -rf ' + wdir
    os.system(cmd)
    cmd = 'mkdir ' + wdir
    os.system(cmd)
    cmd = 'mkdir ' + wdir + '/reports'
    os.system(cmd)

    # Make a file with the site and year, 
    # since this information isn't in any of the files.
    fname = wdir + '/site.txt'
    with np.printoptions(precision=7, suppress=True):
        with open(fname, 'w') as F:
            F.write(S)
    fname = wdir + '/year.txt'

    # executables
    cmd = 'ln -s ../../../ostrich/Source/Ostrich ' + wdir + '/ostrich.exe'  
    os.system(cmd)
    cmd = 'ln -s ../../../model_code/noah_mp.exe ' + wdir + '/noah_mp.exe'  
    os.system(cmd)
    cmd = 'ln -s ../../../setup_dir/periodic_cleanup.sh ' + wdir + '/periodic_cleanup.sh'  
    os.system(cmd)

    cmd = 'cp setup_dir/ostIn.txt ' + wdir + '/ostIn.txt' 
    os.system(cmd) 
    cmd = 'cp setup_dir/cal_parms.tpl ' + wdir
    os.system(cmd) 

    # generic noah-mp input files
    cmd = 'cp ./setup_dir/da_flag_cal.txt ' + wdir + '/da_flag.txt'
    os.system(cmd)
    cmd = 'cp ./setup_dir/init.txt ' + wdir
    os.system(cmd)
    cmd = 'cp ./setup_dir/tbot.txt ' + wdir
    os.system(cmd)

    # Large data files. Forcing and observation.
    cmd = 'ln -s ../../../data/pals/site_data/forcing_' + S +'.txt '+ wdir + '/forcing.txt'  
    os.system(cmd)
    cmd = 'ln -s ../../../data/pals/site_data/obs_' + S +'_'+ Y + '_train.txt ' + wdir + '/obs.txt'  
    os.system(cmd)

    # site-specific noah-mp input files
    ignore_these_parms = [3,4,5,23,25,26,27,30,34,45,46,52,53,54]
    P=[]
    with open(data_dir+'parms/extract_parms/site_data/parms_' + S + '.txt', 'r') as f:
        for i, x in enumerate(f):
            if i not in ignore_these_parms:
                P.append(x)
    with open(wdir+'parms.txt', mode='wt', encoding='utf-8') as f:
          f.writelines("%s" % p for p in P)

    cmd = 'cp '+data_dir+'parms/extract_parms/site_data/cal_parms_' + S + '.txt ' + wdir + '/cal_parms.txt'  
    os.system(cmd)
    cmd = 'cp '+data_dir+'parms/extract_parms/site_data/time_parms_' + S + '.txt ' + wdir + '/time_parms.txt'
    os.system(cmd)

    # get site information
    forcing = np.genfromtxt(wdir + '/forcing.txt')
    pals_sites = np.genfromtxt('data/pals/Sites.txt')
    lat = pals_sites[int(sites[s,0])-1,0]
    lon = pals_sites[int(sites[s,0])-1,1]
    offset = pals_sites[int(sites[s,0])-1,2]

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

    # site-specific noah-mp input files - initial states
    cmd = 'cp ../gpr_fluxnet/initialize/site_data/soil_init_' + S + '.txt ' + wdir + '/soil_init.txt'
    os.system(cmd)
# --- End Script ---------------------------------------------------
