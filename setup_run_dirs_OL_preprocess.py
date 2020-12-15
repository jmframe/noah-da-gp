#!/discover/nobackup/jframe/anaconda3/bin/python
# This is a script to preprocess before 
# Running noah-mp with OL for the PLUMBER-2 data.
# To be run after setup_run_dirs.py and setup_run_dirs_extra.py

import pandas as pd
import numpy as np
import os
from os import path
import sys

###BEGIN ANY EDITS HERE######

proj_dir = '/discover/nobackup/syatheen/Grey/THP/noah-da-gp/'

###END ANY EDITS HERE######

# Since this program will override the initialization files, 
# be careful running it. Check to make sure you meant to run.
print('this preprocess is to be done only after running setup_run_dirs.py .')
continue_setup = str(input('Are you sure that you already ran setup_run_dirs.py ?'))
if not continue_setup in ['yes', 'Yes', 'Y', 'y', 'continue', 'Continue', 'affirmative']:
  print('Exiting program based on your answer.')
  sys.exit()

# Read in list of sites
plum = pd.read_csv(proj_dir+'setup/plumber-2-sites.csv')
plum.set_index('site', inplace=True)
  
# --- Set Up Preprocessing Inside Already Created Runtime Directories ------------------------------------
for s in list(plum.index.values):

  # screen report
  print('Setting up preprocessing in run directory for site {}'.format(s))
  stryr, endyr, source, lat, lon = plum.loc[s,:]
  yrs = str(stryr)+'-'+str(endyr)
  print(stryr, endyr, source, lat, lon)

  # working directory
  wdir = proj_dir+'runs/'+s

  # generic noah-mp input files
  cmd = 'cp setup/da_flag_M99.txt ' + wdir + '/da_flag.txt'
  os.system(cmd)

  # TEMPORARY, just getting sonething that will run. NEED TO EXTRACT REAL VALUES!!!
  cmd = 'cp setup/soil_init.txt ' + wdir + '/soil_init.txt'
  os.system(cmd)
  cmd = 'cp setup/plant_init.txt ' + wdir + '/plant_init.txt'
  os.system(cmd)
  cmd = 'cp setup/soil_init.txt ' + wdir + '/soil_init.txt_StartingPoint'
  os.system(cmd)
  cmd = 'cp setup/plant_init.txt ' + wdir + '/plant_init.txt_StartingPoint'
  os.system(cmd)

  # get site information
  forcing = np.genfromtxt(wdir + '/forcing.txt')

  # site-specific noah-mp input files - startdate
  startdate = str( int(round(forcing[0,0]))*100000000 + int(round(forcing[0,1]))*1000000 + int(round(forcing[0,2]))*10000 + int(round(forcing[0,3]))*100 + int(round(forcing[0,4])) )
  fname = wdir + '/startdate.txt'
  with open(fname,'w') as FF:
    FF.write(startdate)

  if not path.exists(wdir + '/EnsembleOut/'):
    cmd = 'mkdir -p ' + wdir + '/EnsembleOut/'
    os.system(cmd)

  # Output specifications
  cmd = 'cp setup/output_specifications_OL.txt ' + wdir + '/output_specifications.txt'
  os.system(cmd)

  # sig_*txt files
  cmd = 'cp setup/sig_sm.txt ' + wdir + '/sig_sm.txt'
  os.system(cmd)
  cmd = 'cp setup/sig_veg.txt ' + wdir + '/sig_veg.txt'
  os.system(cmd)

#end of for s in list(plum.index.values)
# --- End Script ---------------------------------------------------
