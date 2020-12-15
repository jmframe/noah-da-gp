#!/discover/nobackup/jframe/anaconda3/bin/python
# This is a script to link forcing and obs files that
# didn't get done properly from other user accounts 
# besides Jonathan's, after running setup_run_dirs.py

import pandas as pd
import os

###BEGIN ANY EDITS HERE######

proj_dir = '/discover/nobackup/syatheen/Grey/THP/noah-da-gp/' # Give your working directory path

###END ANY EDITS HERE######

# Read in list of sites
plum = pd.read_csv(proj_dir+'setup/plumber-2-sites.csv')
plum.set_index('site', inplace=True)
  
# --- Set Up Inside Runtime Directories ------------------------------------
for s in list(plum.index.values):

  # screen report
  print('Extra setting up run directory for site {}'.format(s))
  stryr, endyr, source, lat, lon = plum.loc[s,:]
  yrs = str(stryr)+'-'+str(endyr)
  print(stryr, endyr, source, lat, lon)

  # working directory
  wdir = proj_dir+'runs/'+s

  # delete forcing and obs files having ineffective links
  cmd = 'rm -rf ' + wdir + '/forcing.txt'
  os.system(cmd)
  cmd = 'rm -rf ' + wdir + '/obs.txt'
  os.system(cmd)

  # Large data files. Forcing and observation.
  tup = (s,yrs,source,'Met.txt')
  met_file_name = '_'.join(tup)
  print(met_file_name)
  cmd = 'ln -s /discover/nobackup/jframe/data/plumber-2/met-txt/' + met_file_name + ' '+wdir+ '/forcing.txt'
  os.system(cmd)
  tup = (s,yrs,source,'Flux.txt')
  flux_file_name = '_'.join(tup)
  cmd = 'ln -s /discover/nobackup/jframe/data/plumber-2/flux-txt/' + flux_file_name + ' '+wdir+ '/obs.txt'
  os.system(cmd)

# --- End Script ---------------------------------------------------
