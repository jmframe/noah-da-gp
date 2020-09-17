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
gather_dir = parms_dir + 'gather_parms/'

# Read in list of sites
plum = pd.read_csv(proj_dir+'setup/plumber-2-sites.csv')
plum.set_index('site', inplace=True)

for i in range(plum.shape[0]):
    cmd = "mv site_data/parms_"+str(i+1)+".txt site_data/parms_"+plum.index.values[i]+".txt"
    os.system(cmd)
