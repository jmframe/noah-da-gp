# Coded on 2020/11/20 by Soni Yatheendradas 

#!/discover/nobackup/jframe/anaconda3/bin/python
# This is a script to do the actual OL runs for
# noah-mp with the PLUMBER-2 data.


# This is a script to generate sbatch PoDS and corresponding
#   Exec scripts for noah-mp with the PLUMBER-2 data.

import numpy as np
import os
import pandas as pd

# BEGIN any specifications/inputs required for script to run

proj_dir = '/discover/nobackup/syatheen/Grey/THP/noah-da-gp/'

PoDS_FileName = 'OL_sbatchpods1.sh'

nTasksPerNode = 14

Exec_FileName = 'Execfile_OL'  

# END any specifications/inputs required for script to run

# Begin delete of existing scripts
if os.path.exists(PoDS_FileName):
  cmd = '/bin/rm ' + PoDS_FileName
  os.system(cmd)
if os.path.exists(Exec_FileName):
  cmd = '/bin/rm ' + Exec_FileName
  os.system(cmd)
# End delete of existing scripts

# Read in list of sites
plum = pd.read_csv(proj_dir+'setup/plumber-2-sites.csv')
plum.set_index('site', inplace=True)

Num_Sites = (plum.index.values).shape[0]

# copy from PoDS template script file to create PoDS script file 
cmd = 'cp OL_sbatchpods1_template.sh ' + PoDS_FileName
os.system(cmd)

# Calculate ceiling of Num_Sites_TimeRanges divided by nTasksPerNode, and 
#  corresponding total ntasks, for the PoDS file
nNodes_Allocated = -np.floor_divide(-Num_Sites, nTasksPerNode)
Total_nTasks_Allocated = nNodes_Allocated * nTasksPerNode

# Replace text in PoDS script file
cmd = 'sed -i "s/nttasks/' + str(Total_nTasks_Allocated) + '/g" ' + PoDS_FileName
os.system(cmd)
cmd = 'sed -i "s/nTasksPerNode/' + str(nTasksPerNode) + '/g" ' + PoDS_FileName
os.system(cmd)

# Loop through the sites
for s in list(plum.index.values):

  # write information of each site number and time range as each line
  cmd = 'echo ./RunProcess_OL.sh ' + s + ' >> ' + Exec_FileName
  os.system(cmd)

## --- End Script ---------------------------------------------------

