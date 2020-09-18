#!/discover/nobackup/jframe/anaconda3/bin/python

import numpy as np
import csv
import sys
import datetime as dt
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

plot_variable = sys.argv[1]

#########  Set up some lists for the data   ################
run_dir = '/discover/nobackup/jframe/noah-da-gp/runs/ZM-Mon/'
obs_dict = {'NEE':5, 'GPP':6, 'Qle':7, 'Qh':8}
out_dict = {'Qle':24, 'Qh':25, 'NEE':26}
x = []
obs = []
noah = []
p = []

########    import all the data   ########################3## 
with open(run_dir+'obs.txt', 'r') as f:
    d = csv.reader(f, delimiter=' ', skipinitialspace=True)
    for obsRow in d:
        obs.append(float(obsRow[obs_dict[plot_variable]]))
    obs = np.array(obs)
with open(run_dir+'output.noah', 'r') as f:
    d = csv.reader(f, delimiter=' ', skipinitialspace=True )
    for row in d:
        D = dt.datetime(int(row[0]),int(row[1]),int(row[2]),int(row[3]),int(row[4]))
        x.append(D)
        p.append(float(row[10]))
        noah.append(float(row[out_dict[plot_variable]]))
    noah = np.array(noah)
    x = np.array(x)

print('plotting')
fig, ax1 = plt.subplots(figsize=(15,10))
ax1.plot(x,obs, label='Observation', color='k', linewidth=0.5)
ax1.plot(x,noah, label='Noah prediction', color='b', linewidth=0.5)
ax1.set_ylabel(plot_variable)
ax1.set_xlabel('Days')
plt.legend(loc='upper left')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
