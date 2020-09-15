#!/gpfsm/dulocal/sles11/other/SLES11.3/miniconda3/2019.03_py3.7/2019-05-15/bin/python

#Plotting the results of a noah-mp run with calibrated parameters
#Observation vs Model Output

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import csv
import datetime as dt
import os

x = []
xDays = []
y_obs = []
y_mod = []
p = []
count_row = 0
xDays.append(0)
with open('output.out','r') as csvfile:
    model = csv.reader(csvfile, delimiter=' ', skipinitialspace=True )
    for row in model:
        p.append(float(row[9]))
        y_mod.append(float(row[10]))
    simLength = len(y_mod)
with open('obs.txt','r') as csvfile:
    observations = csv.reader(csvfile, delimiter=' ', skipinitialspace=True )
    for row in observations:
        if count_row < simLength:	
            D = dt.datetime(int(row[0]),1,1) + \
                dt.timedelta((int(row[1])-1) + float(row[2])/24)
            x.append(D)
            y_obs.append(float(row[5]))
            if count_row > 0:
                xDays.append(xDays[count_row-1] + 1/48)
            count_row = count_row + 1

sYear = x[0].year
sMonth = x[0].month
sDay = x[0].day

startPlot = dt.date(sYear,sMonth,sDay)
endPlot = startPlot + dt.timedelta(365)
plotDates = [startPlot, endPlot]

fig, ax1 = plt.subplots(figsize=(15,10)) 
ax1.plot(xDays,y_obs, label='Observed')
ax1.plot(xDays,y_mod, label='Test Model output')
#ax1.set_xlim(plotDates)
ax1.set_ylabel('Soil moisture')
ax1.set_xlabel('Days')
plt.legend(loc='upper left')
ax2 = ax1.twinx()
ax2.plot(xDays,p, label='precipitation', linewidth=0.1, color='grey')
ax2.set_ylabel('precip')
#ax2.set_xlim(plotDates)
plt.legend(loc='upper right')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

