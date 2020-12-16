#!/discover/nobackup/jframe/anaconda3/bin/python

# This script is meant to mathc the cumulative distributions
# of the observed soil moisture with the ensemble mean soil moisture. 
# This should be done before data assimilation.
# Probably want to include ensemble generation and cfd Matching in the initialization.
# Adapted from Grey Neaing's code, July 2019.

import csv
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy.random
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# --- CDF-Matching Function ---------------------------------
# Not sure about keeping this as a function, because probably want 
# to run it independently with input from text and output to text.
###################################################################
####   def cdf_matching(plot_indx,A,B,Apred,nQuantiles=1e2):  #####
###################################################################

# Parameter for the cumulative distribution function.
nQuantiles=1e3

# Observed soil moisture, will be changed to match CDF of model ensemble
A = []
# model Ensemble soil moisture.
B = []
count_out = 0
count_obs = 0

with open('open_mean.out','r') as csvfile:
    EnksMean = csv.reader(csvfile, delimiter=' ', skipinitialspace=True )
    for row in EnksMean:
        B.append(float(row[7]))
        count_out = count_out + 1

with open('obs.orig','r') as csvfile:
    observations = csv.reader(csvfile, delimiter=' ', skipinitialspace=True )
    for row in observations:
        A.append(float(row[5]))
        count_obs = count_obs + 1
        if count_obs == count_out:
            break

A = np.asarray(A, dtype = float)
Apred = A
B = np.asarray(B, dtype = float)

# empirical pdfs
A_pdf, A_edges = np.histogram(A, bins=int(nQuantiles))
A_centers = A_edges[1:] - ((A_edges[1] - A_edges[0])/2)
A_cdf = np.cumsum(A_pdf) / len(A)

B_pdf, B_edges = np.histogram(B, bins=int(nQuantiles))
B_centers = B_edges[1:] - ((B_edges[1] - B_edges[0])/2)
B_cdf = np.cumsum(B_pdf) / len(B)

# cubic spline interpolation of the cumulative distribution function
A_centers = A_centers[A_pdf>0] 
A_cdf = A_cdf[A_pdf>0] 
itp_A_cdf = InterpolatedUnivariateSpline(A_centers, A_cdf, k=3)

# the PPF is the inverse of the CDF, so we simply reverse the order of the
# x & y arguments to InterpolatedUnivariateSpline
B_centers = B_centers[B_pdf>0] 
B_cdf = B_cdf[B_pdf>0] 
itp_B_ppf = InterpolatedUnivariateSpline(B_cdf, B_centers, k=3)

# map the original data 
Apred_corrected = itp_B_ppf(np.clip(itp_A_cdf(Apred),0,1))

# stuff for plotting
A_corrected = itp_B_ppf(np.clip(itp_A_cdf(A),0,1))

C_pdf, C_edges = np.histogram(A_corrected, bins=int(nQuantiles))
C_centers = C_edges[1:] - ((C_edges[1] - C_edges[0])/2)
C_cdf = np.cumsum(C_pdf) / len(A)

P_pdf, P_edges = np.histogram(Apred_corrected[~np.isnan(Apred_corrected)], bins=int(nQuantiles))
P_centers = P_edges[1:] - ((P_edges[1] - P_edges[0])/2)
P_cdf = np.cumsum(P_pdf) / len(Apred[~np.isnan(Apred_corrected)])

# plot stuff
#fig, ax = plt.subplots(1, 1)
#ax.plot(A_centers, A_cdf, '-r', lw=3, label='A CDF')
#ax.plot(B_centers, B_cdf, '-k', lw=3, label='B CDF')
#ax.plot(C_centers, C_cdf, '-b', lw=3, label='A Matched to B')
#ax.plot(P_centers, P_cdf, '--g', lw=1, label='A* Matched to B')
#ax.legend(loc=5)
#plt.show()

with open('obs.orig','r') as readObsOrig, open('obs.txt', 'w+') as writeMatchResults:
    observations = csv.reader(readObsOrig, delimiter=' ', skipinitialspace=True )
    for obs_row, iA in zip(observations, Apred_corrected):
        obs_row[5] = iA
        for item in obs_row:
            writeMatchResults.write('%s ' % item)
        writeMatchResults.write('\n')

print('CFDs matched')
##################################################################
###    return Apred_corrected   ##################################
##################################################################
