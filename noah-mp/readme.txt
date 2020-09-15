My (Jonathan Frame's) master copy of the noah-mp code.


to do:


version    Note
------------------------------------------------------------------
20191204   This is the code version used for the results presented
at AGU 2019: Theory Guided Machine Learning to Improve Hydrology Models.

20190717   Implemented a one-step state upade option. To run the
One-step state update set the da_flag.txt to -1.
Added some logic so when the model is run 
withouot Data Assimilation it does not try to read in DA specific 
files for sig_sm, Nlag and obs_cov.txt.

20190716   I added some code to run the Ensemble Kalman Smoother, 
This was mostly implimented by Grey, but I made a few changes, 
like the perturbations to the forcing data.

