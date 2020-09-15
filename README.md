# noah-da-gp
### Multivariate Hydrologic Data Assimilation for Model Structural Learning and Process-Diagnostics
#### This github repository will be specifically for using the PLUMBER-2 data,
#### but let's try to keep it as general as possible, so we can use other datasets.
#### In particular we'll want to assimilate satellite data.
#### 0. Calibrate Noah-MP for all the sites
* Let's hold off on Calibration for now.
* We need to decide how to break up the calibration, if at all.
    * If we calibrate with a 'leave one year out' strategy we will run up against data limitations, as some PLUMBER-2 sites only have one year.
#### 1. Running data assimilation on Noah-MP
* Soni is currently working on the multivariate data assimilation code. Jonathan will start contributing to it once it is on this github repository
* We'll do data assimilation seperately on each site. 
* QUESTION: Should we assimilate data throughout the entire record?
* We should decide how to do this consistently across sites:
* Spinup. We will equilibarate the model runs looping through the whole record X number of times, with each time setting the states at the beginning equal to the states at the previous cycles end.
#### 2. Train a Gaussian Process on the states of the assimilated Noah-MP model
* We will train the Gaussian Process on most sites, leaving out a random set of X validation sites, and one test site. We will cycle through this testing each basin once?
#### 3. Predict model structural error in Noah-MP with the GP to improve model performance
* QUESTION: Do we test the who record at each site, or some number of years? Sites have different years, so I (jmframe) suggest we test the whole record, reguardless of how many years are in it

## WARNINGS:
* Check the .gitignore file, I copied it in from another project, and it might prevent you from pushing something you want.
* Let's be very careful about the number of files generated and saved during runs. My (jmframe) file quota is pretty low.

## Genetal TODO list
* CHECK WHY LOCAL PARAMETERS NOT WORKING
* TROUBLESHOOT THE POOR PERFORMANCE OF QE, QH, NEE. MIGHT BE DUE TO LOCAL PARAMETERS... BAD INPUTS... UNIT CONVERSION...
* Keep this README up do date. 
* Need to extract all the local information for the noah runs, the initial plant/soil states, the general parameters and time offset.
* Finish the multi-objective data assimilation code. Soni has a first version of this. Will try running once the run directories are working.
* Copy in Craig's GP code. The MPI version runs in parallel and is for training. The non-MPI version is for making predictions. 
* Make a run directory, and the ability to fill it with a directory for each site:
    * For the data assimilation runs
        * Add a link to the raw data for processing. Process data as part of the run. Do not save the input data.
        * After all the data assimilation runs have been complete, save all the states in a pickle file, with the site IDs as keys
    * For the GPR training runs
        * Add a link to the raw data for processing. Process data before the run. 
        * Save the Noah-MP input data, since it will be used for training the GP multiple times. But delete after run.
        * Save the GP training data as a pickle file. This will be one file for each site. The file will have dynamic the data from every other site in one array.
        * Save the GP states after training. Copy it into the testing directory.
    * For the run testing the GP dynamic state update
        * Add a link to the raw data for processing. Process data as part of the run. Do not save the input data. 
        * Save the simulated values as output.
* Add a data directory. The run directories will link to this data directory, so we won't need to copy data. Let's avoid copying raw data as much as possible. 

## Completed tasks
* 08-31-2020 [jmframe]: Initialized github repository
* 09-05-2020 [Soni]: Finished test version of the multiparameter EnKS function
* 09-08-2020 [jmframe]: Got directory setup roughly working, first successful test run.
* 09-09-2020 [jmframe]: Extracting local parameters for each site. 
* 09-11-2020 [jmframe]: Setup objective function to calculate for qe, qh & nee
