#!/gpfsm/dulocal/sles11/other/SLES11.3/miniconda3/2019.03_py3.7/2019-05-15/bin/python
import numpy as np
import shutil, errno
import os
import sys
from io import StringIO
site = sys.argv[1]
stepfile = sys.argv[2]

maindir = '/discover/nobackup/jframe/gpr_fluxnet_oss/'
gpr_mpi_dir = maindir+'soil_moisture/gpr-mpi/'+str(site)+'/'

# Define the number of parameiters that will be used for training the GPR
numGPRinputs_wet = 25
numGPRinputs_dry = 24

# Include a lag to stabalize the GPR results
# Until otherwise specified this lag is three timesteps.
lag = 3 

# A list used for generating the GPR target
dif = []

wet_train_input =np.zeros((1,25))
dry_train_input =np.zeros((1,25))
wet_train_target=np.zeros((25))
dry_train_target=np.zeros((25))
print('Shape of the wet training input:', wet_train_input.shape)
print('Shape of the dry training input:', dry_train_input.shape)
print('Shape of the wet training target:', wet_train_target.shape)
print('Shape of the dry training target:', dry_train_target.shape)

sitelist = [1,2,5,7,11,13,14,18]

for iSite in sitelist:
    print('Working on site '+str(iSite))
    if iSite == int(site):
        print('Skipping these data, because this is the testing site')
        continue

    gpr_noah_dir = maindir+'soil_moisture/gpr-noah/train_'+str(iSite)+'/'
    # List to cut down the training data into just what we want to use.
    wet_train_list = []
    dry_train_list = []

    # Grabbing data from two files, and writing to one. Open them at the same time.
    # Note that soil moisture is column 5 (sixth column) of observations and 
      # column 10 (11th column) of ouput.out
    with open(gpr_noah_dir + 'obs.txt', 'r') as Observations_file, \
         open(gpr_noah_dir + 'forcing.txt', 'r') as Forcing_file, \
         open(gpr_noah_dir + stepfile, 'r') as Onestep_file:
    
        # Read the input data as delimited.
        forcing_data = np.genfromtxt(Forcing_file)
        observation_data = np.genfromtxt(Observations_file)
        onestep_data = np.genfromtxt(Onestep_file)
    
        lonlat = np.genfromtxt(gpr_noah_dir + 'lat_lon.txt')
        latitude = lonlat[0]

        # Get the overall number of records (not just useful observations)
        # And loop through them looking for good points to train.
        with open(gpr_noah_dir + 'num_times.txt') as nt:
            total_records = int(nt.readline())
        print('Number of simulation timesteps = '+str(total_records)) 
        # These are the indexes of usable observations 
        # Not only do they have good values, but so do their lags.
        iObs = []
        for i in range(lag, total_records):
            if all(q > 0 for q in observation_data[i-3:i+1, 5]):
                # Add in some criteria to evaluate the truth in the observation
                # one step change (osc)
                osc = observation_data[i - 1] - observation_data[i]
                if abs(osc) < 0.1:
                    iObs.append(i)
    
        all_target = np.zeros(len(iObs)) 
        all_input = np.zeros((len(iObs), numGPRinputs_wet)) 
        #########################################################################################
        #########   MAIN LOOP   ############   MAIN LOOP   ###############   MAIN LOOP  #########
        #########################################################################################
        # Loop through the usable observations. Pull the training data from values associated
        # with these locations. But then place them in the freshly generated 'all_input/target'
        # arrays starting from zero (intarloc)
        intarloc=0 #Input and target array locations
        for i in iObs:
    
            # Make sure all the data are good.
            if all(q > 0 for q in observation_data[i-3:i+1, 5]) is False:
                print('Warning: hit a bad observation in', siteyear)
                print([i-3, i-2, i-1, i])
                print(observation_data[i-3:i+1,5])
    
            #####################################################################################
            #############     GPR TARGET DATA                ####################################
            #####################################################################################
            # Take the difference between observation and 1-step result, 
            # this is what the GPR trains on as 'Target'
            all_target[intarloc] = float(observation_data[i][5]) - float(onestep_data[i][10])
            
            #####################################################################################
            #############     GPR PREDICTION INPUT DATA      ####################################
            #####################################################################################
            # winputs(1) = lagged(3)
            all_input[intarloc][0] = float(observation_data[i-1][5])
           # winputs(2) = lagged(2)
            all_input[intarloc][1] = float(observation_data[i-2][5])
            #  winputs(3) = lagged(1)
            all_input[intarloc][2] = float(observation_data[i-3][5])
            #  winputs(4) = state(t)%smc(1)
            all_input[intarloc][3] = float(onestep_data[i][10])
            #  winputs(5) = forcing(t)%q2
            all_input[intarloc][4] = float(forcing_data[i][6])
            #  winputs(6) = forcing(t)%prcprate                    # DELETE THIS BELOW FOR DRY
            all_input[intarloc][5] = float(forcing_data[i][10])    # DELETE THIS BELOW FOR DRY
            #  winputs(7) = forcing(t)%lwrad
            all_input[intarloc][6] = float(forcing_data[i][9])
            #  winputs(8) = forcing(t)%swrad
            all_input[intarloc][7] = float(forcing_data[i][8])
            #  winputs(9) = forcing(t)%sfcprs
            all_input[intarloc][8] = float(forcing_data[i][7])
            #  winputs(10) = forcing(t)%sfctmp
            all_input[intarloc][9] = float(forcing_data[i][5])
            #  winputs(11)= forcing(t)%sfcspd
            all_input[intarloc][10] = float(forcing_data[i][3])
            #  winputs(12)= state(t)%smc(2)
            all_input[intarloc][11] = float(onestep_data[i][11])
            #  winputs(13)= state(t)%smc(3)
            all_input[intarloc][12] = float(onestep_data[i][12])
            #  winputs(14)= state(t)%smc(4)
            all_input[intarloc][13] = float(onestep_data[i][13])
            #  winputs(15)= state(t)%sh2o(1)
            all_input[intarloc][14] = float(onestep_data[i][14])
            #  winputs(16)= state(t)%sh2o(2)
            all_input[intarloc][15] = float(onestep_data[i][15])
            #  winputs(17)= state(t)%sh2o(3)
            all_input[intarloc][16] = float(onestep_data[i][16])
            #  winputs(18)= state(t)%sh2o(4)
            all_input[intarloc][17] = float(onestep_data[i][17])
            # winputs(19)= state(t)%lai
            all_input[intarloc][18] = float(onestep_data[i][22])
            # winputs(20)= state(t)%fastcp
            all_input[intarloc][19] = float(onestep_data[i][24])
            # winputs(21)= state(t)%stblcp
            all_input[intarloc][20] = float(onestep_data[i][25])
            # winputs(22)= output(t)%qe
            all_input[intarloc][21] = float(onestep_data[i][26])
            # winputs(23)= output(t)%qh
            all_input[intarloc][22] = float(onestep_data[i][27])
            # winputs(24)= output(t)%nee
            all_input[intarloc][23] = float(onestep_data[i][28])
            # winputs(25)= abs(setup%latitude)
            all_input[intarloc][24] = float(abs(latitude))
    
            # USE SOME LOGIC TO SPLIT THE TRAINING SET BETWEEN WET AND DRY
            if float(forcing_data[i][10]) > 0:
                wet_train_list.append(intarloc)
            else:
                dry_train_list.append(intarloc)
            intarloc+=1
    
    #########################################################################################
    ####### END MAIN LOOP   ########## END MAIN LOOP   ############# END MAIN LOOP  #########
    #########################################################################################

    print('total training samples:', len(iObs))
    print('Shape of all the input data:', all_input.shape)
    # Split up the data for wet/dry and training/prediction
    wti = all_input[wet_train_list]
    dti = all_input[dry_train_list]
    wtt = all_target[wet_train_list]
    dtt = all_target[dry_train_list]
    print('wet training data length:', len(wet_train_list))
    print('dry training data length:', len(dry_train_list))
    print('sample training rows:')
    print('wet')
    print(wti[0,:])
    print('dry')
    print(dti[0,:])
    print('Shape of the wet training input:', wti.shape)
    print('Shape of the dry training input:', dti.shape)
    print('Shape of the wet training target:', wtt.shape)
    print('Shape of the dry training target:', dtt.shape)
    # Append the site vector to the total training vector.
    wet_train_input  = np.append(wet_train_input,wti, axis=0)
    dry_train_input  = np.append(dry_train_input, dti, axis=0)
    wet_train_target = np.append(wet_train_target, wtt, axis=0)
    dry_train_target = np.append(dry_train_target, dtt, axis=0)
    
    print('Shape of the wet training input:', wet_train_input.shape)
    print('Shape of the dry training input:', dry_train_input.shape)
    print('Shape of the wet training target:', wet_train_target.shape)
    print('Shape of the dry training target:', dry_train_target.shape)
wet_train_input  = np.delete(wet_train_input, 0, 0) 
dry_train_input  = np.delete(dry_train_input, 0, 0)
wet_train_target = np.delete(wet_train_target, 0, 0)
dry_train_target = np.delete(dry_train_target, 0, 0)
print('Shape of the wet training input:', wet_train_input.shape)
print('Shape of the dry training input:', dry_train_input.shape)
print('Shape of the wet training target:', wet_train_target.shape)
print('Shape of the dry training target:', dry_train_target.shape)
# Delete out the precipitation column.
dry_train_input = np.delete(dry_train_input, 5, 1)
print('Shape of the wet training input:', wet_train_input.shape)
print('Shape of the dry training input:', dry_train_input.shape)
print('Shape of the wet training target:', wet_train_target.shape)
print('Shape of the dry training target:', dry_train_target.shape)

#print('Row one before shuffle:')
#print('dry')
#print(dry_train_input[0,:])
#print('wet')
#print(wet_train_input[0,:])
#np.random.shuffle(dry_train_input)
#np.random.shuffle(wet_train_input)
#print('Row one after shuffle:')
#print('dry')
#print(dry_train_input[0,:])
#print('wet')
#print(wet_train_input[0,:])

# Write targets to binary file
with open(gpr_mpi_dir+'targets_wet.bin','wb') as target_file:
    np.array(wet_train_target, dtype=np.float64).tofile(target_file)
with open(gpr_mpi_dir+'targets_dry.bin','wb') as target_file:
    np.array(dry_train_target, dtype=np.float64).tofile(target_file)

# Write training inputs to binary file
with open(gpr_mpi_dir+'inputs_wet.bin','wb') as input_file:
    np.array(wet_train_input, dtype=np.float64).tofile(input_file)
# Write training inputs to binary file
with open(gpr_mpi_dir+'inputs_dry.bin','wb') as input_file:
    np.array(dry_train_input, dtype=np.float64).tofile(input_file)

## Write out the number of values for the training data. 
with open(gpr_mpi_dir+'num_wet_training_points.txt', 'w+') as n_file:
    n_file.write(str(wet_train_input.shape[0]))
with open(gpr_mpi_dir+'num_dry_training_points.txt', 'w+') as n_file:
    n_file.write(str(dry_train_input.shape[0]))
##################    END PROGRAM    ################################
