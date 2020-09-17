#!/bin/bash

source ../../../environment_setup_june2020.sh
LIB_NETCDF="/discover/nobackup/projects/lis/libs/netcdf/4.7.4_intel-20.0.0.166/lib"
LIB_HDF5="/discover/nobackup/projects/lis/libs/hdf5/1.10.6_intel-20.0.0.166/lib"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIB_NETCDF:$LIB_HDF5
echo $LD_LIBRARY_PATH
