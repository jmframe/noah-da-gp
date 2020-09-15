#!/bin/bash

export LIB_NETCDF="/discover/nobackup/projects/lis/libs/netcdf/4.7.4_intel-20.0.0.166/lib"
export LIB_HDF5="/discover/nobackup/projects/lis/libs/hdf5/1.10.6_intel-20.0.0.166/lib"
export LD_LIBRARY_PATH=$(LIB_NETCDF):$(LIB_HDF5)
source export ../../../environment_setup_june2020.sh
