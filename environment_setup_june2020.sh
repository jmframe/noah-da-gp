#!/bin/bash -f

# modules
module purge
module load comp/intel/19.1.0.166
module load mpi/impi/20.0.0.166

# module load lib/mkl-18.0.0.128 
# module load other/cmake-3.8.2 
# module load other/comp/gcc-4.8.1
# 
# # Libraries
# export FORTRAN_HOME=/usr/local/intel/Composer/composer_xe_2013_sp1.3.174/compiler/lib/intel64
# export BLAS_HOME=/usr/local/other/SLES11/BLAS/intel-14.0.3.174
# export LAPACK_HOME=/usr/local/other/SLES11/lapack/3.5.0/intel-14.0.3.174/lib
# export SCALAPACK_HOME=/usr/local/other/SLES11/ScaLAPACK/2.0.2/intel-14.0.3.174/lib
# export GSL_HOME=/usr/local/other/SLES11/gsl/2.4/intel-14.0.3.174
# export CPPUNIT_HOME=/discover/nobackup/cpelissi/libs/intel/cppunit
# export GXX_NAME=/usr/local/other/SLES11.1/gcc/4.8.1/bin/gcc
# 
# # linking
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/discover/nobackup/projects/lis/libs/netcdf/4.3.3.1_intel-14.0.3.174_sp3/lib
# #export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/discover/nobackup/projects/lis/libs/grib_api/1.15.0_intel-14.0.3.174_sp3/lib
# 

