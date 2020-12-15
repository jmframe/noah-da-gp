#!/bin/bash
args=("$@")
echo ${args[0]} 
cd runs/${args[0]}
. /usr/share/modules/init/bash
#module purge
module load comp/intel/20.0.0.166
module load mpi/impi/20.0.0.166
ulimit -s unlimited
# UNCOMMENT THE FOLLOWING LATER?
#sh periodic_cleanup.sh
./noah_mp.exe >& noah_mp_exe.out_New2_OL
cd ../..

