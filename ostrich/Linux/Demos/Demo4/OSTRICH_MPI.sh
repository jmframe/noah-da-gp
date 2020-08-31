#!/bin/bash

# use parallel version of input file
cp ostIn_Parallel.txt ostIn.txt

# initialize counter used by SaveBest.sh
echo 0 > Counter.txt

# load modules or set path so that it includes desired mpi launcher
module load intel
module load intel-mpi

# match assignment to location of OSTRICH installation
OSTRICH_MPI=/projects/ccrstaff/lsmatott/research/Ostrich/OstrichMPI

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
srun -n 5 $OSTRICH_MPI

# pause for user input
read -n1 -r -p "Press any key to continue..." key

