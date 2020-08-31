#!/bin/bash

cd ..

#module load intel

for ver in 4.1.3 5.0.2 5.1.1 2017.0.1; do
  if [ ! -d ./intelmpi/$ver ]; then
    mkdir ./intelmpi/$ver
  fi

#  module load intel-mpi/$ver

  echo "Building serial application"
  make clean
  make SER 2>&1 | tee make-serial.log
  cp Ostrich ./intelmpi/$ver

  echo "Building parallel application"
  make clean
  make MPI 2>&1 | tee make-mpi.log
  cp OstrichMPI ./intelmpi/$ver

  echo "Building IsoFit application"
  make clean
  make ISO 2>&1 | tee make-isofit.log
  cp IsoFit ./intelmpi/$ver

  echo "Building GNU serial application"
  make clean
  make GCC 2>&1 | tee make-serial.log
  cp OstrichGCC ./intelmpi/$ver

#  module unload intel-mpi/$ver
done

