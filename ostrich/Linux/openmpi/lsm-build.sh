#!/bin/bash

cd ..

module load intel

for ver in 1.6.5 1.7.3 1.8.4 2.0.2; do
  if [ ! -d ./openmpi/$ver ]; then
    mkdir ./openmpi/$ver
  fi

  module load openmpi/gcc-4.8.x/$ver

  echo "Building serial application"
  make clean
  make OMPI_SER 2>&1 | tee make-serial.log
  cp Ostrich ./openmpi/$ver

  echo "Building parallel application"
  make clean
  make OMPI 2>&1 | tee make-mpi.log
  cp OstrichMPI ./openmpi/$ver

  echo "Building IsoFit application"
  make clean
  make OMPI_ISO 2>&1 | tee make-isofit.log
  cp IsoFit ./openmpi/$ver

  module unload openmpi/gcc-4.8.x/$ver
done

