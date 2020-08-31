#!/bin/bash

#for i in Demo*; do
for n in 9; do
  i=Demo$n
  echo "Running $i"
  cd $i

    mkdir -p test/serial
    cp * test/serial
    cd test/serial
      ./OSTRICH.sh 
    cd ../..

    mkdir -p test/parallel
    cp * test/parallel
    cd test/parallel
      ./OSTRICH_MPI.sh
    cd ../..
  cd ..
done

