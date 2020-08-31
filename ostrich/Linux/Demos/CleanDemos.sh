#!/bin/bash

for i in Demo*; do
  echo "Cleaning $i"
  cd $i
    rm -Rf test
  cd ..
done

