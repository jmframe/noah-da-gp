#!/bin/bash

# set path to OSTRICH binary
OSTRICH=/projects/ccrstaff/lsmatott/research/Ostrich/Ostrich

# use the serial version of the input file
cp ostIn_Serial.txt ostIn.txt

# initialize counter used by SaveBest.sh
echo 0 > Counter.txt

$OSTRICH

# pause for user input
read -n1 -r -p "Press any key to continue..." key

