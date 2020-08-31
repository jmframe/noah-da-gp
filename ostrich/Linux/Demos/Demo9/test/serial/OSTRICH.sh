#!/bin/bash

# use serial version of input file
cp OstIn_serial.txt ostIn.txt

# set path to OSTRICH binary
OSTRICH=/projects/ccrstaff/lsmatott/research/Ostrich/Ostrich

$OSTRICH

# pause for user input
read -n1 -r -p "Press any key to continue..." key

