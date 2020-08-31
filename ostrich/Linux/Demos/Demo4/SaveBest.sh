#!/bin/bash

# set a fixed path for the "best" folder
MYROOT=`pwd`

# read in the counter
A=`cat $MYROOT/Counter.txt`

# increment the counter
A=`expr $A + 1`

# create the best directory
mkdir $MYROOT/best 2>/dev/null

# Create a subdir for the latest best solution.
# This is optional, could just overwrite previous best
# if only interest in saving files associated with the
# ulitmate best objective function.
mkdir $MYROOT/best/$A 

# copy model files to new subdir
cp * $MYROOT/best/$A 

# save counter to file for next time
echo $A > $MYROOT/Counter.txt

