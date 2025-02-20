#!/bin/sh
#
# Author: Sohail Rathore
# Date  : 11/02/2025
#

mpirun -np 1 ./EMProb_JVVb -m mesh/OxNanoSys0.mesh -r 2 -o 2 &> log.dat
