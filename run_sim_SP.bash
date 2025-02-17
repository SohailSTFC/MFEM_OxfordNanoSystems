#!/bin/sh
#BSUB -q scafellpikeSKL
#BSUB -o stdout.%J.log
#BSUB -e stderr.%J.err
#BSUB -R "span[ptile=32]"
#BSUB -R "rusage[mem=10000]"
#BSUB -x
#BSUB -n 64
#BSUB -J MyMPI-job
#BSUB -W 00:50


#
# Runs a sample simulation
# on scafell pike, the format
# is suited for the LSF Job
# scheduler
#
# Author: Sohail Rathore
# Date  : 31/01/2025


ulimit -s 10240
mpirun -np 1 ./EMProb_JVVb -m "mesh/OxNanoSys0.mesh"  \ #Input mesh
                           -r  3                      \ #Refinement level
                           -o  2                      \ #Polynomial Order
                           -nbcs  \ #Neumann BC boundary surfaces
                           -dbcs  \ #Dirchlet BC boundary surfaces
                           -dbcv    #Dirchlet BC boundary values

