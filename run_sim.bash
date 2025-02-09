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

module load openmpi-gcc7/4.0.4-cuda11.2
export OMP_NUM_THREADS=32
export OMP_PLACES=cores
export OMP_PROC_BIND=close

ulimit -s 10240
mpirun -np 2 ./mainEMSim &> Sim.log