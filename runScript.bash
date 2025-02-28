#!/bin/sh
#BSUB -q scafellpikeSKL
#BSUB -o stdout.%J.log
#BSUB -e stderr.%J.err
#BSUB -R "span[ptile=16]"
#BSUB -R "rusage[mem=10000]"
#BSUB -n 32
#BSUB -J MyMPI-job
#BSUB -W 02:00

source /lustre/scafellpike/local/HT06187/uxq08/shared/envLoader.sh
ulimit -s 10240
mpirun -np 32 EMsampleProbV   -m mesh/OxNanoSys1.msh -r 0 -p 2 &> Sim_V_USH0.log

mpirun -np 32 EMsampleProbJV  -m mesh/OxNanoSys1.msh -r 0 -p 2 &> Sim_JV_USH0.log

mpirun -np 32 EMsampleProbJVB -m mesh/OxNanoSys1.msh -r 0 -p 2 &> Sim_JVB_USH0.log

