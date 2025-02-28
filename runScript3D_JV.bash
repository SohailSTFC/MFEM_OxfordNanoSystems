#!/bin/sh
#BSUB -q scafellpikeSKL
#BSUB -o stdout.%J.log
#BSUB -e stderr.%J.err
#BSUB -R "span[ptile=128]"
#BSUB -R "rusage[mem=10000]"
#BSUB -n 128
#BSUB -J MyMPI-job
#BSUB -W 02:00

source /lustre/scafellpike/local/HT06187/uxq08/shared/envLoader.sh
ulimit -s 10240


##
## Unshielded Case
##
mpirun -np 128 EMsampleProbJV -m mesh/OxNanoSys3.msh -r 0 -p 2 &> JobData/Sim_3D_JV_USH1.log
cd ParaView
mkdir Sim_3D_USH1
mv EMsampleProbJV   Sim_3D_USH1/EMsampleProbJV
cd ..

##
## Shielded Case
##
mpirun -np 128 EMsampleProbJV -m mesh/OxNanoSys4.msh -r 0 -p 2 &> JobData/Sim_3D_JV_SH1.log
cd ParaView
mkdir Sim_3D_SH1
mv EMsampleProbJV  Sim_3D_SH1/EMsampleProbJV
cd ..
