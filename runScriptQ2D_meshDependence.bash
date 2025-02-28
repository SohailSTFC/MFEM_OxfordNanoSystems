#!/bin/sh
#BSUB -q scafellpikeSKL
#BSUB -o stdout.%J.log
#BSUB -e stderr.%J.err
#BSUB -R "span[ptile=32]"
#BSUB -R "rusage[mem=10000]"
#BSUB -n 32
#BSUB -J MyMPI-job
#BSUB -W 02:00

source /lustre/scafellpike/local/HT06187/uxq08/shared/envLoader.sh
ulimit -s 10240

##
## Agressive routine deletes 
## all data after each use
## comment out lines 20-27
## to remove this behavior
##
rm -rf JobData ParaView
mkdir JobData
mkdir ParaView

##
## Unshielded Cases
##
mpirun -np 32 EMsampleProbV   -m mesh/OxNanoSys1.msh -r 0 -p 2 &> JobData/Sim_V_USH0.log
mpirun -np 32 EMsampleProbJV  -m mesh/OxNanoSys1.msh -r 0 -p 2 &> JobData/Sim_JV_USH0.log
mpirun -np 32 EMsampleProbJVB -m mesh/OxNanoSys1.msh -r 0 -p 2 &> JobData/Sim_JVB_USH0.log
cd ParaView
mkdir Sim_Q2_USH0
mv EMsampleProbV   Sim_Q2_USH0/EMsampleProbV
mv EMsampleProbJV  Sim_Q2_USH0/EMsampleProbJV
mv EMSampleProbJBv Sim_Q2_USH0/EMSampleProbJBv
cd ..

mpirun -np 32 EMsampleProbV   -m mesh/OxNanoSys1.msh -r 1 -p 2 &> JobData/Sim_V_USH1.log
mpirun -np 32 EMsampleProbJV  -m mesh/OxNanoSys1.msh -r 1 -p 2 &> JobData/Sim_JV_USH1.log
mpirun -np 32 EMsampleProbJVB -m mesh/OxNanoSys1.msh -r 1 -p 2 &> JobData/Sim_JVB_USH1.log
cd ParaView
mkdir Sim_Q2_USH1
mv EMsampleProbV   Sim_Q2_USH1/EMsampleProbV
mv EMsampleProbJV  Sim_Q2_USH1/EMsampleProbJV
mv EMSampleProbJBv Sim_Q2_USH1/EMSampleProbJBv
cd ..


##
## Shielded Cases
##
mpirun -np 32 EMsampleProbV   -m mesh/OxNanoSys1s.msh -r 0 -p 2 &> JobData/Sim_V_SH0.log
mpirun -np 32 EMsampleProbJV  -m mesh/OxNanoSys1s.msh -r 0 -p 2 &> JobData/Sim_JV_SH0.log
mpirun -np 32 EMsampleProbJVB -m mesh/OxNanoSys1s.msh -r 0 -p 2 &> JobData/Sim_JVB_SH0.log
cd ParaView
mkdir Sim_Q2_SH0
mv EMsampleProbV   Sim_Q2_SH0/EMsampleProbV
mv EMsampleProbJV  Sim_Q2_SH0/EMsampleProbJV
mv EMSampleProbJBv Sim_Q2_SH0/EMSampleProbJBv
cd ..

mpirun -np 32 EMsampleProbV   -m mesh/OxNanoSys1s.msh -r 1 -p 2 &> JobData/Sim_V_SH1.log
mpirun -np 32 EMsampleProbJV  -m mesh/OxNanoSys1s.msh -r 1 -p 2 &> JobData/Sim_JV_SH1.log
mpirun -np 32 EMsampleProbJVB -m mesh/OxNanoSys1s.msh -r 1 -p 2 &> JobData/Sim_JVB_SH1.log
cd ParaView
mkdir Sim_Q2_SH1
mv EMsampleProbV  Sim_Q2_SH1/EMsampleProbV
mv EMsampleProbJV Sim_Q2_SH1/EMsampleProbJV
mv EMSampleProbJBv Sim_Q2_SH1/EMSampleProbJBv
cd ..