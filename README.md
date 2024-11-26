# MFEM Oxford Nano Systems
A simple electromagnetics case expanded from MFEM
tutorial ex5p solving the Darcy problem using the
same finite element spaces etc.. but different physical
interpertations. The system solves these equations in using mixed H(div),
H1 spaces.
```math
\displaylines{k \vec{J} + \nabla v = \vec{J_e}, \\
- \nabla \cdot \vec{J} = Q,}
```

This application is dependant on MFEM, METIS/PARMETIS, HYPRE, MPI (OpenMPI was used) and a 
C-compiler (GCC was used) however these are straight forward to install. The process of 
installing/getting this library is as follows:

1) Install MFEM using your [preferred installation method](https://mfem.org/building/)

2) Clone this repository alongside MFEM
  `git clone https://github.com/SohailSTFC/MFEM_OxfordNanoSystems.git`

3) Build the application with the following command,
  `cd MFEM_OxfordNanoSystems && make`

4) Run the application

# Sample problem 0:
A sample 2-D problem without any kind of shield around the plate to regularise the local field.
The electrodes are set to 3.0v and the plate is set to 0.0v, the outer boundary has a no-flux
boundary conditions ($\vec{n} \cdot \vec{J} = 0  $). This problem is the default problem and 
will automatically run if no other command line options are specified.
The initial mesh for this problem can be seen bellow:

![image](Imgs/Sample0Mesh.PNG)

This problem is run using the following commmand options if built with MPI 
(the number of processes can be increased to reflect the number of available 
cores however by default its set to 1):

`mpirun -np 1 ./EMsampleProb`.

or

`mpirun -np 1 ./EMsampleProb --mesh mesh/OxNanoSys0.mesh`

# Sample problem 1:
This uses MFEM features for external meshes based off of the [MFEM tutorial](https://mfem.org/tutorial/meshvis/) , it has the additional dependency
of using gmsh however this can be [cloned from gitlab](https://gitlab.onelab.info/gmsh/gmsh) and built for your system or downloaded as a [windows ready binary](https://gmsh.info/).
Gmsh can read in most mesh formats and output them into the Gmsh format which is readable from MFEM, however in this sample problem
the Native geometry format is used and a mesh is generated from that.
