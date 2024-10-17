# MFEM_OxfordNanoSystems
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
C-compiler (GCC was used) however these are straight forward to install.