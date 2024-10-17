# MFEM_OxfordNanoSystems
A simple electromagnetics case expanded from MFEM
tutorial ex5p solving the Darcy problem using the
same finite element spaces etc.. but different physical
interpertations
```math
\displaylines{k \vec{J} + \nabla v = \vec{J_e}, \\
- \nabla \cdot \vec{J} = Q,}
```
The system solves these equations in using mixed H(div),
H1 spaces, using a mixture of Raviart-Thomas and Lagrange
elements.

This Library is dependant on MFEM and HYPRE, MPI and a C-
however these are straight forward to install.