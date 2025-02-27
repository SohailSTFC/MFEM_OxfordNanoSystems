//                       MFEM Electromagnetic Sample - Parallel Version
//
// Compile with: make
//
// Sample runs:  mpirun -np 1 EMSampleProb -m mesh/OxNanoSys0.mesh
//               mpirun -np 1 EMSampleProb -m mesh/OxNanoSys1.msh
//
//
// Description:  This example code solves a simple 2D/3D mixed Darcy problem
//               corresponding to the saddle point system
//
//                                 J + grad v = f
//                                 - div J    = g
//
//               with natural boundary condition -v = <given potenialt>.
//               Here, we use a given exact solution (J,v) and compute the
//               corresponding r.h.s. (f,g).  We discretize with Raviart-Thomas
//               finite elements (Currect flux J) and piecewise discontinuous
//               polynomials (potential v).
//
//               The example demonstrates the use of the BlockOperator class, as
//               well as the collective saving of several grid functions in the
//               ParaView (paraview.org) format.
//
//

#include "mfem.hpp"
#include <fstream>
#include <iostream>

// Define the analytical solution and forcing terms / boundary conditions
#include "include/BoundaryAndInitialSolution.hpp"
#include "include/DarcyEMProblem.hpp"
#include "include/Visualisation.hpp"

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();
   bool verbose = (myid == 0);

   //Boundary condition arrays
   Array<int> dbcs;
   Array<int> nbcs;
   Vector dbcv;

   //Mat props
   double sig = 5.500E-06, MU = 1.257E-06;

   // 2. Parse command-line options.
   const char *mesh_file = "mesh/OxNanoSysU0.msh";
   int ref_levels = -1;
   int order = 2;
   bool pa = false;
   const char *device_config = "cpu"; //"cpu";//"ceed-cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&sig, "-p", "--perm",
                  "The permiability of the electrolyte (Sigma).");


   args.Parse();
   if (!args.Good())
   {
      if (verbose) args.PrintUsage(cout);
      return 1;
   }
   if (verbose) args.PrintOptions(cout);

   // 3. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   // 4. Read the (serial) mesh from the given mesh file on all processors.  We
   //    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   //    and volume meshes with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();



   // 5. Refine the serial mesh on all processors to increase the resolution. In
   //    this example we do 'ref_levels' of uniform refinement. We choose
   //    'ref_levels' to be the largest number that gives a final mesh with no
   //    more than 10,000 elements, unless the user specifies it as input.
   {
      if (ref_levels == -1) ref_levels = (int)floor(log(10000./mesh->GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++) mesh->UniformRefinement();
   }

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   {
      int par_ref_levels = 1;
      for (int l = 0; l < par_ref_levels; l++) pmesh->UniformRefinement();
   }

   // 7. Define a parallel finite element space on the parallel mesh. Here we
   //    use the Raviart-Thomas finite elements of the specified order.
   FiniteElementCollection *RT_coll(new RT_FECollection(order+1, dim));
   FiniteElementCollection *H1_coll(new H1_FECollection(order));

   ParFiniteElementSpace *R_space = new ParFiniteElementSpace(pmesh, RT_coll);
   ParFiniteElementSpace *W_space = new ParFiniteElementSpace(pmesh, H1_coll);

   //Set up the problem
   MemoryType mt = device.GetMemoryType();
   DarcyEMProblem demoProb(R_space, W_space, sig, mt, dim);
    
   //Set the solver and preconditioner
   demoProb.BuildPreconditioner();
   demoProb.Set_Solver(verbose);

   //Solve the equations
   demoProb.Solve(verbose);

   //Visualise the results using ParaView
   double time = 0.0;
   demoProb.SetFields();
   ParaViewVisualise("EMSampleProbJV",demoProb.Fields, demoProb.FieldNames, order, pmesh, time);

   // 20. Free the used memory.
   delete W_space;
   delete R_space;
   delete RT_coll, H1_coll;
   delete pmesh;
   return 0;
}
