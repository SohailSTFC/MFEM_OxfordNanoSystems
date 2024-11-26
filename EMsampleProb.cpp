//                       MFEM Example 5 - Parallel Version
//
// Compile with: make ex5p
//
// Sample runs:  mpirun -np 4 ex5p -m ../data/square-disc.mesh
//               mpirun -np 4 ex5p -m ../data/star.mesh
//               mpirun -np 4 ex5p -m ../data/beam-tet.mesh
//               mpirun -np 4 ex5p -m ../data/beam-hex.mesh
//               mpirun -np 4 ex5p -m ../data/escher.mesh
//               mpirun -np 4 ex5p -m ../data/fichera.mesh
//
//
// Description:  This example code solves a simple 2D/3D mixed Darcy problem
//               corresponding to the saddle point system
//
//                                 u + grad p = f
//                                 - div u      = g
//
//               with natural boundary condition -p = <given pressure>.
//               Here, we use a given exact solution (u,p) and compute the
//               corresponding r.h.s. (f,g).  We discretize with Raviart-Thomas
//               finite elements (velocity u) and piecewise discontinuous
//               polynomials (pressure p).
//
//               The example demonstrates the use of the BlockOperator class, as
//               well as the collective saving of several grid functions in
//               VisIt (visit.llnl.gov) and ParaView (paraview.org) formats.
//               Optional saving with ADIOS2 (adios2.readthedocs.io) streams is
//               also illustrated.
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

   // 2. Parse command-line options.
//   const char *mesh_file = "mesh/star.mesh";
   const char *mesh_file = "mesh/OxNanoSys0.mesh";
   int ref_levels = -1;
   int order = 2;
   bool par_format = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = 1;
   bool adios2 = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&par_format, "-pf", "--parallel-format", "-sf",
                  "--serial-format",
                  "Format to use when saving the results for VisIt.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&adios2, "-adios2", "--adios2-streams", "-no-adios2",
                  "--no-adios2-streams",
                  "Save data using adios2 streams.");
   args.Parse();
   if (!args.Good())
   {
      if (verbose)
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (verbose)
   {
      args.PrintOptions(cout);
   }

   // 3. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   // 4. Read the (serial) mesh from the given mesh file on all processors.  We
   //    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   //    and volume meshes with the same code.


/////void Mesh::ReadTrueGridMesh(std::istream &input)
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();



   // 5. Refine the serial mesh on all processors to increase the resolution. In
   //    this example we do 'ref_levels' of uniform refinement. We choose
   //    'ref_levels' to be the largest number that gives a final mesh with no
   //    more than 10,000 elements, unless the user specifies it as input.
   {
      if (ref_levels == -1)
      {
         ref_levels = (int)floor(log(10000./mesh->GetNE())/log(2.)/dim);
      }

      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
/*
   {
      int par_ref_levels = 2;
      for (int l = 0; l < par_ref_levels; l++)
      {
         pmesh->UniformRefinement();
      }
   }*/

   // 7. Define a parallel finite element space on the parallel mesh. Here we
   //    use the Raviart-Thomas finite elements of the specified order.
   FiniteElementCollection *hdiv_coll(new RT_FECollection(order+1, dim));
   FiniteElementCollection *l2_coll(new H1_FECollection(order));

   ParFiniteElementSpace *R_space = new ParFiniteElementSpace(pmesh, hdiv_coll);
   ParFiniteElementSpace *W_space = new ParFiniteElementSpace(pmesh, l2_coll);

   //Set up the problem
   real_t sig = 1.0;
   MemoryType mt = device.GetMemoryType();
   DarcyEMProblem demoProb(R_space, W_space, sig, mt, dim);
	
   //Set the solver and preconditioner
   //demoProb.BuildPreconditioner();
   demoProb.Set_Solver(verbose);

   //Solve the equations
   demoProb.Solve(verbose);

   //Visualise the results using ParaView
   double time = 0.0;
   demoProb.SetFields();
   ParaViewVisualise(demoProb.Fields, demoProb.FieldNames, order, pmesh, time);

   // 20. Free the used memory.
   delete W_space;
   delete R_space;
   delete l2_coll;
   delete hdiv_coll;
   delete pmesh;
   return 0;
}