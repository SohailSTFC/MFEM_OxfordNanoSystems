
#include "mfem.hpp"
#include <iostream>
#include "include/SNSSolver.hpp"

using namespace mfem;
using namespace std;


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
   const char *device_config = "cpu"; //"cpu";//"ceed-cpu";
   bool visualization = 1;

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

  args.Parse();
  if (!args.Good()and(verbose)) args.PrintUsage(cout);
  if (!args.Good()) return 1;
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
    if (ref_levels == -1) ref_levels = (int)floor(log(10000./mesh->GetNE())/log(2.)/dim);
    for (int l = 0; l < ref_levels; l++) mesh->UniformRefinement();

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
   FiniteElementCollection *vel_coll(new H1_FECollection(order, dim));
   FiniteElementCollection *prs_coll(new H1_FECollection(order-1, dim));

   ParFiniteElementSpace *R_space = new ParFiniteElementSpace(pmesh, vel_coll);
   ParFiniteElementSpace *W_space = new ParFiniteElementSpace(pmesh, prs_coll);

   // 8. Define Gridfunctions and initial conditions, BCS and 
   //    
   //
   ParGridFunction U, P;
   ParGridFunction Zero;


  return 0;
}
