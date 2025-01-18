
#include "mfem.hpp"
#include <iostream>
#include "include/boundaryConditions.hpp"
#include "include/Visualisation.hpp"
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
   const char *mesh_file = "mesh/star.mesh";
   int  ref_levels = -1;
   int  order=3;
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
   //    use the Taylor-Hood finite elements of the specified order.
   vector<int> orders;             
   orders.push_back(order);
   orders.push_back(order-1);
   FiniteElementCollection *vel_coll(new H1_FECollection(orders[0], dim));
   FiniteElementCollection *prs_coll(new H1_FECollection(orders[1], dim));

   vector<ParFiniteElementSpace*> feSpaces;
   feSpaces.push_back(new ParFiniteElementSpace(pmesh, vel_coll, dim));
   feSpaces.push_back(new ParFiniteElementSpace(pmesh, prs_coll, 1));

   Array<int> BOffsets(3), trueBOffsets(3);
   BOffsets[0] = 0;
   for(int I=0; I<feSpaces.size(); I++) BOffsets[I+1] = feSpaces[I]->GetVSize();
   BOffsets.PartialSum();

   trueBOffsets[0] = 0;
   for(int I=0; I<feSpaces.size(); I++) trueBOffsets[I+1] = feSpaces[I]->GetTrueVSize();
   trueBOffsets.PartialSum();


   // 8. Define Gridfunctions and initial conditions, BCS and 
   //
   //
   BlockVector x_vec(BOffsets), tx_vec(trueBOffsets), tb_vec(trueBOffsets);
   vector<ParGridFunction*> U;
   vector<string> VarNames;
 
   for(int I=0; I<feSpaces.size(); I++){
     U.push_back(new ParGridFunction);
     U[I]->MakeRef(feSpaces[I], x_vec.GetBlock(I), 0);
     U[I]->Distribute(&(tx_vec.GetBlock(I)));
   }
   VarNames.push_back("Velocity");
   VarNames.push_back("Pressure");


   // 9. Build  the Problem-Operator
   //    class
   //
   int PSize=0;
   for(int I=0; I<feSpaces.size(); I++) PSize += feSpaces[I]->GetTrueVSize();
   MemoryType mt = device.GetMemoryType();
   SNSSOperator sampleProb(feSpaces, trueBOffsets, dim, mt, PSize);


   // 10. Set the boundary
   //     conditions
   //
   //
   Array<Array<int>*> BDR_markers(feSpaces);
   Array<int> uvel_BDRs(feSpaces[0]->GetMesh()->bdr_attributes.Max());
   Array<int> pres_BDRs(feSpaces[1]->GetMesh()->bdr_attributes.Max());

   uvel_BDRs=0;
   pres_BDRs=0;
   uvel_BDRs[0]=1;
 //  pres_BDRs[0]=0;

   BDR_markers[0] = &uvel_BDRs;
   BDR_markers[1] = &pres_BDRs;

   // 11. Build the Non-linear 
   //     Newton-Rhapson solver
   //
   NewtonSolver newton_solver;/*
   newton_solver.iterative_mode = true;
   newton_solver.SetSolver(*j_solver);
   newton_solver.SetOperator(*this);
   newton_solver.SetPrintLevel(-1);
   newton_solver.SetMonitor(newton_monitor);
   newton_solver.SetRelTol(rel_tol);
   newton_solver.SetAbsTol(abs_tol);
   newton_solver.SetMaxIter(iter);*/
   tb_vec = 1.0;
   tx_vec = 1.0;
   sampleProb.Mult(tb_vec,tx_vec);
   for(int I=0; I<U.size(); I++) U[I]->Distribute(&(tx_vec.GetBlock(I)));

   //12. Output and visualise the data
   //
   //
   //for(int I=0; I<U.size(); I++) U[I]->Distribute(&(tx_vec.GetBlock(I)));
   ParaViewVisualise("NavierStokes", U, VarNames, order, pmesh, 0.0);


   //12. Clean up
   //
   //
  

   return 0;
};