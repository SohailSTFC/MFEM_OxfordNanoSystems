//
// Author: Sohail Rathore
// Date  : 14/01/2025
//
// A Electromagnetic Operator that
// that uses a H(div)-L2-H1 formulation
// of the EM-problem to generate the
// residual and Jacobian
//
// The equation system is:
//  J + sigma*Grad(Vb) = 0
//  Div( J ) = 0
//  V - Vb = 0
//
#include "mfem.hpp"
#include <iostream>
#include "include/boundaryConditions.hpp"
#include "include/Visualisation.hpp"
#include "include/EMProblem_JV.hpp"

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Parse command-line options.
   const char *mesh_file = "mesh/OxNanoSys0.mesh";
   int  ref_levels = -1;
   int  order=1;
   bool par_format = false;
   bool pa = false;
   const char *device_config = "ceed-cpu"; //"cpu";//"ceed-cpu";
   bool visualization = 1;
   bool verbose = (myid == 0);

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
      int par_ref_levels = 3;
      for (int l = 0; l < par_ref_levels; l++) pmesh->UniformRefinement();
   }

   // 7. Define a parallel finite element space on the parallel mesh. Here we
   //    use the Hybrid finite elements of the specified order.
   FiniteElementCollection *J_coll (new RT_FECollection(order, dim));
   FiniteElementCollection *V_coll (new H1_FECollection(order, dim));

   vector<ParFiniteElementSpace*> feSpaces;
   feSpaces.push_back(new ParFiniteElementSpace(pmesh, J_coll));
   feSpaces.push_back(new ParFiniteElementSpace(pmesh, V_coll,1));

   Array<int> BOffsets(3), trueBOffsets(3);
   BOffsets[0] = 0;
   for(int I=0; I<feSpaces.size(); I++) BOffsets[I+1] = feSpaces[I]->GetVSize();
   BOffsets.PartialSum();

   trueBOffsets[0] = 0;
   for(int I=0; I<feSpaces.size(); I++) trueBOffsets[I+1] = feSpaces[I]->GetTrueVSize();
   trueBOffsets.PartialSum();

   // 8. Set the boundary
   //    and initial conditions
   //
   //
   int nTags=feSpaces[0]->GetMesh()->bdr_attributes.Max();
   Array<Array<int>*> BDR_markers(feSpaces.size());
   Array<int> J_BDRs(feSpaces[0]->GetMesh()->bdr_attributes.Max());
   Array<int> V_BDRs(feSpaces[1]->GetMesh()->bdr_attributes.Max());

   J_BDRs=1;
   V_BDRs=0;

   if(nTags>1) J_BDRs[0]=0;
   if(nTags>1) J_BDRs[1]=0;
   if(nTags>2) J_BDRs[2]=0;

   V_BDRs[0]=1;
   if(nTags>1) V_BDRs[1]=1;
   if(nTags>2) V_BDRs[2]=1;

   BDR_markers[0] = &J_BDRs;
   BDR_markers[1] = &V_BDRs;

   // 9. Define Gridfunctions and initial conditions, BCS and 
   //
   //
   BlockVector x_vec(BOffsets), tx_vec(trueBOffsets), tb_vec(trueBOffsets);
   vector<ParGridFunction*> JV;
   vector<string> VarNames;
 
   for(int I=0; I<feSpaces.size(); I++){
     JV.push_back(new ParGridFunction);
     JV[I]->MakeRef(feSpaces[I], x_vec.GetBlock(I), 0);
     JV[I]->Distribute(&(tx_vec.GetBlock(I)));
   }
   VarNames.push_back("J_Field");
   VarNames.push_back("V_H1_Potential");

   tx_vec = 0.0;
   tb_vec = 0.0;
   tx_vec.GetBlock(1) = 0.1;
   
   //Set the BCs from the coefficient funcs
   Vector ZeroVec(dim); ZeroVec = 0.0;
   VectorConstantCoefficient ZeroV(ZeroVec);
   ConstantCoefficient ThreeS(3.0);
   ConstantCoefficient ZeroS(0.0);

   JV[0]->ProjectBdrCoefficient(ZeroV,*BDR_markers[0]);

   Array<int> V_BDRs_tmp(nTags);
   V_BDRs_tmp = 0;
   if(nTags>1) V_BDRs_tmp[1]=1;
   if(nTags>2) V_BDRs_tmp[2]=1;
   JV[1]->ProjectBdrCoefficient(ThreeS,V_BDRs_tmp);
   V_BDRs_tmp = 0;
   V_BDRs_tmp[0]=1;
   JV[1]->ProjectBdrCoefficient(ZeroS,V_BDRs_tmp);


   Array<int> ess_Jtdof, ess_Vtdof;
   feSpaces[0]->GetEssentialTrueDofs(*BDR_markers[0], ess_Jtdof);
   feSpaces[1]->GetEssentialTrueDofs(*BDR_markers[1], ess_Vtdof);

   applyDirchValues(*JV[0], tx_vec.GetBlock(0), ess_Jtdof);
   applyDirchValues(*JV[1], tx_vec.GetBlock(1), ess_Vtdof);

   // 10. Build  the Problem-Operator
   //     class
   //
   int PSize=0;
   for(int I=0; I<feSpaces.size(); I++) PSize += feSpaces[I]->GetTrueVSize();
   MemoryType mt = device.GetMemoryType();
   EMOperator sampleProb(feSpaces, trueBOffsets, BDR_markers, dim, mt, PSize);

   // 11. Build the Non-linear 
   //     system solver
   //
   double rel_tol = 1.0e-15;
   double abs_tol = 1.0e-12;
   int maxIter = 500;
   //FGMRESSolver *solver = new FGMRESSolver(MPI_COMM_WORLD);
   MINRESSolver *solver = new MINRESSolver(MPI_COMM_WORLD);
   solver->SetAbsTol(abs_tol);
   solver->SetRelTol(rel_tol);
   solver->SetMaxIter(maxIter);
   solver->SetPrintLevel(verbose);
   solver->SetOperator(sampleProb);
  // solver->SetPreconditioner(sampleProb.Precon());
   solver->Mult(tb_vec, tx_vec);

   //12. Output and visualise the data
   //
   //
   for(int I=0; I<JV.size(); I++) JV[I]->Distribute(&(tx_vec.GetBlock(I)));
   ParaViewVisualise("EMProblem2", JV, VarNames, order, pmesh, 0.0);

   //13. Clean up and
   //    deallocation
   //

   return 0;
};