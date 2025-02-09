//
// Author: Sohail Rathore
// Date  : 14/01/2025
//
// A Electromagnetic Operator that
// that uses a H1 formulation
// of the EM-problem to generate the
// residual and Jacobian
//
// The equation system is:
//  Div( sigma*Grad(Vb) ) = F
//
#include "mfem.hpp"
#include <iostream>
#include "include/boundaryConditions.hpp"
#include "include/Visualisation.hpp"
#include "include/EMProblem_V.hpp"

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
   int  order=2;
   bool par_format = false;
   bool pa = false;
   const char *device_config = "cpu"; //"cpu";//"ceed-cpu";
   bool visualization = 1;
   bool verbose = (myid == 0);

   //Boundary condition arrays
   Array<int> dbcs;
   Array<int> nbcs;
   Vector dbcv;
   Vector nbcv;


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

   //Boundary condition surface
   //from command line arguments
   //Dirch BCs (V)
   args.AddOption(&dbcs, "-dbcs", "--dirichlet-bc-surf",
                  "Dirichlet Boundary Condition Surfaces");
   args.AddOption(&dbcv, "-dbcv", "--dirichlet-bc-vals",
                  "Dirichlet Boundary Condition Values");

   //Neumann BCs (Dirch on J)
   args.AddOption(&nbcs, "-nbcs", "--neumann-bc-surf",
                  "Neumann Boundary Condition Surfaces");
   args.AddOption(&nbcv, "-nbcv", "--neumann-bc-vals",
                  "Neumann Boundary Condition Values");
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
   int dim  = mesh->Dimension();
   int sdim = mesh->SpaceDimension();

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
      int par_ref_levels = 2;
      for (int l = 0; l < par_ref_levels; l++) pmesh->UniformRefinement();
   }

   // 7. Define a parallel finite element space on the parallel mesh. Here we
   //    use the Hybrid finite elements of the specified order.
   //    Two subsequent Orders are used to estimate error
   FiniteElementCollection *V_c_coll (new H1_FECollection(order,  dim));//Coarse p-Order-element
   FiniteElementCollection *V_f_coll (new H1_FECollection(order+1,dim));//Fine p-Order-element

   vector<ParFiniteElementSpace*> feSpacesF, feSpacesC;
   feSpacesC.push_back(new ParFiniteElementSpace(pmesh, V_c_coll, 1));
   feSpacesF.push_back(new ParFiniteElementSpace(pmesh, V_f_coll, 1));

   //Fine space block offsets
   Array<int> BOffsetsF(2), trueBOffsetsF(2);
   BOffsetsF[0] = 0;
   for(int I=0; I<feSpacesF.size(); I++) BOffsetsF[I+1] = feSpacesF[I]->GetVSize();
   BOffsetsF.PartialSum();

   trueBOffsetsF[0] = 0;
   for(int I=0; I<feSpacesF.size(); I++) trueBOffsetsF[I+1] = feSpacesF[I]->GetTrueVSize();
   trueBOffsetsF.PartialSum();

   //Coarse space Offsets
   Array<int> BOffsetsC(2), trueBOffsetsC(2);
   BOffsetsC[0] = 0;
   for(int I=0; I<feSpacesC.size(); I++) BOffsetsC[I+1] = feSpacesC[I]->GetVSize();
   BOffsetsC.PartialSum();

   trueBOffsetsC[0] = 0;
   for(int I=0; I<feSpacesC.size(); I++) trueBOffsetsC[I+1] = feSpacesC[I]->GetTrueVSize();
   trueBOffsetsC.PartialSum();


   // 9. Set the boundary
   //    and initial conditions
   //
   //
   int nTagsF=feSpacesF[0]->GetMesh()->bdr_attributes.Max();
   int nTagsC=feSpacesC[0]->GetMesh()->bdr_attributes.Max();
   Array<Array<int>*> BDR_markers(feSpacesF.size());
   Array<int> V_BDRs(nTagsF);

   V_BDRs=0;
   V_BDRs[0]=1;
   if(nTagsF>1) V_BDRs[1]=1;
   if(nTagsF>2) V_BDRs[2]=1;

   BDR_markers[0] = &V_BDRs;


   // 8. Define Gridfunctions and initial conditions, BCS and 
   //
   //
   BlockVector x_vec(BOffsetsF), tx_vec(trueBOffsetsF), tb_vec(trueBOffsetsF);
   vector<ParGridFunction*> V;
   vector<string> VarNames;
 
   for(int I=0; I<feSpacesF.size(); I++){
     V.push_back(new ParGridFunction);
     V[I]->MakeRef(feSpacesF[I], x_vec.GetBlock(I), 0);
     V[I]->Distribute(&(tx_vec.GetBlock(I)));
   }

   VarNames.push_back("V-H1-Potential");
   x_vec  = 0.0;
   tx_vec = 0.0;
   tb_vec = 0.0;
   tx_vec.GetBlock(0) = 0.2;
   x_vec.GetBlock(0)  = 0.2;

   //Set the BCs from the coefficient funcs
   ConstantCoefficient ThreeS(3.0);
   ConstantCoefficient ZeroS(0.0);

   Array<int> V_BDRs_tmp(nTagsF);
   V_BDRs_tmp = 0;
   if(nTagsF>1) V_BDRs_tmp[1]=1;
   if(nTagsF>2) V_BDRs_tmp[2]=1;
   V[0]->ProjectBdrCoefficient(ThreeS,V_BDRs_tmp);
   V_BDRs_tmp = 0;
   V_BDRs_tmp[0]=1;
   V[0]->ProjectBdrCoefficient(ZeroS,V_BDRs_tmp);

   Array<int> ess_Jtdof, ess_Vtdof;
   feSpacesF[0]->GetEssentialTrueDofs(*BDR_markers[0], ess_Vtdof);

   applyDirchValues(*V[0], tb_vec.GetBlock(0), ess_Vtdof);
   applyDirchValues(*V[0], tx_vec.GetBlock(0), ess_Vtdof);


   // 10. Build  the Problem-Operator
   //     class
   //
   int PSize=0;
   for(int I=0; I<feSpacesF.size(); I++) PSize += feSpacesF[I]->GetTrueVSize();
   MemoryType mt = device.GetMemoryType();
   EMOperatorV sampleProbF(feSpacesF, BDR_markers, dim, mt, PSize); //Fine problem
  // EMOperatorV sampleProbC(feSpacesC, BDR_markers, dim, mt, PSize); //Coarse problem

   //11.1 Build an error estimator
   //     Using-estimator from
   //     example 6p
   //
   L2_FECollection flux_fec(order, dim);
   ParFiniteElementSpace flux_fes(pmesh, &flux_fec, sdim);
   FiniteElementCollection *smooth_flux_fec = NULL;
   ParFiniteElementSpace *smooth_flux_fes = NULL;

   // A possible option for the smoothed flux space: H1^dim space
   // the AMR technique minimises the L2-norm of a particular grid
   // functional, the way it does this is via a method of projection
   // thus access is needed to Bilinear form that lives inside the
   // problem class
   //
/*
   smooth_flux_fec = new H1_FECollection(order, dim);
   smooth_flux_fes = new ParFiniteElementSpace(pmesh, smooth_flux_fec, dim);
 //  LpErrorEstimator estimator(2, Coefficient &coef, GridFunction &sol); 
 //  L2ZienkiewiczZhuEstimator estimator(*integ, V[0], flux_fes, *smooth_flux_fes);


   //The refiner used here uses a thresholding strategy 
   //The strategy here is to refine elements with errors larger than a
   //fraction of the maximum element error.
   ThresholdRefiner refiner(estimator);
   refiner.SetTotalErrorFraction(0.7);
*/


   //Construct and AMR loop
   int Max_AMR_its=1;
   for(int AMR_it=0; AMR_it<Max_AMR_its; AMR_it++){
     //Reassemble the System
     sampleProbF.ReassembleSystem();

     //Build the solver and
     //solve the linear system
     double rel_tol = 1.0e-15;
     double abs_tol = 1.0e-12;
     int maxIter = 500;
         /**FGMRESSolver *solver = new FGMRESSolver(MPI_COMM_WORLD);**/
     MINRESSolver *solver = new MINRESSolver(MPI_COMM_WORLD);
     solver->SetAbsTol(abs_tol);
     solver->SetRelTol(rel_tol);
     solver->SetMaxIter(maxIter);
     solver->SetPrintLevel(verbose);
     solver->SetOperator(sampleProbF);
     solver->Mult(tb_vec, tx_vec);
     delete solver; //delete once solved

     //Apply the refiner to
     //the finite element mesh
   //  refiner.Apply(*pmesh);


     //Update the feSpaces and grid functions
     //so that the problem is resolved
   }


   //12. Output and visualise the data
   //
   //
   for(int I=0; I<V.size(); I++) V[I]->Distribute(&(tx_vec.GetBlock(I)));
   ParaViewVisualise("EMProblem1", V, VarNames, order, pmesh, 0.0);


   //13. Clean up and
   //    deallocation
   //
   return 0;
};