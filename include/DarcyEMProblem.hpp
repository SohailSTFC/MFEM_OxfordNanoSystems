#ifndef DARCYEMPROBLEM_HPP
#define DARCYEMPROBLEM_HPP 

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

#include "BoundaryAndInitialSolution.hpp"
//
//
//  The problem class
//
//
class DarcyEMProblem
{
  protected:
    //Problem dimension
    int DIM;

    //Finite element spaces
    ParFiniteElementSpace *fespaceRT=NULL; //Raviart Thomas elements
    ParFiniteElementSpace *fespaceL=NULL;  //Lagrange finite element space

    //Block Matrix structure offsets
    Array<int> block_offsets;     // number of variables + 1 (2-variables J and v)
    Array<int> block_trueOffsets; // number of variables + 1 (2-variables J and v)

    //The permiativity of space
    real_t sigma;

    //The Bilinear forms of the block components
    // (Jacobian) (Assuming a symmetric Saddle point problem)
    ParBilinearForm      *JJForm=NULL, *VVForm=NULL;
    ParMixedBilinearForm *JVForm=NULL, *VJForm=NULL;

    //The linear forms of the block components
    // (residual)
    ParLinearForm   *JForm=NULL;
    ParLinearForm   *VForm=NULL;

    //The complete block vectors
    mutable BlockVector x_vec, b_vec;
    mutable BlockVector tx_vec, tb_vec;

    // Form block operators (operates Matrix multiplication)
	// (This aggregates the block components of the forms)
    BlockOperator               *darcyEMOp = NULL;
    BlockDiagonalPreconditioner *darcyEMPr = NULL;

    //The Block hypre matrices and Transposes for Jacobian
    TransposeOperator *Bt = NULL;
    OperatorPtr opM, opB, opC, opD;
    Vector Md_PA;

    //Shared pointer to the solver
    IterativeSolver* solver=NULL;

    //The Preconditioning objects
    HypreParMatrix *MinvBt = NULL;
    HypreParVector *Md = NULL;
    HypreParMatrix *S = NULL;
    Solver *invM=NULL, *invS=NULL;

    //Boundary Conditions
    vector<vector<double>> DirchVal;       //Dirchelet value of BC
    Array<int> ess_tdof_J, ess_tdof_v;     //Dirchelet BC DOF's
    Array<int> ess_bdr_J, ess_bdr_v;       //Dirchelet BC Tag's
  
    //Read in and set the Boundary conditions
    void SetBCsArrays();

    //Read in and set the Boundary conditions
    void SetFieldBCs();

  public:
    //The fields (Needed for postprocessing)
    //Made public for external access (not a great solution
	//but hey it works)
    vector<std::string>      FieldNames;
	vector<ParGridFunction*> Fields;

	//The constructor
    DarcyEMProblem(ParFiniteElementSpace *f1, ParFiniteElementSpace *f2
	             , real_t sig, MemoryType deviceMT, int dim);

    //Build and set a Preconditioner for the solver
    void BuildPreconditioner();

    //Set a linear/non-linear solver
    void Set_Solver(bool verbose);

    //Solve the equation
    void Solve(bool verbose);

    //Make the fields for post-processing
    void SetFields();

	//The destructor
    ~DarcyEMProblem();
};



//
//
//  Implementation of the problem class
//
//

//The constructor
//Constructs the problem and sets-up
//residual+jacobian operators/Forms
DarcyEMProblem::DarcyEMProblem(ParFiniteElementSpace *f1RT
                             , ParFiniteElementSpace *f2L
                             , real_t sig, MemoryType deviceMT, int dim)
: sigma(sig)
{
  DIM = dim;
  fespaceRT = new ParFiniteElementSpace(*f1RT);
  fespaceL  = new ParFiniteElementSpace(*f2L);

  HYPRE_BigInt dimR = fespaceRT->GlobalTrueVSize();
  HYPRE_BigInt dimW = fespaceL->GlobalTrueVSize();

  std::cout << "***********************************************************\n";
  std::cout << "dim(R) = " << dimR << "\n";
  std::cout << "dim(W) = " << dimW << "\n";
  std::cout << "dim(R+W) = " << dimR + dimW << "\n";
  std::cout << "***********************************************************\n";

  //Setting the boundary conditions
  SetBCsArrays();


  // Get the block offsets and true block offsets
  // to construct the block structured vectors
  //and matrix operators
  Array<int> bofs(3); // number of variables + 1
  bofs[0] = 0;
  bofs[1] = fespaceRT->GetVSize();
  bofs[2] = fespaceL->GetVSize();
  bofs.PartialSum();

  Array<int> btofs(3); // number of variables + 1
  btofs[0] = 0;
  btofs[1] = fespaceRT->TrueVSize();
  btofs[2] = fespaceL->TrueVSize();
  btofs.PartialSum();
 
  block_offsets     = Array<int>(bofs);
  block_trueOffsets = Array<int>(btofs);


  // 10. Define the parallel grid function and parallel linear forms, solution
  //     vector and rhs.
  x_vec.Update  (block_offsets, deviceMT);
  b_vec.Update  (block_offsets, deviceMT);
  tx_vec.Update (block_trueOffsets, deviceMT);
  tb_vec.Update (block_trueOffsets, deviceMT);


  // 9. Define the coefficients, analytical solution, and rhs of the PDE.
  // the coefficients and functions
  ConstantCoefficient k(1.0);
  VectorFunctionCoefficient fcoeff(dim, fFun);
  FunctionCoefficient fnatcoeff(f_natural);
  FunctionCoefficient gcoeff(gFun);

  //
  //  The linear forms (The RHS/Residual forms)
  //
  //v-Linear-form of the equation J + grad v = Je
  VForm = new ParLinearForm;
  VForm->Update(fespaceRT, b_vec.GetBlock(0), 0);
  VForm->AddDomainIntegrator(new VectorFEDomainLFIntegrator(fcoeff));
  VForm->AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator(fnatcoeff));
  VForm->Assemble();
  VForm->SyncAliasMemory(b_vec);
  VForm->ParallelAssemble(tb_vec.GetBlock(0));
  tb_vec.GetBlock(0).SyncAliasMemory(tb_vec);

  //J-Linear-form of the equation  div J = q
  JForm = new ParLinearForm;
  JForm->Update(fespaceL, b_vec.GetBlock(1), 0);
  JForm->AddDomainIntegrator(new DomainLFIntegrator(gcoeff));
  JForm->Assemble();
  JForm->SyncAliasMemory(b_vec);
  JForm->ParallelAssemble(tb_vec.GetBlock(1));
  tb_vec.GetBlock(1).SyncAliasMemory(tb_vec);


  //
  // The Bilinear forms (matrix/jacobian forms)
  //
  Array<int> ess_tdof_empty;

  //The Bilinear block forms
  JJForm = new ParBilinearForm(fespaceRT);
  JVForm = new ParMixedBilinearForm(fespaceRT, fespaceL);
  VJForm = new ParMixedBilinearForm(fespaceRT, fespaceL);
  VVForm = new ParBilinearForm(fespaceL);

  //Set the assembly level
  JJForm->SetAssemblyLevel(AssemblyLevel::PARTIAL);
  JVForm->SetAssemblyLevel(AssemblyLevel::PARTIAL);
  VJForm->SetAssemblyLevel(AssemblyLevel::PARTIAL);
  VVForm->SetAssemblyLevel(AssemblyLevel::PARTIAL);

  //Set the integrators/integral forms and assemble the block matrices/bilinear forms
  JJForm->AddDomainIntegrator(new VectorFEMassIntegrator(k));
  JJForm->Assemble();

  JVForm->AddDomainIntegrator(new VectorFEDivergenceIntegrator);
  JVForm->Assemble();

  VJForm->AddDomainIntegrator(new VectorFEDivergenceIntegrator);
  VJForm->Assemble();

  VVForm->Assemble();

  //Set the BCs and finalize the bilinear forms
  JJForm->FormSystemMatrix( ess_tdof_J, opM);
  VJForm->FormRectangularSystemMatrix( ess_tdof_J, ess_tdof_empty, opB);
  JVForm->FormRectangularSystemMatrix( ess_tdof_empty, ess_tdof_v, opC);
  Bt = new TransposeOperator(opB.Ptr());
  VVForm->FormSystemMatrix( ess_tdof_v, opD);

  //Set the block matrix operator
  darcyEMOp = new BlockOperator(block_trueOffsets);
  darcyEMOp->SetBlock(0,0, opM.Ptr());
  darcyEMOp->SetBlock(0,1, Bt, -1.0);
  darcyEMOp->SetBlock(1,0, opC.Ptr(), -1.0);
  darcyEMOp->SetBlock(1,1, opD.Ptr());
};

//Sets the natural and essential boundary
//conditions
void DarcyEMProblem::SetBCsArrays(){
  //Essential boundary condition tags
  ess_bdr_J = Array<int>(fespaceRT->GetMesh()->bdr_attributes.Max());
  ess_bdr_v = Array<int>(fespaceL->GetMesh()->bdr_attributes.Max());

  //initialise the arrays
  ess_bdr_J = 0;
  ess_bdr_v = 0;

  //fixed J
  ess_bdr_J[3] = 1;

  //fixed v
  ess_bdr_v[0] = 1;
  ess_bdr_v[1] = 1;
  ess_bdr_v[2] = 1;

  //Find the True Dofs
  fespaceRT->GetEssentialTrueDofs(ess_bdr_J, ess_tdof_J);
  fespaceL->GetEssentialTrueDofs(ess_bdr_v, ess_tdof_v);

  cout << setw(10) << "RT elements: " << setw(10) << ess_tdof_J.Size() << "\n";
  cout << setw(10) << "H1 elements: " << setw(10) << ess_tdof_v.Size() << "\n";

  //The Dirchelet BC values
  vector<double> DirchVal_tmp;
  DirchVal.clear();
  DirchVal_tmp.clear();
  
  //J-Field BC-values
  DirchVal_tmp.push_back(0.00); // N/A
  DirchVal_tmp.push_back(0.00); // N/A
  DirchVal_tmp.push_back(0.00); // N/A
  DirchVal_tmp.push_back(0.00); // n.J = 0
  DirchVal.push_back(DirchVal_tmp);

  //v-Field BC-values
  DirchVal_tmp[0] = 0.00; // Fixed v = c
  DirchVal_tmp[1] = 3.00; // Fixed v = c
  DirchVal_tmp[2] = 3.00; // Fixed v = c
  DirchVal_tmp[3] = 0.00; // N/A
  DirchVal.push_back(DirchVal_tmp);
  DirchVal_tmp.clear();
};

//Set the BCs by setting the solution vectors and
void DarcyEMProblem::SetFieldBCs(){
  //Set the boundary Solution function for
  //the current field
  int nJ_tags = fespaceRT->GetMesh()->bdr_attributes.Max();
  int nv_tags = fespaceL->GetMesh()->bdr_attributes.Max();

  tb_vec = 0.0;
  b_vec  = 0.0;
  tx_vec = 0.0;
  x_vec  = 0.0;
  x_vec.GetBlock(1) = 0.2;
  tx_vec.GetBlock(1) = 0.2;

  //Set the J-Field BCs by looping over the
  //active boundaries
  cout << setw(10) << "RT elements: " << setw(10) << ess_tdof_J.Size() << "\n";
  for(int I=0; I<nJ_tags; I++){
	int K = ess_bdr_J[I];
    if(K == 1){ //Checks if boundary is active
      Array<int> ess_tdof, ess_bdr_tmp(nJ_tags);
      ess_bdr_tmp = 0;
      ess_bdr_tmp[I] = 1;
      fespaceRT->GetEssentialTrueDofs(ess_bdr_tmp, ess_tdof);
      cout << setw(10) << I << setw(10) << ess_tdof.Size() << "\n";
      x_vec.GetBlock(0).SetSubVector( ess_tdof, DirchVal[0][I] );
      b_vec.GetBlock(0).SetSubVector( ess_tdof, DirchVal[0][I] );
      tx_vec.GetBlock(0).SetSubVector( ess_tdof, DirchVal[0][I] );
      tb_vec.GetBlock(0).SetSubVector( ess_tdof, DirchVal[0][I] );
    }
  }

  //Set the J-Field BCs by looping over the
  //active boundaries
  cout << setw(10) << "H1 elements: " << setw(10) << ess_tdof_v.Size() << "\n";
  for(int I=0; I<nv_tags; I++){
	int K = ess_bdr_v[I];
    if(K == 1){ //Checks if boundary is active
      Array<int> ess_tdof, ess_bdr_tmp(nv_tags);
      ess_bdr_tmp = 0;
      ess_bdr_tmp[I] = 1;
      fespaceL->GetEssentialTrueDofs(ess_bdr_tmp, ess_tdof);
      cout << setw(10) << I << setw(10) << ess_tdof.Size() << "\n";
      x_vec.GetBlock(1).SetSubVector( ess_tdof, DirchVal[1][I] );
      b_vec.GetBlock(1).SetSubVector( ess_tdof, DirchVal[1][I] );
      tx_vec.GetBlock(1).SetSubVector( ess_tdof, DirchVal[1][I] );
      tb_vec.GetBlock(1).SetSubVector( ess_tdof, DirchVal[1][I] );
    }
  }
};


//Builds a preconditioner needed to
//accelerate the Darcy problem solver (Optional)
void DarcyEMProblem::BuildPreconditioner()
{
   //Construct the a Schurr Complement
   //Gauss-Seidel block Preconditioner
   Md_PA.SetSize(fespaceRT->GetTrueVSize());
   JJForm->AssembleDiagonal(Md_PA);
   auto Md_host = Md_PA.HostRead();
   Vector invMd(Md_PA.Size());
   for (int i=0; i<Md_PA.Size(); ++i) invMd(i) = 1.0 / Md_host[i];

   Vector BMBt_diag(fespaceL->GetTrueVSize());
   BMBt_diag = 1.0;
   VJForm->AssembleDiagonal_ADAt(invMd, BMBt_diag);
   invM = new OperatorJacobiSmoother(Md_PA, ess_tdof_J);
   invS = new OperatorJacobiSmoother(BMBt_diag, ess_tdof_v);

   invM->iterative_mode = false;
   invS->iterative_mode = false;
   darcyEMPr = new BlockDiagonalPreconditioner(block_trueOffsets);
   darcyEMPr->SetDiagonalBlock(0, invM);
   darcyEMPr->SetDiagonalBlock(1, invS);
}

//Sets the linear/non-linear solver
//for the Darcy problem
void DarcyEMProblem::Set_Solver( bool verbosity){
  int maxIter(500);
  real_t rtol(1.e-6);
  real_t atol(1.e-10);
  solver = new MINRESSolver(MPI_COMM_WORLD);
  solver->SetAbsTol(atol);
  solver->SetRelTol(rtol);
  solver->SetMaxIter(maxIter);
  solver->SetPrintLevel(verbosity);
  if(darcyEMOp != NULL) solver->SetOperator(*darcyEMOp);
  if(darcyEMPr != NULL) solver->SetPreconditioner(*darcyEMPr);
};

//Solves the system of equations
//for the Darcy problem
void DarcyEMProblem::Solve(bool verbosity){
  if(darcyEMOp != NULL){
    StopWatch chrono;
    chrono.Clear();
    chrono.Start();
    SetFieldBCs();
    solver->Mult(tb_vec, tx_vec);
    chrono.Stop();

    if (verbosity)
    {
      std::cout << "MINRES ended in "                     << solver->GetNumIterations()
                << " iterations with a residual norm of " << solver->GetFinalNorm() << ".\n";
      std::cout << "MINRES solver took "                  << chrono.RealTime()      << "s. \n";
    }
  }else{
    if (verbosity) std::cout << "Error Darcy operator not built" << ".\n";
  }
};

//Setting-up/unpacking the fields for the Darcy problem
//these are needed before post-processing
void DarcyEMProblem::SetFields(){
  FieldNames.push_back("Current");
  Fields.push_back(new ParGridFunction);
  Fields[0]->MakeRef(fespaceRT, x_vec.GetBlock(0), 0);
  Fields[0]->Distribute(&(tx_vec.GetBlock(0)));

  FieldNames.push_back("Potential");
  Fields.push_back(new ParGridFunction);
  Fields[1]->MakeRef(fespaceL, x_vec.GetBlock(1), 0);
  Fields[1]->Distribute(&(tx_vec.GetBlock(1)));
};


DarcyEMProblem::~DarcyEMProblem(){
   if(JForm     != NULL) delete JForm;
   if(VForm     != NULL) delete VForm;
   if(darcyEMOp != NULL) delete darcyEMOp;
   if(darcyEMPr != NULL) delete darcyEMPr;
   if(JJForm    != NULL) delete JJForm;
   if(JVForm    != NULL) delete JVForm;
   if(VJForm    != NULL) delete VJForm;
   if(VVForm    != NULL) delete VVForm;
};
#endif