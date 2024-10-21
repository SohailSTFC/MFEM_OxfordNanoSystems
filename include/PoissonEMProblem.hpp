#ifndef DARCYEMPROBLEM_HPP
#define DARCYEMPROBLEM_HPP 

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

#include "AnalyticSolution.hpp"
//
//
//  The problem class
//
//
class PoissonEMProblem
{
  protected:
    //Finite element spaces
    ParFiniteElementSpace *fespaceL=NULL; //Raviart Thomas elements

    //Block Matrix structure offsets
    Array<int> block_offsets;     // number of variables + 1 (2-variables J and v)
    Array<int> block_trueOffsets; // number of variables + 1 (2-variables J and v)

    //The permiativity of space
    real_t sigma;

    //The Bilinear forms of the block components
    // (Jacobian) (Assuming a symmetric Saddle point problem)
    ParBilinearForm *VVForm=NULL;

    //The linear forms of the block components
    // (residual)
    ParLinearForm   *VForm=NULL;

    //The complete block vectors
    mutable BlockVector x_vec, b_vec;
    mutable BlockVector tx_vec, tb_vec;

    // Form block operators (operates Matrix multiplication)
	// (This aggregates the block components of the forms)
    BlockOperator               *PoissonEMOp = NULL;
    BlockDiagonalPreconditioner *PoissonEMPr = NULL;

    //The Block hypre matrices and Transposes for Jacobian
    HypreParMatrix *M = NULL;

    //Shared pointer to the solver
    IterativeSolver* solver=NULL;

    //The Preconditioning objects
    HypreParMatrix *MinvBt = NULL;
    HypreParVector *Md = NULL;
    HypreParMatrix *S = NULL;
    Solver *invM=NULL, *invS=NULL;

  public:
    //The fields (Needed for postprocessing)
    //Made public for external access (not a great solution
	//but hey it works)
    vector<std::string>      FieldNames;
	vector<ParGridFunction*> Fields;


	//The constructor
    PoissonEMProblem(ParFiniteElementSpace *f1, ParFiniteElementSpace *f2
	             , real_t sig, MemoryType deviceMT, int dim);


    //Build and set a Preconditioner for the solver
    void BuildPreconditioner();

    //Read in and set the Boundary conditions

    //Set a linear/non-linear solver
    void Set_Solver(bool verbose);

    //Solve the equation
    void Solve(bool verbose);

    //Make the fields for post-processing
    void SetFields();

	//The destructor
    ~PoissonEMProblem();
};



//
//
//  Implementation of the problem class
//
//

//The constructor
//Constructs the problem and sets-up
//residual+jacobian operators/Forms
PoissonEMProblem::PoissonEMProblem(ParFiniteElementSpace *fL
                             , real_t sig, MemoryType deviceMT, int dim)
: sigma(sig)
{
  fespaceL  = new ParFiniteElementSpace(*fL);

  HYPRE_BigInt dimW = fespaceL->GlobalTrueVSize();

  std::cout << "***********************************************************\n";
  std::cout << "dim(R) = " << dimW << "\n";
  std::cout << "***********************************************************\n";


  // 8. Define the two BlockStructure of the problem.  block_offsets is used
  //    for Vector based on dof (like ParGridFunction or ParLinearForm),
  //    block_trueOffstes is used for Vector based on trueDof (HypreParVector
  //    for the rhs and solution of the linear system).  The offsets computed
  //    here are local to the processor.
  Array<int> bofs(2); // number of variables + 1
  bofs[0] = 0;
  bofs[1] = fespaceL->GetVSize();
  bofs.PartialSum();

  Array<int> btofs(2); // number of variables + 1
  btofs[0] = 0;
  btofs[1] = fespaceL->TrueVSize();
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
  VectorFunctionCoefficient ucoeff(dim, uFun_ex);
  FunctionCoefficient pcoeff(pFun_ex);


  //J-Linear-form of the equation  div grad V = q
  VForm = new ParLinearForm;
  VForm->Update(fespaceL, b_vec.GetBlock(1), 0);
  VForm->AddDomainIntegrator(new DiffusionIntegrator(gcoeff));
  VForm->Assemble();
  VForm->SyncAliasMemory(b_vec);
  VForm->ParallelAssemble(tb_vec.GetBlock(1));
  tb_vec.GetBlock(1).SyncAliasMemory(tb_vec);


  //The Bilinear block forms
  VVForm = new ParBilinearForm(fespaceL);

  //Set the integrators/integral forms and assemble the block matrices/bilinear forms
  VVForm->AddDomainIntegrator(new VectorFEDivergenceIntegrator);
  VVForm->Assemble();
  VVForm->Finalize();

  //Set the block bilinear/matrix operator
  Array<int> empty_tdof_list;  // empty

  PoissonEMOp = new BlockOperator(block_trueOffsets);
  M = VVForm->ParallelAssemble();
  PoissonEMOp->SetBlock(0,0, M);
};

//Builds a preconditioner needed to
//accelerate the Darcy problem solver (Optional)
void DarcyEMProblem::BuildPreconditioner()
{
   //Construct the a Schurr Complement
   //Gauss-Seidel block Preconditioner
   Md = new HypreParVector(MPI_COMM_WORLD, M->GetGlobalNumRows(),M->GetRowStarts());
   M->GetDiag(*Md);
   invM = new HypreDiagScale(*M);
   invM->iterative_mode = false;

   //Set the block diagonal bilinear/matrix preconditioning operator
   PoissonEMPr = new BlockDiagonalPreconditioner(block_trueOffsets);
   PoissonEMPr->SetDiagonalBlock(0, invM);


*PoissonEMOp = NULL;
    BlockDiagonalPreconditioner *PoissonEMPr = NULL;

}

//Sets the linear/non-linear solver
//for the Darcy problem
void PoissonEMProblem::Set_Solver( bool verbosity){
   int maxIter(500);
   real_t rtol(1.e-6);
   real_t atol(1.e-10);
   solver = new MINRESSolver(MPI_COMM_WORLD);
   solver->SetAbsTol(atol);
   solver->SetRelTol(rtol);
   solver->SetMaxIter(maxIter);
   solver->SetPrintLevel(verbosity);
  if(PoissonEMOp != NULL) solver->SetOperator(*PoissonEMOp);
  if(PoissonEMPr != NULL) solver->SetPreconditioner(*PoissonEMPr);
};

//Solves the system of equations
//for the Darcy problem
void PoissonEMProblem::Solve(bool verbosity){
  if(darcyEMOp != NULL){
    StopWatch chrono;
    chrono.Clear();
    chrono.Start();
    tx_vec = 0.0;
    solver->Mult(tb_vec, tx_vec);
    chrono.Stop();

    if (verbosity)
    {
       if(solver->GetConverged())
         std::cout << "MINRES converged in " << solver->GetNumIterations()
                   << " iterations with a residual norm of " << solver->GetFinalNorm() << ".\n";
       else
         std::cout << "MINRES did not converge in " << solver->GetNumIterations()
                   << " iterations. Residual norm is " << solver->GetFinalNorm() << ".\n";
         std::cout << "MINRES solver took " << chrono.RealTime() << "s. \n";
    }
  }else{
    if (verbosity) std::cout << "Error Darcy operator not built" << ".\n";
  }
};

//Setting-up/unpacking the fields for the Darcy problem
//these are needed before post-processing
void PoissonEMProblem::SetFields(){
  FieldNames.push_back("Potential");
  Fields.push_back(new ParGridFunction);
  Fields[0]->MakeRef(fespaceL, x_vec.GetBlock(0), 0);
  Fields[0->Distribute(&(tx_vec.GetBlock(0)));
};

	
PoissonEMProblem::~PoissonEMProblem(){
   if(VForm     != NULL) delete VForm;
   if(darcyEMOp != NULL) delete darcyEMOp;
   if(darcyEMPr != NULL) delete darcyEMPr;
   if(VVForm    != NULL) delete JJForm;
   if(M         != NULL) delete M;
   if(B         != NULL) delete B;
   if(Bt        != NULL) delete Bt;
   if(MinvBt    != NULL) delete MinvBt;
   if(Md        != NULL) delete Md;
   if(S         != NULL) delete S;
};
#endif