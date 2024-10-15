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
class DarcyEMProblem
{
  protected:
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
    ParBilinearForm      *JJForm=NULL;
    ParMixedBilinearForm *JVForm=NULL;

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
    HypreParMatrix *M = NULL;
    HypreParMatrix *B = NULL;
    TransposeOperator *Bt = NULL;

    //Shared pointer to the solver
    IterativeSolver* solver=NULL;

    //The Preconditioning objects
    HypreParMatrix *MinvBt = NULL;
    HypreParVector *Md = NULL;
    HypreParMatrix *S = NULL;
    Solver *invM=NULL, *invS=NULL;

  public:
	//The constructor
    DarcyEMProblem(ParFiniteElementSpace *f1, ParFiniteElementSpace *f2
	             , real_t sig, MemoryType deviceMT, int dim);


    //Build and set a Preconditioner for the solver
    void BuildPreconditioner();

    //Set a linear/non-linear solver
    void Set_Solver(bool verbose);

    //Solve the equation
    void Solve(bool verbose);

	//The destructor
    ~DarcyEMProblem();
};



//
//
//  Implementation of the problem class
//
//

//The constructor
DarcyEMProblem::DarcyEMProblem(ParFiniteElementSpace *f1RT
                             , ParFiniteElementSpace *f2L
                             , real_t sig, MemoryType deviceMT, int dim)
: fespaceRT(f1RT), fespaceL(f2L), sigma(sig)
{
  

  HYPRE_BigInt dimR = fespaceRT->GlobalTrueVSize();
  HYPRE_BigInt dimW = fespaceL->GlobalTrueVSize();

  std::cout << "***********************************************************\n";
  std::cout << "dim(R) = " << dimR << "\n";
  std::cout << "dim(W) = " << dimW << "\n";
  std::cout << "dim(R+W) = " << dimR + dimW << "\n";
  std::cout << "***********************************************************\n";


  // 8. Define the two BlockStructure of the problem.  block_offsets is used
  //    for Vector based on dof (like ParGridFunction or ParLinearForm),
  //    block_trueOffstes is used for Vector based on trueDof (HypreParVector
  //    for the rhs and solution of the linear system).  The offsets computed
  //    here are local to the processor.
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
  VectorFunctionCoefficient ucoeff(dim, uFun_ex);
  FunctionCoefficient pcoeff(pFun_ex);


  //DarcyEMProblem demoProb(fespaceRT, W_space, sig);
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


  //The Bilinear block forms
  JJForm = new ParBilinearForm(fespaceRT);
  JVForm = new ParMixedBilinearForm(fespaceRT, fespaceL);

  //Set the integrators/integral forms and assemble the block matrices/bilinear forms
  JJForm->AddDomainIntegrator(new VectorFEMassIntegrator(k));
  JJForm->Assemble();
  JJForm->Finalize();

  JVForm->AddDomainIntegrator(new VectorFEDivergenceIntegrator);
  JVForm->Assemble();
  JVForm->Finalize();

  //Set the block bilinear/matrix operator
  Array<int> empty_tdof_list;  // empty

  darcyEMOp = new BlockOperator(block_trueOffsets);
  M = JJForm->ParallelAssemble();
  B = JVForm->ParallelAssemble();
  (*B) *= -1;
  Bt = new TransposeOperator(B);

  darcyEMOp->SetBlock(0,0, M);
  darcyEMOp->SetBlock(0,1, Bt);
  darcyEMOp->SetBlock(1,0, B);
};


void DarcyEMProblem::BuildPreconditioner()
{
   //Construct the a Schurr Complement
   //Gauss-Seidel block Preconditioner
   Md = new HypreParVector(MPI_COMM_WORLD, M->GetGlobalNumRows(),M->GetRowStarts());
   M->GetDiag(*Md);

   MinvBt = B->Transpose();
   MinvBt->InvScaleRows(*Md);
   S = ParMult(B, MinvBt);

   invM = new HypreDiagScale(*M);
   invS = new HypreBoomerAMG(*S);

   invM->iterative_mode = false;
   invS->iterative_mode = false;

   //Set the block diagonal bilinear/matrix preconditioning operator
   darcyEMPr = new BlockDiagonalPreconditioner(block_trueOffsets);
   darcyEMPr->SetDiagonalBlock(0, invM);
   darcyEMPr->SetDiagonalBlock(1, invS);
}


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


void DarcyEMProblem::Solve(bool verbosity){

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
};


DarcyEMProblem::~DarcyEMProblem(){
   if(JForm     != NULL) delete JForm;
 /*  if(VForm     != NULL) delete VForm;
   if(darcyEMOp != NULL) delete darcyEMOp;
   if(darcyEMPr != NULL) delete darcyEMPr;
   if(JJForm    != NULL) delete JJForm;
   if(JVForm    != NULL) delete JVForm;
   if(M         != NULL) delete M;
   if(B         != NULL) delete B;
   if(Bt        != NULL) delete Bt;
   if(MinvBt    != NULL) delete MinvBt;
   if(Md        != NULL) delete Md;
   if(S         != NULL) delete S;
   if(invM      != NULL) delete invM;
   if(invS      != NULL) delete invS;
   if(fespaceRT != NULL) delete fespaceRT;
   if(fespaceL  != NULL) delete fespaceL;*/
};
#endif