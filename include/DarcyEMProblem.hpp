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
    ParFiniteElementSpace *fespaceRT; //Raviart Thomas elements
    ParFiniteElementSpace *fespaceL;  //Lagrange finite element space

    //Block Matrix structure offsets
    Array<int> block_offsets;     // number of variables + 1 (2-variables J and v)
    Array<int> block_trueOffsets; // number of variables + 1 (2-variables J and v)

    //The permiativity of space
    real_t sigma;

    //The linear and Bilinear forms of the block components
    // (residual and jacobian)
    ParBilinearForm JJForm;
    ParLinearForm   *JForm;
    ParLinearForm   *VForm;

    //The complete block vectors
    BlockVector x, rhs;
    BlockVector trueX, trueRhs;


    // Form block operators 
	// (This aggregates the block components of the forms)
    BlockOperator *darcyOp;


 //  BlockOperator *darcyOp = new BlockOperator(block_trueOffsets);
 //  Array<int> empty_tdof_list;  // empty
  // OperatorPtr opM, opB;
  // TransposeOperator *Bt = NULL;


  public:
	//The constructor
    DarcyEMProblem(ParFiniteElementSpace *f1, ParFiniteElementSpace *f2
	             , real_t sig, MemoryType deviceMT, int dim);


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
: fespaceRT(f1RT), fespaceL(f2L), JJForm(fespaceRT), sigma(sig)
{

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
  x       = BlockVector(block_offsets, deviceMT);
  rhs     = BlockVector(block_offsets, deviceMT);
  trueX   = BlockVector(block_trueOffsets, deviceMT);
  trueRhs = BlockVector(block_trueOffsets, deviceMT);

   // 9. Define the coefficients, analytical solution, and rhs of the PDE.
   ConstantCoefficient k(1.0);

   VectorFunctionCoefficient fcoeff(dim, fFun);
   FunctionCoefficient fnatcoeff(f_natural);
   FunctionCoefficient gcoeff(gFun);

   VectorFunctionCoefficient ucoeff(dim, uFun_ex);
   FunctionCoefficient pcoeff(pFun_ex);


   //DarcyEMProblem demoProb(fespaceRT, W_space, sig);
   //v-Linear-form of the equation J + grad v = Je
   VForm = new ParLinearForm;
   VForm->Update(fespaceRT, rhs.GetBlock(0), 0);
   VForm->AddDomainIntegrator(new VectorFEDomainLFIntegrator(fcoeff));
   VForm->AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator(fnatcoeff));
   VForm->Assemble();
   VForm->SyncAliasMemory(rhs);
   VForm->ParallelAssemble(trueRhs.GetBlock(0));
   trueRhs.GetBlock(0).SyncAliasMemory(trueRhs);

   //J-Linear-form of the equation  div J = q
   JForm = new ParLinearForm;
   JForm->Update(fespaceL, rhs.GetBlock(1), 0);
   JForm->AddDomainIntegrator(new DomainLFIntegrator(gcoeff));
   JForm->Assemble();
   JForm->SyncAliasMemory(rhs);
   JForm->ParallelAssemble(trueRhs.GetBlock(1));
   trueRhs.GetBlock(1).SyncAliasMemory(trueRhs);
};

/*
BlockOperator  JacobianOperator()
{
	
}

BlockOperator  PreconditionerOperator()
{
	
}*/


DarcyEMProblem::~DarcyEMProblem(){
   //delete darcyOp;
   delete JForm;
   delete VForm;
};
#endif