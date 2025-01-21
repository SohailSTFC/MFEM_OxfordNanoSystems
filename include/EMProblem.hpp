#ifndef EMPROBLEM_HPP
#define EMPROBLEM_HPP

#include "mfem.hpp"
#include "boundaryConditions.hpp"
#include <iostream>


using namespace mfem;
using namespace std;
//
// Author: Sohail Rathore
// Date  : 14/01/2025
//
// A Electromagnetic Operator that
// that uses the H1-H(Curl) formulation
// of the EM-problem to generate the
// residual and Jacobian
// because H(Curl) elements are being
// used this Operator is restricted to
// 3-D
//
class EMOperator : public Operator{
  private:
    int Dim;
    bool PA=true;

	//FE spaces and Forms
    vector<ParFiniteElementSpace*> feSpaces;
    vector<ParLinearForm*>  LForms;  //Applied forces
    vector<int>             L_ints;  //Integrated L-Forms

    ParBilinearForm      *BBForm=NULL;             //Linear type forms
    ParMixedBilinearForm *BUForm=NULL,*UBForm=NULL;//Mixed components

    //Matrix Operators 
    OperatorPtr OpBB, OpBU, OpUB;
    BlockOperator *EMSolverOp=NULL;

    Array<int> BlockOffsets;
    Array<Array<int> *> ess_bdr;

    mutable BlockVector w_vec, z_vec;
    mutable BlockVector b_vec, x_DirchRef;

    //Electromagnetic Coefficients
    ConstantCoefficient muInv, sigma, one;
    VectorGridFunctionCoefficient u_vel;

    //Assemble the linear forms
    void AssembleLinearForms();

  public:
    EMOperator(vector<ParFiniteElementSpace*> feSpaces_
               , Array<int> BOffsets
			   , int dim
               , MemoryType mt
			   , int OpSize);

    ~EMOperator();

    void SetBCs(Array<Array<int>*> ess_bdr);

    void SetDirchRefVector(const Vector x_Ref);

    virtual void Mult(const Vector &k, Vector &y) const;

    virtual Operator &GetGradient(const Vector &xp) const;
};

//
// The class constructor
//
EMOperator::EMOperator(vector<ParFiniteElementSpace*> feSpaces_
                     , Array<int> BOffsets
		             , ParGridFunction U
		             , int dim
                     , MemoryType mt
		             , int OpSize): Operator(OpSize),
		             muInv(1.0), rho(1.0), one(1.0), u_vel(&U)
{
  Dim = dim;
  feSpaces.clear();
  LForms.clear();
  L_ints.clear();
  for(int I=0; I<feSpaces_.size(); I++) feSpaces.push_back(new ParFiniteElementSpace(*feSpaces_[I]) );
  for(int I=0; I<feSpaces.size();  I++) LForms.push_back (new ParLinearForm(feSpaces[I]));
  for(int I=0; I<feSpaces.size();  I++) L_ints.push_back(0);

  BBForm = new ParBilinearForm(feSpaces[0]);
  BUForm = new ParMixedBilinearForm(feSpaces[1],feSpaces[0]);
  UBForm = new ParMixedBilinearForm(feSpaces[0],feSpaces[1]);

  BlockOffsets = Array<int>(BOffsets);

  //Set the Vectors (zero them)
  z_vec.Update(BlockOffsets);      z_vec = 0.0;
  w_vec.Update(BlockOffsets);      w_vec = 0.0;
  b_vec.Update(BlockOffsets);      b_vec = 0.0;
  x_DirchRef.Update(BlockOffsets); x_DirchRef = 0.0;

  //Linear Form assembly
  AssembleLinearForms();

  //Bilinear Forms
  Array<int> empty_tdofs;

  BBForm->AddDomainIntegrator( new VectorFECurlIntegrator(muInv) );
  BBForm->AddDomainIntegrator( new MixedCrossProductIntegrator(u_vel) );
  BBForm->Assemble();
  BBForm->FormSystemMatrix(empty_tdofs, OpBB);

  BUForm->AddDomainIntegrator( new VectorDiffusionIntegrator(one) );
  BUForm->Assemble();
  BUForm->FormRectangularSystemMatrix( empty_tdofs, empty_tdofs, OpBU);

  UBForm->AddDomainIntegrator( new GradientIntegrator(one));
  UBForm->Assemble();
  UBForm->FormRectangularSystemMatrix( empty_tdofs, empty_tdofs, OpUB);

  //Block Operator/Matrix
  EMSolverOp = new BlockOperator(BlockOffsets);
  EMSolverOp->SetBlock(0,0, OpBB.Ptr(), 1.0);
  EMSolverOp->SetBlock(0,1, OpBU.Ptr(), 1.0);
  EMSolverOp->SetBlock(1,0, OpUB.Ptr(),-1.0);
}


//
// The class destructor
//
EMOperator::~EMOperator(){
  for(int I=0; I<LForms.size();   I++) delete LForms[I];
  for(int I=0; I<feSpaces.size(); I++) delete feSpaces[I];
  LForms.clear();
  feSpaces.clear();
};

//
// Assemble the linear forms
//
void EMOperator::AssembleLinearForms(){
  for(int I=0; I<LForms.size(); I++){
    if(L_ints[I] != 0){
      LForms[I]->Assemble();
      LForms[I]->SyncAliasMemory(b_vec);
      LForms[I]->ParallelAssemble(b_vec.GetBlock(I));
      b_vec.GetBlock(I).SyncAliasMemory(b_vec);
    }
  }
};

//
// Set the essential BCs
//
void EMOperator::SetBCs(Array<Array<int>*> ess_bdr_){ess_bdr = Array<Array<int>*>(ess_bdr_);};

//
// Set the reference dirchelet BCs
//
void EMOperator::SetDirchRefVector(const Vector x_Ref){setValues(x_Ref,x_DirchRef);};


//
// The residual calculation
//
void EMOperator::Mult(const Vector &k, Vector &y) const
{
  //Set the boundary values to non-Homogenous Dirch values
  setValues(k, w_vec);
  if(ess_bdr.Size() != 0){
    for(int I=0; I<feSpaces.size(); I++){
      Array<int> ess_tDofs;
      feSpaces[I]->GetEssentialTrueDofs(*ess_bdr[I], ess_tDofs);
      applyDirchValues(x_DirchRef.GetBlock(I), w_vec.GetBlock(I), ess_tDofs);
    }
  }
  //Calculate the unconstrained residual
  SNSSolverOp->Mult(w_vec,z_vec);
  z_vec.Add(-1.0,b_vec);

  //Eliminate the residuals corresponding to the constrained boundary
  //DOFs
  if(ess_bdr.Size() != 0){
    for(int I=0; I<feSpaces.size(); I++){
      Array<int> ess_tDofs;
      feSpaces[I]->GetEssentialTrueDofs(*ess_bdr[I], ess_tDofs);
      applyDirchElimination(z_vec.GetBlock(I), ess_tDofs);
    }
  }
  setValues(z_vec, y);
};


//
// Getting the Jacobian matrix
//
Operator &EMOperator::GetGradient(const Vector &xp) const
{
  return NLForm->GetGradient(xp);
};
#endif

