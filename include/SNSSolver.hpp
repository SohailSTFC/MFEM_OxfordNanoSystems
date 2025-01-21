#ifndef SNSSOLVER_HPP
#define SNSSOLVER_HPP

#include "mfem.hpp"
#include "boundaryConditions.hpp"
#include <iostream>


using namespace mfem;
using namespace std;
//
// Author: Sohail Rathore
// Date  : 14/01/2025
//
// A steady state Navier-Stokes
// Operator that calculates the
// residual and Jacobian using
// the Taylor-Hood element
// H1(order)-H1(order-1)
//
//
class SNSSOperator : public Operator{
  private:
    int Dim;
    bool PA=true;

	//FE spaces and Forms
    vector<ParFiniteElementSpace*> feSpaces;
    ParBlockNonlinearForm  *NLForm;  //Nonlinear stuff
    vector<ParLinearForm*>  LForms;  //Applied forces
    vector<int>             L_ints;  //Integrated L-Forms

    ParBilinearForm      *UUForm=NULL, *PPForm=NULL; //Linear type forms
    ParMixedBilinearForm *UPForm=NULL; //Mixed components

    //Matrix Operators 
    OperatorPtr OpUU, OpUP, OpPP;
    BlockOperator *SNSSolverOp=NULL;
//////OperatorHandle Jacobian;

    Array<int> BlockOffsets;
    Array<Array<int> *> ess_bdr;

    mutable BlockVector w_vec, z_vec;
    mutable BlockVector b_vec, x_DirchRef;

    //Fluid coefficients
    ConstantCoefficient mu, rho, one;

    //Assemble the linear forms
    void AssembleLinearForms();

  public:
    SNSSOperator(vector<ParFiniteElementSpace*> feSpaces_
               , Array<int> BOffsets
			   , int dim
               , MemoryType mt
			   , int OpSize);

    ~SNSSOperator();

    void SetBCs(Array<Array<int>*> ess_bdr);

    void SetDirchRefVector(const Vector x_Ref);

    virtual void Mult(const Vector &k, Vector &y) const;

    virtual Operator &GetGradient(const Vector &xp) const;

    virtual void SetOperator(const Operator &op){};
};

//
// The class constructor
//
SNSSOperator::SNSSOperator(vector<ParFiniteElementSpace*> feSpaces_
           , Array<int> BOffsets
		   , int dim
           , MemoryType mt
		   , int OpSize): Operator(OpSize),
		   mu(1.0), rho(1.0), one(1.0)
{
  Dim = dim;
  feSpaces.clear();
  LForms.clear();
  L_ints.clear();
  for(int I=0; I<feSpaces_.size(); I++) feSpaces.push_back(new ParFiniteElementSpace(*feSpaces_[I]) );
  for(int I=0; I<feSpaces.size();  I++) LForms.push_back (new ParLinearForm(feSpaces[I]));
  for(int I=0; I<feSpaces.size();  I++) L_ints.push_back(0);

  UUForm = new ParBilinearForm(feSpaces[0]);
  PPForm = new ParBilinearForm(feSpaces[1]);
  UPForm = new ParMixedBilinearForm(feSpaces[1],feSpaces[0]);

  BlockOffsets = Array<int>(BOffsets);

  //Set the Vectors (zero them)
  z_vec.Update(BlockOffsets);      z_vec = 0.0;
  w_vec.Update(BlockOffsets);      w_vec = 0.0;
  b_vec.Update(BlockOffsets);      b_vec = 0.0;
  x_DirchRef.Update(BlockOffsets); x_DirchRef = 0.0;

  //Linear Form assembly
  Vector RhoG(dim);
  RhoG = 0.0;
  RhoG[dim-1] = 9.81;
  VectorConstantCoefficient Rhog(RhoG);
  LForms[0]->AddDomainIntegrator  ( new VectorDomainLFIntegrator(Rhog) );
  LForms[0]->AddBoundaryIntegrator( new VectorBoundaryLFIntegrator(Rhog) );
  L_ints[0] = 1;
  AssembleLinearForms();

  //Bilinear Forms
  Array<int> empty_tdofs;

  UUForm->AddDomainIntegrator( new VectorDiffusionIntegrator(mu) );
  UUForm->Assemble();
  UUForm->FormSystemMatrix(empty_tdofs, OpUU);

  PPForm->AddDomainIntegrator( new VectorDiffusionIntegrator(one) );
  PPForm->Assemble();
  PPForm->FormSystemMatrix(empty_tdofs, OpPP);

  UPForm->AddDomainIntegrator( new GradientIntegrator(one));
  UPForm->Assemble();
  UPForm->FormRectangularSystemMatrix( empty_tdofs, empty_tdofs, OpUP);

  //Block Operator/Matrix
  SNSSolverOp = new BlockOperator(BlockOffsets);
  SNSSolverOp->SetBlock(0,0, OpUU.Ptr(), 1.0);
  SNSSolverOp->SetBlock(1,1, OpPP.Ptr(), 1.0);
  SNSSolverOp->SetBlock(0,1, OpUP.Ptr(),-1.0);
}


//
// The class destructor
//
SNSSOperator::~SNSSOperator(){
  for(int I=0; I<LForms.size();   I++) delete LForms[I];
  for(int I=0; I<feSpaces.size(); I++) delete feSpaces[I];
  LForms.clear();
  feSpaces.clear();
};


//
// Assemble the linear forms
//
void SNSSOperator::AssembleLinearForms(){
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
void SNSSOperator::SetBCs(Array<Array<int>*> ess_bdr_){ess_bdr = Array<Array<int>*>(ess_bdr_);};

//
// Set the reference dirchelet BCs
//
void SNSSOperator::SetDirchRefVector(const Vector x_Ref){setValues(x_Ref,x_DirchRef);};


//
// The residual calculation
//
void SNSSOperator::Mult(const Vector &k, Vector &y) const
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
Operator &SNSSOperator::GetGradient(const Vector &xp) const
{
  return NLForm->GetGradient(xp);
};
#endif















