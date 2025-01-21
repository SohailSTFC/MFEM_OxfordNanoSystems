#ifndef EMPROBLEM2_HPP
#define EMPROBLEM2_HPP

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
// that uses a H(div)-H(Curl)-L2-H1 formulation
// of the EM-problem to generate the
// residual and Jacobian
// because H(Curl) elements are being
// used this Operator is restricted to
// 3-D
//
// The equation system is:
//  J + sigma*Grad(Vb) = 0
//  Div(J) = 0
//  V - Vb = 0
//
class EMOperator : public Operator{
  private:
    int Dim;
    bool PA=true;

	//FE spaces and Forms
    vector<ParFiniteElementSpace*> feSpaces;
    vector<ParLinearForm*>  LForms;  //Applied forces
    vector<int>             L_ints;  //Integrated L-Forms

    ParBilinearForm      *JJForm=NULL,*VVForm=NULL;  //Bilinear type forms
    ParMixedBilinearForm *JVbForm=NULL;              //Mixed components
    ParMixedBilinearForm *VbJForm=NULL,*VVbForm=NULL;//Mixed components

    //Matrix Operators 
    OperatorPtr OpJJ, OpVV;
    OperatorPtr OpJVb, OpVbJ, OpVVb;
    TransposeOperator *OpJVb_T, *OpVbJ_T;

    BlockOperator *EMSolverOp=NULL;

    Array<int> BlockOffsets;
    Array<Array<int> *> ess_bdr;

    mutable BlockVector w_vec, z_vec;
    mutable BlockVector b_vec, x_DirchRef;

    //Electromagnetic Coefficients
    ConstantCoefficient mu, sigma, one;

    //Assemble the linear forms
    void AssembleLinearForms();

  public:
    EMOperator(vector<ParFiniteElementSpace*> feSpaces_
               , Array<int> BOffsets
               , Array<Array<int>*> ess_bdrs
               , int dim
               , MemoryType mt
			   , int OpSize);

    ~EMOperator();

    void SetDirchRefVector(const Vector x_Ref);

    virtual void Mult(const Vector &k, Vector &y) const;

    virtual Operator &GetGradient(const Vector &xp) const{};
};

//
// The class constructor
//
EMOperator::EMOperator(vector<ParFiniteElementSpace*> feSpaces_
                     , Array<int> BOffsets
					 , Array<Array<int>*> ess_bdrs
		             , int dim
                     , MemoryType mt
		             , int OpSize): Operator(OpSize),
		             mu(1.0), sigma(1.0), one(1.0)
{	
  Dim = dim;
  feSpaces.clear();
  LForms.clear();
  L_ints.clear();
  for(int I=0; I<feSpaces_.size(); I++) feSpaces.push_back(new ParFiniteElementSpace(*feSpaces_[I]) );
  for(int I=0; I<feSpaces.size();  I++) LForms.push_back (new ParLinearForm(feSpaces[I]));
  for(int I=0; I<feSpaces.size();  I++) L_ints.push_back(0);

  //Make the Bilinear forms from the FESpaces
  JJForm = new ParBilinearForm(feSpaces[0]);
  VVForm = new ParBilinearForm(feSpaces[2]);

  JVbForm = new ParMixedBilinearForm(feSpaces[0],feSpaces[1]);
  VbJForm = new ParMixedBilinearForm(feSpaces[1],feSpaces[0]);
  VVbForm = new ParMixedBilinearForm(feSpaces[1],feSpaces[2]);

  BlockOffsets = Array<int>(BOffsets);
  ess_bdr = Array<Array<int>*>(ess_bdrs);

  //Set the Vectors (zero them)
  z_vec.Update(BlockOffsets);      z_vec = 0.0;
  w_vec.Update(BlockOffsets);      w_vec = 0.0;
  b_vec.Update(BlockOffsets);      b_vec = 0.0;
  x_DirchRef.Update(BlockOffsets); x_DirchRef = 0.0;

  //Linear Form assembly
  AssembleLinearForms();

  //Bilinear Forms
  Array<int> empty_tdofs;

  JJForm->AddDomainIntegrator( new VectorFEMassIntegrator(one) );
  JJForm->Assemble();
  JJForm->FormSystemMatrix(empty_tdofs, OpJJ);

  VVForm->AddDomainIntegrator( new MassIntegrator(one) );
  VVForm->Assemble();
  VVForm->FormSystemMatrix(empty_tdofs, OpVV);

  //Mixed Bilinear Forms
  JVbForm->AddDomainIntegrator( new VectorFEDivergenceIntegrator);
  JVbForm->Assemble();
  JVbForm->FormRectangularSystemMatrix(empty_tdofs, empty_tdofs, OpJVb);

  VbJForm->AddDomainIntegrator( new MixedScalarWeakGradientIntegrator);
  VbJForm->Assemble();
  VbJForm->FormRectangularSystemMatrix(empty_tdofs, empty_tdofs, OpVbJ);

  VVbForm->AddDomainIntegrator( new MixedScalarMassIntegrator(one));
  VVbForm->Assemble();
  VVbForm->FormRectangularSystemMatrix(empty_tdofs, empty_tdofs, OpVVb);

  //Transpose certain operators
  OpJVb_T = new TransposeOperator(OpJVb.Ptr());
  OpVbJ_T = new TransposeOperator(OpVbJ.Ptr());

  //Block Operator/Matrix
  EMSolverOp = new BlockOperator(BlockOffsets);

  //Row 0
  EMSolverOp->SetBlock(0,0, OpJJ.Ptr(), 1.0);
  EMSolverOp->SetBlock(0,1, OpJVb_T,    1.0);

  //Row 1
  EMSolverOp->SetBlock(1,0, OpVbJ_T,    1.0);

  //Row 2
  EMSolverOp->SetBlock(2,2, OpVV.Ptr(),  1.0);
  EMSolverOp->SetBlock(2,1, OpVVb.Ptr(),-1.0);
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
// Set the reference dirchelet BC vector
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
      if(ess_bdr[I]->Size() != 0){
        Array<int> ess_tDofs;
        feSpaces[I]->GetEssentialTrueDofs(*ess_bdr[I], ess_tDofs);
        applyDirchValues(x_DirchRef.GetBlock(I), w_vec.GetBlock(I), ess_tDofs);
      }
    }
  }
  //Calculate the unconstrained residual
  EMSolverOp->Mult(w_vec,z_vec);

  //Eliminate the residuals corresponding to the constrained boundary
  //DOFs
  if(ess_bdr.Size() != 0){
    for(int I=0; I<feSpaces.size(); I++){
      if(ess_bdr[I]->Size() != 0){
        Array<int> ess_tDofs;
        feSpaces[I]->GetEssentialTrueDofs(*ess_bdr[I], ess_tDofs);
        applyDirchValues(x_DirchRef.GetBlock(I), z_vec.GetBlock(I), ess_tDofs);
      }
    }
  }
  setValues(z_vec, y);
};
#endif

