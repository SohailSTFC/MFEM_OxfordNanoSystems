#ifndef SNSSOLVER_HPP
#define SNSSOLVER_HPP

#include "mfem.hpp"
#include "boundaryConditions.hpp"
#include <iostream>


using namespace mfem;
using namespace std;

class SNSSOperator : public Operator{
  private:
    int Dim;
    vector<FiniteElementSpace*> feSpaces;
    BlockNonlinearForm        *NLForm;  //Nonlinear stuff
    vector<BilinearForm*>      BLForms; //Linear type forms
    vector<MixedBilinearForm*> MBLForms;//Mixed components
    vector<LinearForm*>        LForms;  //Applied forces

    Array<int> BlockOffsets;
    Array<Array<int> *> ess_bdr;

    mutable BlockVector w_vec, z_vec;
    mutable BlockVector b_vec, x_DirchRef;

    //Fluid coefficients
    ConstantCoefficient mu, rho;

    //Update the 

  public:
    SNSSOperator(vector<FiniteElementSpace*> feSpaces_
               , Array<int> BOffsets
			   , int dim
               , MemoryType mt);

    ~SNSSOperator();

    virtual void Mult(const Vector &k, Vector &y) const;

    virtual Operator &GetGradient(const Vector &xp) const;
};

//
// The class constructor
//
SNSSOperator::SNSSOperator(vector<FiniteElementSpace*> feSpaces_
           , Array<int> BOffsets
		   , int dim
           , MemoryType mt):
		   mu(1.0), rho(1.0)
{
  Dim = dim;
  feSpaces.clear();
  LForms.clear();
  for(int I=0; I<feSpaces_.size(); I++) feSpaces.push_back(new FiniteElementSpace(*feSpaces_[I]) );
  for(int I=0; I<feSpaces.size();  I++) LForms.push_back (new LinearForm(feSpaces[I]));
  BlockOffsets = Array<int>(BOffsets);

  //Set the Vectors
  z_vec.Update(BOffsets);      z_vec = 0.0;
  w_vec.Update(BOffsets);      w_vec = 0.0;
  b_vec.Update(BOffsets);      b_vec = 0.0;
  x_DirchRef.Update(BOffsets); x_DirchRef = 0.0;

  //Linear Forms
  Vector RhoG(dim);
  RhoG = 0.0;
  RhoG[dim-1] = 9.81;
  VectorConstantCoefficient Rhog(RhoG);
  LForms[0]->AddDomainIntegrator(new VectorDomainLFIntegrator(Rhog));


//ParBilinearForm *HH_form = new ParBilinearForm(feSpaces[0]);
//HH_form.AddDomainIntegrator(new VectorDiffusionIntegrator(mu));
//vector<*>;  //Applied forces
}

//
// The class destructor
//
SNSSOperator::~SNSSOperator(){
  for(int I=0; I<feSpaces.size(); I++) delete feSpaces[I];
  for(int I=0; I<LForms.size();   I++) delete LForms[I];
  feSpaces.clear();
  LForms.clear();
};

//
// The residual calculation
//
void SNSSOperator::Mult(const Vector &k, Vector &y) const
{
  //Set the bloundary values to non-Homogenous Dirch values
  setValues(k, w_vec);
  for(int I=0; I<feSpaces.size(); I++){
    Array<int> ess_tDofs;
    feSpaces[I]->GetEssentialTrueDofs(*ess_bdr[I], ess_tDofs);
    applyDirchValues(x_DirchRef.GetBlock(I), w_vec.GetBlock(I), ess_tDofs);
  }

  //Calculate the unconstrained residual
  //  z_vec

  //Eliminate the residuals corresponding to the constrained boundary
  //DOFs
  for(int I=0; I<feSpaces.size(); I++){
    Array<int> ess_tDofs;
    feSpaces[I]->GetEssentialTrueDofs(*ess_bdr[I], ess_tDofs);
    applyDirchElimination(z_vec.GetBlock(I), ess_tDofs);
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















