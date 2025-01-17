#ifndef SNSSOLVER_HPP
#define SNSSOLVER_HPP

#include "mfem.hpp"
#include "DirchBCs.hpp"
#include <iostream>


using namespace mfem;
using namespace std;

class SNSSOperator : public Operator{
  private:
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
  public:
    SNSSOperator(vector<FiniteElementSpace*> feSpaces
               , Array<int> BOffsets
               , MemoryType mt);

    void SetNLForm(NonlinearForm *NLForm, int I);

    void SetNLForm(LinearForm *LForm, int I);

    virtual void Mult(const Vector &k, Vector &y) const;

    virtual Operator &GetGradient(const Vector &xp) const;
};



SNSSOperator::SNSSOperator(vector<FiniteElementSpace*> feSpaces
           , Array<int> BOffsets
           , MemoryType mt): mu(1.0), rho(1.0)
{

//   ParBilinearForm *HH_form = new ParBilinearForm(feSpaces[0]);
//   HH_form.AddDomainIntegrator(new VectorDiffusionIntegrator(mu));
}

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

Operator &SNSSOperator::GetGradient(const Vector &xp) const
{
  return NLForm->GetGradient(xp);
};
#endif















