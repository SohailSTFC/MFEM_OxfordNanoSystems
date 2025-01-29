#ifndef EMPROBLEM_V_HPP
#define EMPROBLEM_V_HPP

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
// that uses a H1 formulation
// of the EM-problem to generate the
// residual and Jacobian
//
// The equation system is:
//  Div( sigma*Grad(Vb) ) = F
//
class EMOperator : public Operator{
  private:
    int Dim;
    bool PA=true;

	//FE spaces and Forms
    vector<ParFiniteElementSpace*> feSpaces;
    vector<ParLinearForm*>  LForms;  //Applied forces
    vector<int>             L_ints;  //Integrated L-Forms

    ParBilinearForm *VVForm=NULL;   //Linear type forms

    //Matrix Operators 
    OperatorPtr OpVV;

    BlockOperator *EMSolverOp=NULL;

    Array<int> ess_Jtdof, ess_Btdof, ess_Vtdof;
    Array<int> BlockOffsets;
    Array<Array<int> *> ess_bdr;

    mutable BlockVector b_vec;

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

    virtual void Mult(const Vector &k, Vector &y) const;
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
  VVForm = new ParBilinearForm(feSpaces[0]);

  BlockOffsets = Array<int>(BOffsets);
  ess_bdr = Array<Array<int>*>(ess_bdrs);

  feSpaces[0]->GetEssentialTrueDofs(*ess_bdr[0], ess_Jtdof);

  //Set the Vectors (zero them)
  b_vec.Update(BlockOffsets);      b_vec = 0.0;

  //Linear Form assembly
  AssembleLinearForms();

  //Bilinear Forms
  Array<int> empty_tdofs;

  VVForm->AddDomainIntegrator( new DiffusionIntegrator(one) );
  VVForm->Assemble();
  VVForm->FormSystemMatrix(ess_Vtdof, OpVV);

  //Block Operator/Matrix
  EMSolverOp = new BlockOperator(BlockOffsets);

  //Row 0
  EMSolverOp->SetBlock(0,0, OpVV.Ptr(), 1.0);
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
// The residual calculation
//
void EMOperator::Mult(const Vector &k, Vector &y) const
{
  //Calculate the unconstrained residual
  EMSolverOp->Mult(k,y);
};
#endif

