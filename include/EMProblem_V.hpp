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
class EMOperatorV : public Operator{
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

    Array<int> ess_Vtdof;
    Array<int> BlockOffsets;
    Array<Array<int> *> ess_bdr;

    mutable BlockVector b_vec;

    //Electromagnetic Coefficients
    ConstantCoefficient mu, sigma, one;

    //Assemble the linear forms
    void AssembleLinearForms();

  public:
    EMOperatorV(vector<ParFiniteElementSpace*> feSpaces_
               , Array<Array<int>*> ess_bdrs
               , int dim
               , MemoryType mt
			   , int OpSize);

    ~EMOperatorV();

    void ReassembleSystem();

    virtual void Mult(const Vector &k, Vector &y) const;
};

//
// The class constructor
//
EMOperatorV::EMOperatorV(vector<ParFiniteElementSpace*> feSpaces_
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
  ess_bdr = Array<Array<int>*>(ess_bdrs);

  //Bilinear Forms
  VVForm->AddDomainIntegrator( new DiffusionIntegrator(one) );

  //Assemble the linear system
  ReassembleSystem();
}


//
// The class destructor
//
EMOperatorV::~EMOperatorV(){
  for(int I=0; I<LForms.size();   I++) delete LForms[I];
  for(int I=0; I<feSpaces.size(); I++) delete feSpaces[I];
  LForms.clear();
  feSpaces.clear();
};

//
// Assemble the linear forms
//
void EMOperatorV::AssembleLinearForms(){
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
// Reassembles all the operators
// used for AMR
//
void EMOperatorV::ReassembleSystem(){
  Array<int> BOfSets_tmp(feSpaces.size()+1);
  BOfSets_tmp[0] = 0;
  for(int I=0; I<feSpaces.size(); I++) BOfSets_tmp[I+1] = feSpaces[I]->GetTrueVSize();
  BOfSets_tmp.PartialSum();
  BlockOffsets = Array<int>(BOfSets_tmp);

  feSpaces[0]->GetEssentialTrueDofs(*ess_bdr[0], ess_Vtdof);

  //Set the Vectors (zero them)
  b_vec.Update(BlockOffsets); b_vec = 0.0;

  //Linear Form assembly
  AssembleLinearForms();

  //Bilinear Forms
  Array<int> empty_tdofs;
  VVForm->Update();
  VVForm->Assemble();
  VVForm->FormSystemMatrix(ess_Vtdof, OpVV);

  //Block Operator/Matrix
  if(EMSolverOp != NULL) delete EMSolverOp;
  EMSolverOp = new BlockOperator(BlockOffsets);

  //Row 0
  EMSolverOp->SetBlock(0,0, OpVV.Ptr(), 1.0);
};


//
// The residual calculation
//
void EMOperatorV::Mult(const Vector &k, Vector &y) const
{
  //Calculate the unconstrained residual
  EMSolverOp->Mult(k,y);
};
#endif

