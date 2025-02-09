#ifndef EMPROBLEM_JV_HPP
#define EMPROBLEM_JV_HPP

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
// that uses a H(div)-H1 formulation
// of the EM-problem to generate the
// residual and Jacobian
//
// The equation system is:
//  J + sigma*Grad(V) = 0
//  Div(J) = 0
//

//
// Problem class
//
class EMOperatorJV : public Operator{
  private:
    int Dim;
    bool PA=true;

	//FE spaces and Forms
    vector<ParFiniteElementSpace*> feSpaces;
    vector<ParLinearForm*>  LForms;  //Applied forces
    vector<int>             L_ints;  //Integrated L-Forms

    ParBilinearForm      *JJForm=NULL; //Bilinear type forms
    ParMixedBilinearForm *JVForm=NULL; //Mixed components

    //Matrix Operators 
    OperatorPtr OpJJ, OpJV;
    HypreParMatrix *JJ_mat=NULL, *JV_mat=NULL;
    TransposeOperator *OpJV_T;

    BlockOperator *EMSolverOp=NULL;
    BlockDiagonalPreconditioner *EMPreconOp=NULL;

    Array<int> ess_Jtdof, ess_Vtdof;
    Array<int> BlockOffsets;
    Array<Array<int> *> ess_bdr;

    mutable BlockVector b_vec;

    //Electromagnetic Coefficients
    ConstantCoefficient mu, sigma, one;

    //Assemble the linear forms
    void AssembleLinearForms();

  public:
    EMOperatorJV(vector<ParFiniteElementSpace*> feSpaces_
               , Array<int> BOffsets
               , Array<Array<int>*> ess_bdrs 
               , int dim
               , MemoryType mt
			   , int OpSize);

    ~EMOperatorJV();

    Solver& Precon() const {return *EMPreconOp;};

    virtual void Mult(const Vector &k, Vector &y) const;
};

//
// The class constructor
//
EMOperatorJV::EMOperatorJV(vector<ParFiniteElementSpace*> feSpaces_
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
  JVForm = new ParMixedBilinearForm(feSpaces[0],feSpaces[1]);

  BlockOffsets = Array<int>(BOffsets);
  ess_bdr = Array<Array<int>*>(ess_bdrs);

  //Set the Vectors (zero them)
  b_vec.Update(BlockOffsets); b_vec = 0.0;

  //Linear Form assembly
  AssembleLinearForms();

  //Bilinear Forms
  Array<int> empty_tdofs;
  feSpaces[0]->GetEssentialTrueDofs(*ess_bdr[0], ess_Jtdof);
  feSpaces[1]->GetEssentialTrueDofs(*ess_bdr[1], ess_Vtdof);

  JJForm->AddDomainIntegrator( new VectorFEMassIntegrator(one) );
  JJForm->Assemble();
  JJForm->FormSystemMatrix(ess_Jtdof, OpJJ);

  //Mixed Bilinear Forms
  JVForm->AddDomainIntegrator( new VectorFEDivergenceIntegrator	(sigma));
  JVForm->Assemble();
  JVForm->FormRectangularSystemMatrix(ess_Jtdof,ess_Vtdof, OpJV);

  //Transpose certain operators
  OpJV_T = new TransposeOperator(OpJV.Ptr());

  //Block Operator/Matrix
  EMSolverOp = new BlockOperator(BlockOffsets);

  //Row 0
  EMSolverOp->SetBlock(0,0, OpJJ.Ptr(), 1.0);
  //EMSolverOp->SetBlock(0,1, OpJV.Ptr(),-1.0);
  EMSolverOp->SetBlock(0,1, OpJV_T,-1.0);

  //Row 1
 // EMSolverOp->SetBlock(1,0, OpJV_T, -1.0);
  EMSolverOp->SetBlock(1,0, OpJV.Ptr(), -1.0);
}


//
// The class destructor
//
EMOperatorJV::~EMOperatorJV(){
  for(int I=0; I<LForms.size();   I++) delete LForms[I];
  for(int I=0; I<feSpaces.size(); I++) delete feSpaces[I];
  LForms.clear();
  feSpaces.clear();
};

//
// Assemble the linear forms
//
void EMOperatorJV::AssembleLinearForms(){
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
void EMOperatorJV::Mult(const Vector &k, Vector &y) const
{
  EMSolverOp->Mult(k,y);
};
#endif