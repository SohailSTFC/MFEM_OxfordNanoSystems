#ifndef EQUATIONSYSTEM_HPP
#define EQUATIONSYSTEM_HPP 

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <functional>
#include <array>
#include <map>

using namespace std;
using namespace mfem;

#include "BoundaryAndInitialSolution.hpp"
#include "boundaryConditions.hpp"
#include "referenceTypes.hpp"

//
//
//  The problem class/equation system
//
//

//The facade class is the interface to construct MFEM equation system
template<class FacadeClass>
class equationSystem
{
  protected:
    //Problem dimension
    int DIM;

    //Finite element spaces
    I_map<ParFiniteElementSpace*> feSpaces;

    //Block Matrix structure offsets
    Array<int> block_offsets;     // Offsets for the Block system
    Array<int> block_trueOffsets; // True Offsets for the Block system

    //The Bilinear forms of the block components
    // (Jacobian) (Assuming a symmetric Saddle point problem)
    I_map<ParBilinearForm*>        BLForms;  //Bilinear forms (used for unconstrained Precons)
    IJ_map<ParMixedBilinearForm*>  MBLForms; //Mixed Bilinear Forms

    //The linear forms of the block components
    // (residual)
    I_map<ParLinearForm*> LForms;

    //The complete block vectors
    mutable BlockVector x_vec, b_vec;
    mutable BlockVector tx_vec, tb_vec;
    vector<ParGridFunction*> RefFields;

    // Form block operators (operates Matrix multiplication)
	// (This aggregates the block components of the forms)
    BlockOperator               *ProblemOp = NULL;
    BlockDiagonalPreconditioner *ProblemPr = NULL;

    //The Block hypre matrices and Transposes for Jacobian
    IJ_map<HypreParMatrix*>    HypreMats; //Hypre Matrices for Hypre Precons
    IJ_map<OperatorPtr*>       OpPtrs;    //Standard Operators
    IJ_map<TransposeOperator*> TrOps;     //Transpose Operators (0 - No Tr, 1 - Tr)
    IJ_map<int>                TrFlag;    //Transpose operator flag

    //Shared pointer to the solver
    IterativeSolver* solver=NULL;

    //Essential Boundary Conditions
    I_map<Array<int>*>                essBCtag_markers;
    I_map<Array<int>*>                essBCtag_tdofs;
    I_map<Array<Coefficient*>*>       ess_bdr_sCoeffs; //Scalar Coeff
    I_map<Array<VectorCoefficient*>*> ess_bdr_vCoeffs; //Vector Coeff
    vector<int>                       sv_flag; //0=scalar 1=Vector
//I want to implement my own general arbitrary rank/dim 
//differentiable Tensor Coefficients however It is a
//project for later, so this dirty implementation 
//will have todo (multiple arrays of coeffs yuck)


    //Read in and set the Boundary conditions
    void SetBCsArrays();

    //Read in and set the Boundary conditions
    void SetFieldBCs();

  public:
    //The fields (Needed for postprocessing)
    //Made public for external access (not a great solution
    //but hey it works)
    vector<std::string>      FieldNames;
    vector<ParGridFunction*> Fields;

    //The constructor
    equationSystem(FacadeClass FC, MemoryType deviceMT, int dim);

    //Build and set a Preconditioner for the solver
    void BuildPreconditioner();

    //Set a linear/non-linear solver
    void Set_Solver(bool verbose);

    //Solve the equation
    void Solve(bool verbose);

    //Make the fields for post-processing
    void SetFields();

	//The destructor
    ~equationSystem();
};


//
//
//  Implementation of the problem class
//
//

//The constructor
//Constructs the problem and sets-up
//residual+jacobian operators/Forms
template<class FacadeClass>
equationSystem<FacadeClass>::equationSystem(FacadeClass *FC, MemoryType deviceMT, int dim)
{
  DIM = dim;
  FC->SetFESpaces(&feSpaces, &FieldNames, &sv_flag);
  FC->SetBCs(&essBCtag_markers, &ess_bdr_sCoeffs, &ess_bdr_vCoeffs);
  FC->SetBlockBilinearForms(&MBLForms,&TrFlag);

  // Get the block offsets and true block offsets
  // to construct the block structured vectors
  // and matrix operators
  int feSpaceSize = feSpaces.size();
  Array<int> bofs(feSpaceSize+1); // number of variables + 1
  bofs[0] = 0;
  I_map<ParFiniteElementSpace*> 
  for(int I=0; I<feSpaceSize; I++) bofs[I+1] = feSpaces[I]->GetVSize();
  bofs.PartialSum();

  Array<int> btofs(feSpaceSize+1); // number of variables + 1
  btofs[0] = 0;
  for(int I=0; I<feSpaceSize; I++) btofs[I+1] = feSpaces[I]->TrueVSize();
  btofs.PartialSum();
 
  block_offsets     = Array<int>(bofs);
  block_trueOffsets = Array<int>(btofs);


  // Define the solution vectors
  // for a given problem
  x_vec.Update  (block_offsets, deviceMT);
  b_vec.Update  (block_offsets, deviceMT);
  tx_vec.Update (block_trueOffsets, deviceMT);
  tb_vec.Update (block_trueOffsets, deviceMT);
  tb_vec = 0.0;   b_vec  = 0.0;
  tx_vec = 0.0;   x_vec  = 0.0;

  //Setting the boundary conditions
  SetBCsArrays();
  Array<int> ess_tdof_empty;

  //Construct the constrained Operators
  //for the block structured system
  for(auto const& x : MBLForms){
    int I = (x.first)[0];
    int J = (x.first)[1];
    Array<int> essBC_tdofs = *essBCtag_tdofs[I];
    OpPtrs[{I,J}] = new OperatorPtr;
    if(TrFlag[{I,J}]==0) (x.second)->FormRectangularSystemMatrix(ess_tdof_empty, essBC_tdofs, OpPtrs[{I,J}]);
    if(TrFlag[{I,J}]==1) (x.second)->FormRectangularSystemMatrix(essBC_tdofs, ess_tdof_empty, OpPtrs[{I,J}]);
    if(TrFlag[{I,J}]==1) TrOps[{I,J}] = new TransposeOperator(OpPtrs[{I,J}]);
  }

  //Set the block matrix operator
  ProblemOp = new BlockOperator(block_trueOffsets);
  for(auto const& x : OpPtrs){
    int I = (x.first)[0];
    int J = (x.first)[1];
    if(TrFlag[{I,J}]==0) ProblemOp->SetBlock(I,J, (x.second)->Ptr(),1.0);
    if(TrFlag[{I,J}]==1) ProblemOp->SetBlock(I,J, TrOps[{I,J}],1.0);
  }
};



//Sets the natural and essential boundary
//conditions
template<class FacadeClass>
void equationSystem<FacadeClass>::SetBCsArrays(){
  //Build the reference field Grid functions
  //set the Dirch BCs via Coefficients etc..
  for(int I=0; I <feSpaces.size(); I++ ){
    //Construct a grid function
    RefFields.push_back(new ParGridFunction(feSpaces[I]));

    //Essential boundary condition tags
    int nTags = feSpaces[I]->GetMesh()->bdr_attributes.Max();
    if((nTags > 0 )and(essBCtag_markers[I]->Size() > 0)){
      //Get the TrueDof of the boundary tags
      Array<int> ess_bc_tdofs;
      ess_bc_tags = Array<int>(*essBCtag_markers[I]);
      feSpaces[I]->GetEssentialTrueDofs(*essBCtag_markers[I], ess_bc_tdofs);
      essBCtag_tdofs[I] = new Array<int>(ess_bc_tdofs);

      //Construct the reference configuration
      for(int J=0; J<nTags; J++){
        int K = (*essBCtag_markers[I])[J];
        if(K != 0){
          Array<int> bc_tags_tmp(nTags);
		  bc_tags_tmp = 0;
          bc_tags_tmp[J] = 1;

          //I hate doing this but we need both s and v Coeffs for now
          //then we zero them if they are undefined otherwise
          //they should be defined in an array
          Vector ZeroV=0.0;
          Coefficient sCoeff;
          VectorCoefficient vCoeff;

          sCoeff = ((ess_bdr_sCoeffs[I] != ess_bdr_sCoeffs.end())and(ess_bdr_sCoeffs[I]->Size() != 0))
		           ? Coefficient(*(*ess_bdr_sCoeffs[I])[J]) : Coefficient(0.00);

          vCoeff = ((ess_bdr_vCoeffs[I] != ess_bdr_vCoeffs.end())and(ess_bdr_vCoeffs[I]->Size() != 0))
		           ? VectorCoefficient(*(*ess_bdr_vCoeffs[I])[J]) : VectorCoefficient(ZeroV);
	
	      //Checks whether Field is Scalar or Vector
          if(sv_flag[J]==0) RefFields[I]->ProjectBdrCoefficient(sCoeff, bc_tags_tmp);
          if(sv_flag[J]==1) RefFields[I]->ProjectBdrCoefficient(vCoeff, bc_tags_tmp);
        }
      }
    }
  }
};

//Set the BCs by setting the solution vectors and
template<class FacadeClass>
void equationSystem<FacadeClass>:::SetFieldBCs(){
  for(int I=0; I<RefFields.size(); I++){
    if(essBCtag_tdofs[I] != essBCtag_tdofs.end()){
      applyDirchValues(*RefFields[I], (x_vec.GetBlock(I)),  *essBCtag_tdofs[I]);
      applyDirchValues(*RefFields[I], (b_vec.GetBlock(I)),  *essBCtag_tdofs[I]);
      applyDirchValues(*RefFields[I], (tx_vec.GetBlock(I)), *essBCtag_tdofs[I]);
      applyDirchValues(*RefFields[I], (tb_vec.GetBlock(I)), *essBCtag_tdofs[I]);
    }
  }
};

//Builds a preconditioner needed to
//accelerate the Darcy problem solver (Optional)
template<class FacadeClass>
void equationSystem<FacadeClass>::BuildPreconditioner()
{
/*
  //Construct the a Schurr Complement
  //Gauss-Seidel block Preconditioner
  matUU = static_cast<HypreParMatrix*>( OpUU1.Ptr());
  matUP = static_cast<HypreParMatrix*>( OpUP.Ptr() );
  matPU = static_cast<HypreParMatrix*>( OpPU.Ptr() );
  (*matUP) *= -1.0; 
  (*matPU) *= -1.0;
  Md = new HypreParVector(MPI_COMM_WORLD, matUU->GetGlobalNumRows(),matUU->GetRowStarts());
  matUU->GetDiag(*Md);
  MinvBt = matPU->Transpose();
  MinvBt->InvScaleRows(*Md);
  matS1  = ParMult(matUP,MinvBt);

//  invM  = new HypreDiagScale(*matUU);
  invM1 = new HypreBoomerAMG(*matUU);
  invS1 = new HypreBoomerAMG(*matS1);

  invM1->SetInterpolation(6);
  invM1->SetCoarsening(8);
  invM1->SetRelaxType(6);
  invM1->SetCycleNumSweeps(1,1);
  invM1->SetCycleType(2);

  invS1->SetInterpolation(6);
  invS1->SetCoarsening(8);
  invS1->SetRelaxType(6);
  invS1->SetCycleNumSweeps(1,1);
  invS1->SetCycleType(2);

  StokesUPPr = new BlockDiagonalPreconditioner(block_trueOffsets);
  StokesUPPr->SetDiagonalBlock(0,invM1);
  StokesUPPr->SetDiagonalBlock(1,invS1);
*/
}

//Sets the linear/non-linear solver
//for the Darcy problem
template<class FacadeClass>
void equationSystem<FacadeClass>::Set_Solver(bool verbosity){
  int maxIter(1500);
  double rtol(1.e-7);
  double atol(1.e-10);
  solver = new MINRESSolver(MPI_COMM_WORLD);
  solver->SetAbsTol(atol);
  solver->SetRelTol(rtol);
  solver->SetMaxIter(maxIter);
  solver->SetPrintLevel(verbosity);
  if(ProblemOp != NULL) solver->SetOperator(*ProblemOp);
  if(ProblemPr != NULL) solver->SetPreconditioner(*ProblemPr);
};

 //Solves the system of equations
//for the Darcy problem
template<class FacadeClass>
void equationSystem<FacadeClass>::Solve(bool verbosity){
  if(ProblemOp != NULL){
    StopWatch chrono;
    chrono.Clear();
    chrono.Start();
    SetFieldBCs();
    solver->Mult(tb_vec, tx_vec);
    chrono.Stop();
    if (verbosity)
    {
      std::cout << "MINRES ended in "                     << solver->GetNumIterations()
                << " iterations with a residual norm of " << solver->GetFinalNorm() << ".\n";
      std::cout << "MINRES solver took "                  << chrono.RealTime()      << "s. \n";
    }
  }else{
    if (verbosity) std::cout << "Error Stokes operator not built" << ".\n";
  }
};

//Setting-up/unpacking the fields for the problem
//these are needed before post-processing
template<class FacadeClass>
void equationSystem<FacadeClass>::SetFields(){
  //Construct the Field Grid functions
  for(int I=0; I <feSpaces.size(); I++ ){
    //Construct a grid function 
    Fields.push_back(new ParGridFunction);
    Fields[I]->MakeRef(feSpaces[I], x_vec.GetBlock(I), 0);
    Fields[I]->Distribute(&(tx_vec.GetBlock(I)));
  }
};


template<class FacadeClass>
equationSystem<FacadeClass>::~equationSystem(){
   if(ProblemOp != NULL) delete ProblemOp;
   if(ProblemPr != NULL) delete ProblemPr;
};
#endif
