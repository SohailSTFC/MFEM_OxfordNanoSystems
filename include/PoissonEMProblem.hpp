#ifndef poissonEMprOBLEM_HPP
#define poissonEMprOBLEM_HPP 

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

#include "BoundaryAndInitialSolution.hpp"
//
//
//  The problem class
//
//
class poissonEMproblem
{
  protected:
    //Problem dimension
    int DIM;

    //Finite element spaces
    ParFiniteElementSpace *fespaceH1=NULL;  //Lagrange finite element space

    //Block Matrix structure offsets
    Array<int> block_offsets;     // number of variables + 1 (2-variables J and v)
    Array<int> block_trueOffsets; // number of variables + 1 (2-variables J and v)

    //The permiativity of space
    double sigma;

    //The Bilinear forms of the block components
    ParMixedBilinearForm *VV_Form=NULL;

    //The linear forms of the block components
    // (residual)
    ParLinearForm   *VForm=NULL;

    //The complete block vectors
    mutable BlockVector x_vec, b_vec;
    mutable BlockVector tx_vec, tb_vec;

    // Form block operators (operates Matrix multiplication)
    // (This aggregates the block components of the forms)
    BlockOperator               *poissonEMOp = NULL;
    BlockDiagonalPreconditioner *poissonEMpr = NULL;

    //The Block hypre matrices and Transposes for Jacobian
    TransposeOperator *Bt = NULL;
    HypreParMatrix *matM=NULL;
    OperatorPtr opM;

    //Shared pointer to the solver
    IterativeSolver* solver=NULL;

    //The Preconditioning objects
    HypreBoomerAMG *invM=NULL;

    //Boundary Conditions
    vector<vector<double>> DirchVal;  //Dirchelet value of BC
    Array<int>  ess_tdof_v;           //Dirchelet BC DOF's
    Array<int>  ess_bdr_v;            //Dirchelet BC Tag's
  
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
    poissonEMproblem(ParFiniteElementSpace *f1, double sig, MemoryType deviceMT, int dim);

    //Read in and set the Boundary conditions
    void UpdateArrayBCs(const int & BCsdTags, const double & BCdVals
	                  , const Array<int> & BCsTags, const Vector & BCVals);

    //Build the problem operator
    void BuildOperator();

    //Build and set a Preconditioner for the solver
    void BuildPreconditioner();

    //Set a linear/non-linear solver
    void Set_Solver(bool verbose);

    //Solve the equation
    void Solve(bool verbose);

    //Make the fields for post-processing
    void SetFields();

	//The destructor
    ~poissonEMproblem();
};


//
//
//  Implementation of the problem class
//
//

//The constructor
//Constructs the problem and sets-up
//residual+jacobian operators/Forms
poissonEMproblem::poissonEMproblem(ParFiniteElementSpace *fL
                             , double sig, MemoryType deviceMT, int dim)
: sigma(sig)
{
  DIM = dim;
  fespaceH1  = new ParFiniteElementSpace(*fL);
  HYPRE_BigInt dimW = fespaceH1->GlobalTrueVSize();

  std::cout << "***********************************************************\n";
  std::cout << "dim(W) = " << dimW << "\n";
  std::cout << "***********************************************************\n";

  // Get the block offsets and true block offsets
  // to construct the block structured vectors
  // and matrix operators
  Array<int> bofs(2); // number of variables + 1
  bofs[0] = 0;
  bofs[1] = fespaceH1->GetVSize();
  bofs.PartialSum();

  Array<int> btofs(2); // number of variables + 1
  btofs[0] = 0;
  btofs[1] = fespaceH1->TrueVSize();
  btofs.PartialSum();
 
  block_offsets     = Array<int>(bofs);
  block_trueOffsets = Array<int>(btofs);


  //Setting the boundary conditions
  SetBCsArrays();

  // 10. Define the parallel grid function and parallel linear forms, solution
  //     vector and rhs.
  x_vec.Update  (block_offsets, deviceMT);
  b_vec.Update  (block_offsets, deviceMT);
  tx_vec.Update (block_trueOffsets, deviceMT);
  tb_vec.Update (block_trueOffsets, deviceMT);
  tb_vec = 0.0;
  b_vec  = 0.0;
  tx_vec = 0.0;
  x_vec  = 0.0;

  // 9. Define the coefficients, analytical solution, and rhs of the PDE.
  // the coefficients and functions
  ConstantCoefficient k(sigma);

  //
  // The Bilinear forms (matrix/jacobian forms)
  //
  Array<int> ess_tdof_empty;
  VV_Form = new ParMixedBilinearForm(fespaceH1, fespaceH1);
  VV_Form->AddDomainIntegrator(new DiffusionIntegrator(k));
  VV_Form->Assemble();
};


//Build the problem operator
void poissonEMproblem::BuildOperator(){
  Array<int> ess_tdof_empty;
  VV_Form->FormRectangularSystemMatrix(ess_tdof_empty, ess_tdof_v, opM);

  //Set the block matrix operator
  poissonEMOp = new BlockOperator(block_trueOffsets);
  poissonEMOp->SetBlock(0,0, opM.Ptr());
};

//Sets the default boundary
//conditions
void poissonEMproblem::SetBCsArrays(){
  //Boundary condition tags
  int nTagsMax = fespaceH1->GetMesh()->bdr_attributes.Max();
  ess_bdr_v = Array<int>(nTagsMax);

  //initialise the arrays
  ess_bdr_v = 0;

  //fixed v
  if(nTagsMax > 3){
    ess_bdr_v[0] = 1;
    ess_bdr_v[1] = 1;
    ess_bdr_v[2] = 1;
  }

  //Find the True Dofs
  fespaceH1->GetEssentialTrueDofs(ess_bdr_v, ess_tdof_v);

  //The Dirchelet BC values
  vector<double> DirchVal_tmp;
  DirchVal.clear();
  DirchVal_tmp.clear();
  for(int I=0; I<nTagsMax; I++) DirchVal_tmp.push_back( 0.00);
  
  //v-Field BC-values
  if(nTagsMax > 3){
    DirchVal_tmp[1] = 3.00; // Fixed v = c
    DirchVal_tmp[2] = 3.00; // Fixed v = c
  }
  DirchVal.push_back(DirchVal_tmp);
  DirchVal_tmp.clear();
};


//Read in and set the Boundary conditions
void poissonEMproblem::UpdateArrayBCs(const int & BCsdTags, const double & BCdVals
	                              , const Array<int>& BCsTags, const Vector & BCVals)
{
  //Number of boundary tags
  int nTagsMax = fespaceH1->GetMesh()->bdr_attributes.Max();

  //Initialise the arrays with the default tag
  ess_bdr_v = BCsdTags;

  //Set the Tag value that opposes the default
  int OppTagValV = ((BCsdTags == 1) ? 0 : 1);

  //Set the tags
  if((BCsTags.Size() > 0)and(BCsTags != NULL)){//Makes sure array if not empty/NULL
    for(int J=0; J<BCsTags.Size(); J++){
      int K = BCsTags[J];
      if((K > nTagsMax)and(K < 0)) continue;
      ess_bdr_v[K] = OppTagValV;
    }
  }

  //Find the True Dofs
  fespaceH1->GetEssentialTrueDofs(ess_bdr_v, ess_tdof_v);

  //v-Field BC-values Defaults
  for(int I=0; I<nTagsMax; I++) DirchVal[0][I] = BCdVals;

  //v-Field BC-values Unique
  if(BCVals.Size() == BCsTags.Size()){//Check whether the number of Array values correspond
    if((BCsTags.Size() > 0)and(BCsTags != NULL)){//Makes sure array if not empty/NULL
      for(int J=0; J<BCsTags.Size(); J++){
        int K = BCsTags[J];
        if((K > nTagsMax)or(K < 0)) continue;
        DirchVal[0][K] =  BCVals[J];
      }
    }
  }
  //End of function
};


//Set the BCs by setting the solution vectors and
void poissonEMproblem::SetFieldBCs(){
  //Set the boundary Solution function for
  //the current field
  int nv_tags = fespaceH1->GetMesh()->bdr_attributes.Max();
  cout << setw(10) << "H1 element Tags: " << setw(10) << nv_tags << "\n";
  x_vec.GetBlock(0)  = 0.2;
  tx_vec.GetBlock(0) = 0.2;

  //Set the V-Field BCs by looping over the
  //active boundaries
  cout << setw(10) << "H1 elements: " << setw(10) << ess_tdof_v.Size() << "\n";
  for(int I=0; I<nv_tags; I++){
    int K = ess_bdr_v[I];
    if(K == 1){ //Checks if boundary is active
      Array<int> ess_tdof, ess_bdr_tmp(nv_tags);
      ess_bdr_tmp = 0;
      ess_bdr_tmp[I] = 1;
      fespaceH1->GetEssentialTrueDofs(ess_bdr_tmp, ess_tdof);
      cout << setw(10) << I << setw(10) << ess_tdof.Size() << "\n";
      x_vec.GetBlock(0).SetSubVector(  ess_tdof, DirchVal[0][I] );
      b_vec.GetBlock(0).SetSubVector(  ess_tdof, DirchVal[0][I] );
      tx_vec.GetBlock(0).SetSubVector( ess_tdof, DirchVal[0][I] );
      tb_vec.GetBlock(0).SetSubVector( ess_tdof, DirchVal[0][I] );
    }
  }
};


//Builds a preconditioner needed to
//accelerate the Darcy problem solver (Optional)

void poissonEMproblem::BuildPreconditioner()
{
  //Construct the a Schurr Complement
  //Gauss-Seidel block Preconditioner
  matM = static_cast<HypreParMatrix*>( opM.Ptr() );
  invM = new HypreBoomerAMG(*matM);
  invM->SetInterpolation(6);
  invM->SetCoarsening(8);
  invM->SetRelaxType(6);
  invM->SetCycleNumSweeps(2,2);
  invM->SetCycleType(2);

  poissonEMpr = new BlockDiagonalPreconditioner(block_trueOffsets);
  poissonEMpr->SetDiagonalBlock(0,invM);
}

//Sets the linear/non-linear solver
//for the Darcy problem
void poissonEMproblem::Set_Solver(bool verbosity){
  int maxIter(2500);
  double rtol(1.e-9);
  double atol(1.e-12);
  solver = new MINRESSolver(MPI_COMM_WORLD);
  solver->SetAbsTol(atol);
  solver->SetRelTol(rtol);
  solver->SetMaxIter(maxIter);
  solver->SetPrintLevel(verbosity);
  if(poissonEMOp != NULL) solver->SetOperator(*poissonEMOp);
  if(poissonEMpr != NULL) solver->SetPreconditioner(*poissonEMpr);
};

 //Solves the system of equations
//for the Darcy problem
void poissonEMproblem::Solve(bool verbosity){
  if(poissonEMOp != NULL){
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
    if (verbosity) std::cout << "Error Darcy operator not built" << ".\n";
  }
};

//Setting-up/unpacking the fields for the Darcy problem
//these are needed before post-processing
void poissonEMproblem::SetFields(){
  FieldNames.push_back("Potential");
  Fields.push_back(new ParGridFunction);
  Fields[0]->MakeRef(fespaceH1, x_vec.GetBlock(0), 0);
  Fields[0]->Distribute(&(tx_vec.GetBlock(0)) );
};


poissonEMproblem::~poissonEMproblem(){
   if(VForm       != NULL) delete VForm;
   if(poissonEMOp != NULL) delete poissonEMOp;
   if(poissonEMpr != NULL) delete poissonEMpr;
   if(VV_Form     != NULL) delete VV_Form;
};
#endif
