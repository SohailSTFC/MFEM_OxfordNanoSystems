#ifndef StokesProblem_HPP
#define StokesProblem_HPP 

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

#include "BoundaryAndInitialSolution.hpp"
#include "boundaryConditions.hpp"

//
// 1-D Parabolic Inlet function
//
void InletFunc(const Vector & x, Vector & f){
   const double Vmax = 0.001;
   const double L    = 0.500;
   const double W    = 0.500;
   const double D    = 0.000;
   const double xO   = 0.000;
   const double yO   = 0.000;
   const double zO   = 0.000;

   const double xC = (x(0) - xO)/L;
   const double yC = (x.Size() > 1) ? ((x(1) - yO)/W) : 0.00;
   const double zC = (x.Size() > 2) ? ((x(2) - zO)/D) : 0.00;
   f = 0.0;
   f(0) = 4.0*yC*(1.0 - yC)*Vmax;
  // if(x.Size() > 1) f(1) = 4.0*yC*(1.0 - yC)*Vmax;
  // if(x.Size() > 2) f(2) = 0.00;
};

//
//
//  The problem class
//
//
class StokesUPProblem
{
  protected:
    //Problem dimension
    int DIM;

    //Finite element spaces
    ParFiniteElementSpace *fespaceH1_U=NULL; //H1 elements of Order p+1
    ParFiniteElementSpace *fespaceH1_P=NULL; //H1 elements of Order p

    //Block Matrix structure offsets
    Array<int> block_offsets;     // number of variables + 1 (2-variables J and v)
    Array<int> block_trueOffsets; // number of variables + 1 (2-variables J and v)

    //The permiativity of space
    double sigma, mu;

    //The Bilinear forms of the block components
    // (Jacobian) (Assuming a symmetric Saddle point problem)
    ParBilinearForm      *UUForm=NULL;
    ParMixedBilinearForm *UPForm=NULL, *PUForm=NULL; //Mixed components

    //The linear forms of the block components
    // (residual)
    ParLinearForm   *UForm=NULL;
    ParLinearForm   *PForm=NULL;

    //The complete block vectors
    mutable BlockVector x_vec, b_vec;
    mutable BlockVector tx_vec, tb_vec;

    // Form block operators (operates Matrix multiplication)
	// (This aggregates the block components of the forms)
    BlockOperator               *StokesUPOp = NULL;
    BlockDiagonalPreconditioner *StokesUPPr = NULL;


    //The Block hypre matrices and Transposes for Jacobian
    HypreParMatrix *matUU=NULL, *matUP=NULL, *matPU=NULL;
    OperatorPtr OpUU, OpUP, OpPU;
    TransposeOperator *UPt = NULL;

    //Shared pointer to the solver
    IterativeSolver* solver=NULL;

    //The Preconditioning objects
    HypreParMatrix *MinvBt = NULL, *matS1 = NULL;
    HypreParVector *Md = NULL;
    HypreBoomerAMG *invS1=NULL, *invM1=NULL;

    //Boundary Conditions
    vector<double> DirchVal; //Dirchelet value of BC
    Array<int> ess_tdof_U, ess_tdof_P;  //Dirchelet BC DOF's
    Array<int> ess_bdr_U, ess_bdr_P;    //Dirchelet BC Tag's
  
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
    StokesUPProblem(ParFiniteElementSpace *f1, ParFiniteElementSpace *f2
	              , double *sig, double *MU, MemoryType deviceMT, int dim);

    //Build and set a Preconditioner for the solver
    void BuildPreconditioner();

    //Set a linear/non-linear solver
    void Set_Solver(bool verbose);

    //Solve the equation
    void Solve(bool verbose);

    //Make the fields for post-processing
    void SetFields();

	//The destructor
    ~StokesUPProblem();
};


//
//
//  Implementation of the problem class
//
//

//The constructor
//Constructs the problem and sets-up
//residual+jacobian operators/Forms
StokesUPProblem::StokesUPProblem(ParFiniteElementSpace *f1RT, ParFiniteElementSpace *f2L
                               , double *sig, double *MU, MemoryType deviceMT, int dim)
: sigma(*sig), mu(*MU)
{
  DIM = dim;
  fespaceH1_U = new ParFiniteElementSpace(*f1RT);
  fespaceH1_P  = new ParFiniteElementSpace(*f2L);

  HYPRE_BigInt dimR = fespaceH1_U->GlobalTrueVSize();
  HYPRE_BigInt dimW = fespaceH1_P->GlobalTrueVSize();

  std::cout << "***********************************************************\n";
  std::cout << "dim(R)   = " << dimR << "\n";
  std::cout << "dim(W)   = " << dimW << "\n";
  std::cout << "dim(R+W) = " << dimR + dimW << "\n";
  std::cout << "***********************************************************\n";

  // Get the block offsets and true block offsets
  // to construct the block structured vectors
  // and matrix operators
  Array<int> bofs(3); // number of variables + 1
  bofs[0] = 0;
  bofs[1] = fespaceH1_U->GetVSize();
  bofs[2] = fespaceH1_P->GetVSize();
  bofs.PartialSum();

  Array<int> btofs(3); // number of variables + 1
  btofs[0] = 0;
  btofs[1] = fespaceH1_U->TrueVSize();
  btofs[2] = fespaceH1_P->TrueVSize();
  btofs.PartialSum();
 
  block_offsets     = Array<int>(bofs);
  block_trueOffsets = Array<int>(btofs);


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
  ConstantCoefficient One(1.00);
  ConstantCoefficient Sigma(sigma);
  ConstantCoefficient Mu(mu);

  //
  // The Bilinear forms (matrix/jacobian forms)
  //
  //The Bilinear block forms
  UUForm = new ParBilinearForm(fespaceH1_U); 
  UPForm = new ParMixedBilinearForm(fespaceH1_U, fespaceH1_P); 
  PUForm = new ParMixedBilinearForm(fespaceH1_U, fespaceH1_P);

  //Setting the boundary conditions
  SetBCsArrays();
  Array<int> ess_tdof_empty;

  //Set the integrators/integral forms and assemble the block matrices/bilinear forms
  //Mixed Bilinear forms (for solution)
  UUForm->AddDomainIntegrator(new VectorDiffusionIntegrator(One));
  UUForm->Assemble();

  UPForm->AddDomainIntegrator(new VectorDivergenceIntegrator(One));
  UPForm->Assemble();

  PUForm->AddDomainIntegrator(new VectorDivergenceIntegrator(One));
  PUForm->Assemble();

  //Set the BCs and finalize the bilinear forms (PA)
  UUForm->FormSystemMatrix(ess_tdof_U, OpUU);
  UPForm->FormRectangularSystemMatrix(ess_tdof_U,ess_tdof_empty, OpUP);
  PUForm->FormRectangularSystemMatrix(ess_tdof_empty, ess_tdof_P, OpPU);
  UPt = new TransposeOperator(OpUP.Ptr());
  
  //Set the block matrix operator
  StokesUPOp = new BlockOperator(block_trueOffsets);

  //Row 0
  StokesUPOp->SetBlock(0,0, OpUU.Ptr(), 1.0);
  StokesUPOp->SetBlock(0,1, UPt, -1.0);

  //Row 1
  StokesUPOp->SetBlock(1,0, OpPU.Ptr(), -1.0);
};

//Sets the natural and essential boundary
//conditions
void StokesUPProblem::SetBCsArrays(){
  //Essential boundary condition tags
  int nUTags = fespaceH1_U->GetMesh()->bdr_attributes.Max();
  int nPTags = fespaceH1_P->GetMesh()->bdr_attributes.Max();
  ess_bdr_U = Array<int>(nUTags);
  ess_bdr_P = Array<int>(nPTags);

  //initialise the arrays
  ess_bdr_U = 1;

  //not fixed U
  ess_bdr_U[7] = 0;
  ess_bdr_U[5] = 0;
/*
  ess_bdr_U[0] = 0;
  ess_bdr_U[1] = 0;
  ess_bdr_U[2] = 0;
*/
  ess_bdr_P = 0;
  ess_bdr_P[4] = 1;
  fespaceH1_P->GetEssentialTrueDofs(ess_bdr_P, ess_tdof_P);

  //Find the True Dofs
  fespaceH1_U->GetEssentialTrueDofs(ess_bdr_U, ess_tdof_U);
  cout << setw(10) << "RT elements: " << setw(10) << ess_tdof_U.Size() << "\n";
};
/*
Physical Surface(7) = {9};
Physical Surface(6) = {8};
Physical Surface(5) = {7};
Physical Surface(4) = {6};
*/

//Set the BCs by setting the solution vectors and
void StokesUPProblem::SetFieldBCs(){
  //Set the boundary Solution function for
  //the current field
  int nP_tags = fespaceH1_P->GetMesh()->bdr_attributes.Max();
  int nU_tags = fespaceH1_U->GetMesh()->bdr_attributes.Max();
  cout << setw(10) << "RT element Tags: " << setw(10) << nU_tags << "\n";
  x_vec.GetBlock(0)  = 0.00;
  tx_vec.GetBlock(0) = 0.00;
  x_vec.GetBlock(1)  = 0.00;
  tx_vec.GetBlock(1) = 0.00;


  //Applied inlet condition
  Array<int> ess_tdof0, ess_bdr_tmp0(nU_tags);
  ess_bdr_tmp0 = 0;
  ess_bdr_tmp0[4] = 1;
  ParGridFunction *velocity  = new ParGridFunction(fespaceH1_U);


  VectorFunctionCoefficient vel(DIM, InletFunc);
  fespaceH1_U->GetEssentialTrueDofs(ess_bdr_tmp0, ess_tdof0);
  velocity->ProjectBdrCoefficient(vel, ess_bdr_tmp0);
  applyDirchValues(*velocity, (x_vec.GetBlock(0)),  ess_tdof0);
  applyDirchValues(*velocity, (b_vec.GetBlock(0)),  ess_tdof0);
  applyDirchValues(*velocity, (tx_vec.GetBlock(0)), ess_tdof0);
  applyDirchValues(*velocity, (tb_vec.GetBlock(0)), ess_tdof0);
  delete velocity;


  //Applied inlet pressure
  x_vec.GetBlock(1).SetSubVector( ess_tdof_P, 0.00 );
  b_vec.GetBlock(1).SetSubVector( ess_tdof_P, 0.00 );
  tx_vec.GetBlock(1).SetSubVector(ess_tdof_P, 0.00 );
  tb_vec.GetBlock(1).SetSubVector(ess_tdof_P, 0.00 );



  //Set the J-Field BCs by looping over the
  //active boundaries
  cout << setw(10) << "RT elements: " << setw(10) << ess_tdof_U.Size() << "\n";
  for(int I=0; I<nU_tags; I++){
    int K = ess_bdr_U[I];
    if((K == 1)and(I != 4)){ //Checks if boundary is active
      Array<int> ess_tdof, ess_bdr_tmp(nU_tags);
      ess_bdr_tmp = 0;
      ess_bdr_tmp[I] = 1;
      fespaceH1_U->GetEssentialTrueDofs(ess_bdr_tmp, ess_tdof);
      cout << setw(10) << I << setw(10) << ess_tdof.Size() << "\n";
      x_vec.GetBlock(0).SetSubVector( ess_tdof, 0.00 );
      b_vec.GetBlock(0).SetSubVector( ess_tdof, 0.00 );
      tx_vec.GetBlock(0).SetSubVector(ess_tdof, 0.00 );
      tb_vec.GetBlock(0).SetSubVector(ess_tdof, 0.00 );
    }
  }
};


//Builds a preconditioner needed to
//accelerate the Darcy problem solver (Optional)

void StokesUPProblem::BuildPreconditioner()
{
  //Construct the a Schurr Complement
  //Gauss-Seidel block Preconditioner
  matUU = static_cast<HypreParMatrix*>( OpUU.Ptr());
  matUP = static_cast<HypreParMatrix*>( OpUP.Ptr() );
  matPU = static_cast<HypreParMatrix*>( OpPU.Ptr() );
  (*matUP) *= -1.0; 
  (*matPU) *= -1.0;
  Md = new HypreParVector(MPI_COMM_WORLD, matUU->GetGlobalNumRows(),matUU->GetRowStarts());
  matUU->GetDiag(*Md);
  MinvBt = matPU->Transpose();
  MinvBt->InvScaleRows(*Md);
  matS1  = ParMult(matUP,MinvBt);

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
}

//Sets the linear/non-linear solver
//for the Darcy problem
void StokesUPProblem::Set_Solver(bool verbosity){
  int maxIter(1500);
  double rtol(1.e-7);
  double atol(1.e-10);
  solver = new MINRESSolver(MPI_COMM_WORLD);
  solver->SetAbsTol(atol);
  solver->SetRelTol(rtol);
  solver->SetMaxIter(maxIter);
  solver->SetPrintLevel(verbosity);
  if(StokesUPOp != NULL) solver->SetOperator(*StokesUPOp);
  if(StokesUPPr != NULL) solver->SetPreconditioner(*StokesUPPr);
};

 //Solves the system of equations
//for the Darcy problem
void StokesUPProblem::Solve(bool verbosity){
  if(StokesUPOp != NULL){
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

//Setting-up/unpacking the fields for the Darcy problem
//these are needed before post-processing
void StokesUPProblem::SetFields(){
  FieldNames.push_back("Velocity");
  Fields.push_back(new ParGridFunction);
  Fields[0]->MakeRef(fespaceH1_U, x_vec.GetBlock(0), 0);
  Fields[0]->Distribute(&(tx_vec.GetBlock(0)));

  FieldNames.push_back("Pressure");
  Fields.push_back(new ParGridFunction);
  Fields[1]->MakeRef(fespaceH1_P, x_vec.GetBlock(1), 0);
  Fields[1]->Distribute(&(tx_vec.GetBlock(1)));
};


StokesUPProblem::~StokesUPProblem(){
  // if(JForm   != NULL) delete JForm;
  // if(VForm   != NULL) delete VForm;
   if(StokesUPOp != NULL) delete StokesUPOp;
   if(StokesUPPr != NULL) delete StokesUPPr;
};
#endif
