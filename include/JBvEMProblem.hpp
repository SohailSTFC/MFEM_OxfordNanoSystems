#ifndef JBvEMProblem_HPP
#define JBvEMProblem_HPP 

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
class JBvEMProblem
{
  protected:
    //Problem dimension
    int DIM;

    //Finite element spaces
    ParFiniteElementSpace *fespaceRT=NULL; //Raviart Thomas elements
    ParFiniteElementSpace *fespaceN=NULL;  //Nedlec elements
    ParFiniteElementSpace *fespaceL=NULL;  //Lagrange finite element space

    //Block Matrix structure offsets
    Array<int> block_offsets;     // number of variables + 1 (2-variables J and v)
    Array<int> block_trueOffsets; // number of variables + 1 (2-variables J and v)

    //The permiativity of space
    double sigma, mu;

    //The Bilinear forms of the block components
    // (Jacobian) (Assuming a symmetric Saddle point problem)
    ParBilinearForm      *JJ_Form=NULL,*BB_Form=NULL; //Bilinear type forms
    ParMixedBilinearForm *JJForm=NULL,*BBForm=NULL;   //Bilinear type forms
    ParMixedBilinearForm *JBForm=NULL,*BJForm=NULL;   //Mixed components
    ParMixedBilinearForm *JVForm=NULL,*VJForm=NULL;   //Mixed components

    //The linear forms of the block components
    // (residual)
    ParLinearForm   *JForm=NULL;
    ParLinearForm   *BForm=NULL;
    ParLinearForm   *VForm=NULL;

    //The complete block vectors
    mutable BlockVector x_vec, b_vec;
    mutable BlockVector tx_vec, tb_vec;

    // Form block operators (operates Matrix multiplication)
	// (This aggregates the block components of the forms)
    BlockOperator               *JBvEMOp = NULL;
    BlockDiagonalPreconditioner *JBvEMPr = NULL;


    //The Block hypre matrices and Transposes for Jacobian
    HypreParMatrix *matJJ=NULL, *matJV=NULL, *matJB=NULL;
    HypreParMatrix *matBB=NULL, *matBJ=NULL;
    HypreParMatrix *matVJ=NULL;
    OperatorPtr OpJJ1, OpBB1;
    OperatorPtr OpJJ, OpBB, OpJB, OpJV, OpBJ, OpVJ;
    TransposeOperator *JVt = NULL,*JBt = NULL, *BJt = NULL;

    //Shared pointer to the solver
    IterativeSolver* solver=NULL;

    //The Preconditioning objects
    HypreParMatrix *MinvBt = NULL, *matS1 = NULL, *matS2 = NULL;
    HypreParVector *Md = NULL, *Md1 = NULL;
    Solver *invM=NULL;
    HypreAMS *invS1=NULL;
    HypreBoomerAMG *invS2=NULL;

    //Boundary Conditions
    vector<vector<double>> DirchVal;               //Dirchelet value of BC
    Array<int> ess_tdof_J, ess_tdof_B, ess_tdof_v; //Dirchelet BC DOF's
    Array<int> ess_bdr_J, ess_bdr_B, ess_bdr_v;    //Dirchelet BC Tag's
  
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
    JBvEMProblem(ParFiniteElementSpace *f1, ParFiniteElementSpace *f2, ParFiniteElementSpace *f3
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
    ~JBvEMProblem();
};


//
//
//  Implementation of the problem class
//
//

//The constructor
//Constructs the problem and sets-up
//residual+jacobian operators/Forms
JBvEMProblem::JBvEMProblem(ParFiniteElementSpace *f1RT
                         , ParFiniteElementSpace *f2N
                         , ParFiniteElementSpace *f3L
                         , double *sig, double *MU, MemoryType deviceMT, int dim)
: sigma(*sig), mu(*MU)
{
  DIM = dim;
  fespaceRT = new ParFiniteElementSpace(*f1RT);
  fespaceN  = new ParFiniteElementSpace(*f2N);
  fespaceL  = new ParFiniteElementSpace(*f3L);

  HYPRE_BigInt dimR = fespaceRT->GlobalTrueVSize();
  HYPRE_BigInt dimN = fespaceN->GlobalTrueVSize();
  HYPRE_BigInt dimW = fespaceL->GlobalTrueVSize();

  std::cout << "***********************************************************\n";
  std::cout << "dim(R) = " << dimR << "\n";
  std::cout << "dim(N) = " << dimN << "\n";
  std::cout << "dim(W) = " << dimW << "\n";
  std::cout << "dim(R+N+W) = " << dimR + dimN + dimW << "\n";
  std::cout << "***********************************************************\n";

  // Get the block offsets and true block offsets
  // to construct the block structured vectors
  // and matrix operators
  Array<int> bofs(4); // number of variables + 1
  bofs[0] = 0;
  bofs[1] = fespaceRT->GetVSize();
  bofs[2] = fespaceN->GetVSize();
  bofs[3] = fespaceL->GetVSize();
  bofs.PartialSum();

  Array<int> btofs(4); // number of variables + 1
  btofs[0] = 0;
  btofs[1] = fespaceRT->TrueVSize();
  btofs[2] = fespaceN->TrueVSize();
  btofs[3] = fespaceL->TrueVSize();
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
  VectorConstantCoefficient U_func(Vector({0.0, 0.0, 0.1}));

  //
  // The Bilinear forms (matrix/jacobian forms)
  //
  //The Bilinear block forms
  JJ_Form  = new ParBilinearForm(fespaceRT); 
  BB_Form  = new ParBilinearForm(fespaceN); 

  JJForm = new ParMixedBilinearForm(fespaceRT, fespaceRT); 
  JBForm = new ParMixedBilinearForm(fespaceRT, fespaceN ); 
  JVForm = new ParMixedBilinearForm(fespaceRT, fespaceL );

  BBForm = new ParMixedBilinearForm(fespaceN , fespaceN);
  BJForm = new ParMixedBilinearForm(fespaceRT, fespaceN);

  VJForm = new ParMixedBilinearForm(fespaceRT, fespaceL);


  //Setting the boundary conditions
  SetBCsArrays();
  Array<int> ess_tdof_empty;

  //Set the integrators/integral forms and assemble the block matrices/bilinear forms
  //Square Bilinear forms (for preconditioning)
  JJ_Form->AddDomainIntegrator(new VectorFEMassIntegrator(One));
  JJ_Form->Assemble();

  BB_Form->AddDomainIntegrator(new VectorFECurlIntegrator(One));
  BB_Form->Assemble();

  //Mixed Bilinear forms (for solution)
  JJForm->AddDomainIntegrator(new MixedVectorMassIntegrator(One));
  JJForm->Assemble();

  JBForm->AddDomainIntegrator(new MixedCrossProductIntegrator(U_func));
  JBForm->Assemble();

  JVForm->AddDomainIntegrator(new MixedVectorWeakDivergenceIntegrator(Sigma));
  JVForm->Assemble();

  BBForm->AddDomainIntegrator(new MixedVectorWeakCurlIntegrator(One));
  BBForm->Assemble();

  BJForm->AddDomainIntegrator(new MixedVectorMassIntegrator(Mu));
  BJForm->Assemble();

  VJForm->AddDomainIntegrator(new MixedVectorWeakDivergenceIntegrator(Sigma));
  VJForm->Assemble();


  //Set the BCs and finalize the bilinear forms (PA)
  JJ_Form->FormSystemMatrix(ess_tdof_empty, OpJJ1);
  BB_Form->FormSystemMatrix(ess_tdof_empty, OpBB1);

  JJForm->FormRectangularSystemMatrix(ess_tdof_empty, ess_tdof_J, OpJJ);
  BBForm->FormRectangularSystemMatrix(ess_tdof_empty, ess_tdof_B, OpBB);

  JBForm->FormRectangularSystemMatrix( ess_tdof_J, ess_tdof_empty, OpJB);
  JVForm->FormRectangularSystemMatrix( ess_tdof_J, ess_tdof_empty, OpJV);
  BJForm->FormRectangularSystemMatrix( ess_tdof_empty, ess_tdof_B, OpBJ);
  VJForm->FormRectangularSystemMatrix( ess_tdof_empty, ess_tdof_v, OpVJ);

  JVt = new TransposeOperator(OpJV.Ptr());
  JBt = new TransposeOperator(OpJB.Ptr());
  
  //Set the block matrix operator
  JBvEMOp = new BlockOperator(block_trueOffsets);

  //Row 0
  JBvEMOp->SetBlock(0,0, OpJJ.Ptr(), 1.0);
//  JBvEMOp->SetBlock(0,1, JBt, 1.0);
  JBvEMOp->SetBlock(0,2, JVt, -1.0);

  //Row 1
  JBvEMOp->SetBlock(1,0, OpBJ.Ptr(), -1.0);
  JBvEMOp->SetBlock(1,1, OpBB.Ptr(),  1.0);

  //Row 2
  JBvEMOp->SetBlock(2,0, OpVJ.Ptr(), -1.0);
};

//Sets the natural and essential boundary
//conditions
void JBvEMProblem::SetBCsArrays(){
  //Essential boundary condition tags
  int nJTags = fespaceRT->GetMesh()->bdr_attributes.Max();
  int nBTags = fespaceN->GetMesh()->bdr_attributes.Max();
  int nVTags = fespaceL->GetMesh()->bdr_attributes.Max();
  
  ess_bdr_J = Array<int>(nJTags);
  ess_bdr_B = Array<int>(nBTags);
  ess_bdr_v = Array<int>(nVTags);

  //initialise the arrays
  ess_bdr_J = 1;
  ess_bdr_B = 0;
  ess_bdr_v = 0;

  //fixed J
  ess_bdr_J[0] = 0;
  ess_bdr_J[1] = 0;
  ess_bdr_J[2] = 0;

  //fixed B
//  ess_bdr_B[0] = 1;
//  ess_bdr_B[1] = 1;
//  ess_bdr_B[2] = 1;

  //fixed v
  ess_bdr_v[0] = 1;
  ess_bdr_v[1] = 1;
  ess_bdr_v[2] = 1;

  //Find the True Dofs
  fespaceRT->GetEssentialTrueDofs(ess_bdr_J, ess_tdof_J);
  fespaceRT->GetEssentialTrueDofs(ess_bdr_B, ess_tdof_B);
  fespaceL->GetEssentialTrueDofs(ess_bdr_v, ess_tdof_v);

  cout << setw(10) << "RT elements: " << setw(10) << ess_tdof_J.Size() << "\n";
  cout << setw(10) << "ND elements: " << setw(10) << ess_tdof_B.Size() << "\n";
  cout << setw(10) << "H1 elements: " << setw(10) << ess_tdof_v.Size() << "\n";

  //The Dirchelet BC values
  vector<double> DirchVal_tmp;
  DirchVal.clear();
  DirchVal_tmp.clear();
  
  //J-Field BC-values
  DirchVal_tmp.push_back(0.00); // N/A
  DirchVal_tmp.push_back(0.00); // N/A
  DirchVal_tmp.push_back(0.00); // N/A
  DirchVal_tmp.push_back(0.00); // n.J = 0
  DirchVal_tmp.push_back(0.00); // N/A
  DirchVal.push_back(DirchVal_tmp);

  //B-Field BC-values
  DirchVal_tmp[0] = 0.00; // Fixed B = c
  DirchVal_tmp[1] = 0.00; // Fixed B = c
  DirchVal_tmp[2] = 0.00; // Fixed B = c
  DirchVal_tmp[3] = 0.00; // N/A
  DirchVal_tmp[4] = 0.00; // N/A
  DirchVal.push_back(DirchVal_tmp);

  //v-Field BC-values
  DirchVal_tmp[0] = 0.00; // Fixed v = c
  DirchVal_tmp[1] = 3.00; // Fixed v = c
  DirchVal_tmp[2] = 3.00; // Fixed v = c
  DirchVal_tmp[3] = 0.00; // N/A
  DirchVal_tmp[4] = 0.00; // N/A
  DirchVal.push_back(DirchVal_tmp);
  DirchVal_tmp.clear();
};

//Set the BCs by setting the solution vectors and
void JBvEMProblem::SetFieldBCs(){
  //Set the boundary Solution function for
  //the current field
  int nJ_tags = fespaceRT->GetMesh()->bdr_attributes.Max();
  int nB_tags = fespaceN->GetMesh()->bdr_attributes.Max();
  int nv_tags = fespaceL->GetMesh()->bdr_attributes.Max();

  cout << setw(10) << "RT element Tags: " << setw(10) << nJ_tags << "\n";
  cout << setw(10) << "H1 element Tags: " << setw(10) << nv_tags << "\n";
  x_vec.GetBlock(2) = 0.2;
  tx_vec.GetBlock(2) = 0.2;

  //Set the J-Field BCs by looping over the
  //active boundaries
  cout << setw(10) << "RT elements: " << setw(10) << ess_tdof_J.Size() << "\n";
  for(int I=0; I<nJ_tags; I++){
    int K = ess_bdr_J[I];
    if(K == 1){ //Checks if boundary is active
      Array<int> ess_tdof, ess_bdr_tmp(nJ_tags);
      ess_bdr_tmp = 0;
      ess_bdr_tmp[I] = 1;
      fespaceRT->GetEssentialTrueDofs(ess_bdr_tmp, ess_tdof);
      cout << setw(10) << I << setw(10) << ess_tdof.Size() << "\n";
      x_vec.GetBlock(0).SetSubVector( ess_tdof, DirchVal[0][I] );
      b_vec.GetBlock(0).SetSubVector( ess_tdof, DirchVal[0][I] );
      tx_vec.GetBlock(0).SetSubVector( ess_tdof, DirchVal[0][I] );
      tb_vec.GetBlock(0).SetSubVector( ess_tdof, DirchVal[0][I] );
    }
  }

  //Set the B-Field BCs by looping over the
  //active boundaries
  cout << setw(10) << "ND elements: " << setw(10) << ess_tdof_B.Size() << "\n";
  for(int I=0; I<nB_tags; I++){
    int K = ess_bdr_B[I];
    if(K == 1){ //Checks if boundary is active
      Array<int> ess_tdof, ess_bdr_tmp(nB_tags);
      ess_bdr_tmp = 0;
      ess_bdr_tmp[I] = 1;
      fespaceN->GetEssentialTrueDofs(ess_bdr_tmp, ess_tdof);
      cout << setw(10) << I << setw(10) << ess_tdof.Size() << "\n";
      x_vec.GetBlock(1).SetSubVector( ess_tdof, DirchVal[1][I] );
      b_vec.GetBlock(1).SetSubVector( ess_tdof, DirchVal[1][I] );
      tx_vec.GetBlock(1).SetSubVector( ess_tdof, DirchVal[1][I] );
      tb_vec.GetBlock(1).SetSubVector( ess_tdof, DirchVal[1][I] );
    }
  }


  //Set the v-Field BCs by looping over the
  //active boundaries
  cout << setw(10) << "H1 elements: " << setw(10) << ess_tdof_v.Size() << "\n";
  for(int I=0; I<nv_tags; I++){
    int K = ess_bdr_v[I];
    if(K == 1){ //Checks if boundary is active
      Array<int> ess_tdof, ess_bdr_tmp(nv_tags);
      ess_bdr_tmp = 0;
      ess_bdr_tmp[I] = 1;
      fespaceL->GetEssentialTrueDofs(ess_bdr_tmp, ess_tdof);
      cout << setw(10) << I << setw(10) << ess_tdof.Size() << "\n";
      x_vec.GetBlock(2).SetSubVector( ess_tdof, DirchVal[2][I] );
      b_vec.GetBlock(2).SetSubVector( ess_tdof, DirchVal[2][I] );
      tx_vec.GetBlock(2).SetSubVector( ess_tdof, DirchVal[2][I] );
      tb_vec.GetBlock(2).SetSubVector( ess_tdof, DirchVal[2][I] );
    }
  }
};


//Builds a preconditioner needed to
//accelerate the Darcy problem solver (Optional)

void JBvEMProblem::BuildPreconditioner()
{
  //Construct the a Schurr Complement
  //Gauss-Seidel block Preconditioner
  matBB = static_cast<HypreParMatrix*>( OpBB1.Ptr() );
  matJJ = static_cast<HypreParMatrix*>( OpJJ1.Ptr() );
  matJV = static_cast<HypreParMatrix*>( OpJV.Ptr() );
  matVJ = static_cast<HypreParMatrix*>( OpVJ.Ptr() );
  (*matJV) *= -1.0; 
  (*matVJ) *= -1.0;
  Md = new HypreParVector(MPI_COMM_WORLD, matJJ->GetGlobalNumRows(),matJJ->GetRowStarts());
  matJJ->GetDiag(*Md);
  MinvBt = matVJ->Transpose();
  MinvBt->InvScaleRows(*Md);
  matS1  = ParMult(matJV, MinvBt);

/*
    OperatorPtr OpJJ1, OpBB1;
    OperatorPtr OpJJ, OpBB, OpJB, OpJV, OpBJ, OpVJ;
    TransposeOperator *JVt = NULL,*JBt = NULL, *BJt = NULL;

    HypreParMatrix *matJJ=NULL, *matJV=NULL, *matJB=NULL;
    HypreParMatrix *matBB=NULL, *matBJ=NULL;
    HypreParMatrix *matVJ=NULL;

    HypreParMatrix *MinvBt = NULL, *matS1 = NULL, *matS2 = NULL;
    HypreParVector *Md = NULL, *Md1 = NULL;
    Solver *invM=NULL;
    Solver *invS1=NULL;
    HypreBoomerAMG *invS2=NULL;

    invS1 = new OperatorJacobiSmoother(Md_PA, ess_tdof_list);
*/
  invM  = new HypreADS(*matJJ, fespaceRT);
  invS1 = new HypreAMS(*matBB, fespaceN);
  invS2 = new HypreBoomerAMG(*matS1);

  invS2->SetInterpolation(6);
  invS2->SetCoarsening(8);
  invS2->SetRelaxType(6);
  invS2->SetCycleNumSweeps(1,1);
  invS2->SetCycleType(2);

  JBvEMPr = new BlockDiagonalPreconditioner(block_trueOffsets);
  JBvEMPr->SetDiagonalBlock(0,invM);
  JBvEMPr->SetDiagonalBlock(1,invS1);
  JBvEMPr->SetDiagonalBlock(2,invS2);
}

//Sets the linear/non-linear solver
//for the Darcy problem
void JBvEMProblem::Set_Solver(bool verbosity){
  int maxIter(75);
  double rtol(1.e-7);
  double atol(1.e-10);
  solver = new MINRESSolver(MPI_COMM_WORLD);
  solver->SetAbsTol(atol);
  solver->SetRelTol(rtol);
  solver->SetMaxIter(maxIter);
  solver->SetPrintLevel(verbosity);
  if(JBvEMOp != NULL) solver->SetOperator(*JBvEMOp);
  if(JBvEMPr != NULL) solver->SetPreconditioner(*JBvEMPr);
};

 //Solves the system of equations
//for the Darcy problem
void JBvEMProblem::Solve(bool verbosity){
  if(JBvEMOp != NULL){
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
void JBvEMProblem::SetFields(){
  FieldNames.push_back("Current");
  Fields.push_back(new ParGridFunction);
  Fields[0]->MakeRef(fespaceRT, x_vec.GetBlock(0), 0);
  Fields[0]->Distribute(&(tx_vec.GetBlock(0)));

  FieldNames.push_back("B-Field");
  Fields.push_back(new ParGridFunction);
  Fields[1]->MakeRef(fespaceN, x_vec.GetBlock(1), 0);
  Fields[1]->Distribute(&(tx_vec.GetBlock(1)));

  FieldNames.push_back("Potential");
  Fields.push_back(new ParGridFunction);
  Fields[2]->MakeRef(fespaceL, x_vec.GetBlock(2), 0);
  Fields[2]->Distribute(&(tx_vec.GetBlock(2)));
};


JBvEMProblem::~JBvEMProblem(){
   if(JForm   != NULL) delete JForm;
   if(VForm   != NULL) delete VForm;
   if(JBvEMOp != NULL) delete JBvEMOp;
   if(JBvEMPr != NULL) delete JBvEMPr;
   if(JJ_Form != NULL) delete JJ_Form;
   if(JJForm  != NULL) delete JJForm;
   if(JVForm  != NULL) delete JVForm;
   if(VJForm  != NULL) delete VJForm;
};
#endif
