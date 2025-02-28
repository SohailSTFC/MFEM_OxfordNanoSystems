#ifndef DARCYEMPROBLEM_HPP
#define DARCYEMPROBLEM_HPP 

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
class DarcyEMProblem
{
  protected:
    //Problem dimension
    int DIM;

    //Finite element spaces
    ParFiniteElementSpace *fespaceRT=NULL; //Raviart Thomas elements
    ParFiniteElementSpace *fespaceH1=NULL;  //Lagrange finite element space

    //Block Matrix structure offsets
    Array<int> block_offsets;     // number of variables + 1 (2-variables J and v)
    Array<int> block_trueOffsets; // number of variables + 1 (2-variables J and v)

    //The permiativity of space
    double sigma;

    //The Bilinear forms of the block components
    // (Jacobian) (Assuming a symmetric Saddle point problem)
    ParBilinearForm      *JJForm=NULL;
    ParMixedBilinearForm *JJ_Form=NULL, *JVForm=NULL, *VJForm=NULL;

    //The linear forms of the block components
    // (residual)
    ParLinearForm   *JForm=NULL;
    ParLinearForm   *VForm=NULL;

    //The complete block vectors
    mutable BlockVector x_vec, b_vec;
    mutable BlockVector tx_vec, tb_vec;

    // Form block operators (operates Matrix multiplication)
	// (This aggregates the block components of the forms)
    BlockOperator               *darcyEMOp = NULL;
    BlockDiagonalPreconditioner *darcyEMPr = NULL;


    //The Block hypre matrices and Transposes for Jacobian
    TransposeOperator *Bt = NULL;
    HypreParMatrix *matM=NULL, *matB=NULL, *matC=NULL;
    OperatorPtr opM1, opM, opB, opC;

    //Shared pointer to the solver
    IterativeSolver* solver=NULL;

    //The Preconditioning objects
    HypreParMatrix *MinvBt = NULL, *matS = NULL;;
    HypreParVector *Md = NULL;
    Solver *invM=NULL;
    HypreBoomerAMG *invS=NULL;

    //Boundary Conditions
    vector<vector<double>> DirchVal;       //Dirchelet value of BC
    Array<int> ess_tdof_J, ess_tdof_v;     //Dirchelet BC DOF's
    Array<int> ess_bdr_J, ess_bdr_v;       //Dirchelet BC Tag's
  
    //Set the default Boundary conditions
    void SetBCsArrays();

    //Set the Boundary conditions before solving
    void SetFieldBCs();
  public:
    //The fields (Needed for postprocessing)
    //Made public for external access (not a great solution
    //but hey it works)
    vector<std::string>      FieldNames;
    vector<ParGridFunction*> Fields;

    //The constructor
    DarcyEMProblem(ParFiniteElementSpace *f1, ParFiniteElementSpace *f2
	             , double sig, MemoryType deviceMT, int dim);

    //Read in and set the Boundary conditions
    void UpdateArrayBCs(const Array<int>& BCsdTags, const Vector & BCdVals
	                  , const Array<Array<int>*>& BCsTags, const Array<Vector*>& BCVals);

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
    ~DarcyEMProblem();
};


//
//
//  Implementation of the problem class
//
//

//The constructor
//Constructs the problem and sets-up
//residual+jacobian operators/Forms
DarcyEMProblem::DarcyEMProblem(ParFiniteElementSpace *f1RT
                             , ParFiniteElementSpace *f2L
                             , double sig, MemoryType deviceMT, int dim)
: sigma(sig)
{
  DIM = dim;
  fespaceRT = new ParFiniteElementSpace(*f1RT);
  fespaceH1  = new ParFiniteElementSpace(*f2L);

  HYPRE_BigInt dimR = fespaceRT->GlobalTrueVSize();
  HYPRE_BigInt dimW = fespaceH1->GlobalTrueVSize();

  std::cout << "***********************************************************\n";
  std::cout << "dim(R) = " << dimR << "\n";
  std::cout << "dim(W) = " << dimW << "\n";
  std::cout << "dim(R+W) = " << dimR + dimW << "\n";
  std::cout << "***********************************************************\n";

  // Get the block offsets and true block offsets
  // to construct the block structured vectors
  // and matrix operators
  Array<int> bofs(3); // number of variables + 1
  bofs[0] = 0;
  bofs[1] = fespaceRT->GetVSize();
  bofs[2] = fespaceH1->GetVSize();
  bofs.PartialSum();

  Array<int> btofs(3); // number of variables + 1
  btofs[0] = 0;
  btofs[1] = fespaceRT->TrueVSize();
  btofs[2] = fespaceH1->TrueVSize();
  btofs.PartialSum();
 
  block_offsets     = *(new Array<int>(bofs));
  block_trueOffsets = *(new Array<int>(btofs));


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
  ConstantCoefficient k(1.0/sigma);

  //
  // The Bilinear forms (matrix/jacobian forms)
  //

  //The Bilinear block forms
  JJForm  = new ParBilinearForm(fespaceRT); 
  JJ_Form = new ParMixedBilinearForm(fespaceRT, fespaceRT); 
  JVForm  = new ParMixedBilinearForm(fespaceRT, fespaceH1);
  VJForm  = new ParMixedBilinearForm(fespaceRT, fespaceH1);

  //Setting the boundary conditions
  SetBCsArrays();

  //Set the integrators/integral forms and assemble the block matrices/bilinear forms
  JJForm->AddDomainIntegrator(new VectorFEMassIntegrator(k));
  JJForm->Assemble();

  JJ_Form->AddDomainIntegrator(new MixedVectorMassIntegrator(k));
  JJ_Form->Assemble();

  JVForm->AddDomainIntegrator(new MixedVectorWeakDivergenceIntegrator);
  JVForm->Assemble();

  VJForm->AddDomainIntegrator(new MixedVectorWeakDivergenceIntegrator);
  VJForm->Assemble();
};


//Build the problem operator
void DarcyEMProblem::BuildOperator(){
  Array<int> ess_tdof_empty;

  //Set the BCs and finalize the bilinear forms (PA)
  JJForm->FormSystemMatrix(ess_tdof_empty, opM1);
  JJ_Form->FormRectangularSystemMatrix(ess_tdof_empty, ess_tdof_J, opM);
  JVForm->FormRectangularSystemMatrix( ess_tdof_J, ess_tdof_empty, opB);
  VJForm->FormRectangularSystemMatrix( ess_tdof_empty, ess_tdof_v, opC);
  Bt = new TransposeOperator(opB.Ptr());

  //Set the block matrix operator
  darcyEMOp = new BlockOperator(block_trueOffsets);
  darcyEMOp->SetBlock(0,0, opM.Ptr());
  darcyEMOp->SetBlock(0,1, Bt, -1.0);
  darcyEMOp->SetBlock(1,0, opC.Ptr(), -1.0);
};


//Set the default Boundary conditions
void DarcyEMProblem::SetBCsArrays(){ //Default case
  int nJTags = fespaceRT->GetMesh()->bdr_attributes.Max();
  int nVTags = fespaceH1->GetMesh()->bdr_attributes.Max();
  int nTagsMax = std::max(nVTags,nJTags);

  //Essential boundary condition tags
  ess_bdr_J = Array<int>(nJTags);
  ess_bdr_v = Array<int>(nVTags);

  //initialise the arrays
  ess_bdr_J = 1; //assume fixed (no-Flux)
  ess_bdr_v = 0; //assume unfixed (Free potential)

  if(nTagsMax > 3){
    //un-fixed J
    ess_bdr_J[0] = 0;
    ess_bdr_J[1] = 0;
    ess_bdr_J[2] = 0;

    //fixed v
    ess_bdr_v[0] = 1;
    ess_bdr_v[1] = 1;
    ess_bdr_v[2] = 1;
  }

  //Find the True Dofs
  fespaceRT->GetEssentialTrueDofs(ess_bdr_J, ess_tdof_J);
  fespaceH1->GetEssentialTrueDofs(ess_bdr_v, ess_tdof_v);

  //The Dirchelet BC values
  vector<double> DirchVal_tmp;
  DirchVal.clear();
  DirchVal_tmp.clear();
  for(int I=0; I<nTagsMax; I++) DirchVal_tmp.push_back(0.00);

  //J-Field BC-values
  for(int I=0; I<nTagsMax; I++) DirchVal_tmp[I] = 0.00;
  DirchVal.push_back(DirchVal_tmp);

  //v-Field BC-values
  for(int I=0; I<nTagsMax; I++) DirchVal_tmp[I] = 0.00;
  if(nTagsMax > 3){
    DirchVal_tmp[1] = 3.00; // Fixed v = c
    DirchVal_tmp[2] = 3.00; // Fixed v = c
  }
  DirchVal.push_back(DirchVal_tmp);

  //clear the tmp arrays
  DirchVal_tmp.clear();
};



//Read in and set the Boundary conditions
void DarcyEMProblem::UpdateArrayBCs(const Array<int>& BCsdTags, const Vector & BCdVals
	                              , const Array<Array<int>*>& BCsTags, const Array<Vector*>& BCVals)
{
  int nJTags = fespaceRT->GetMesh()->bdr_attributes.Max();
  int nVTags = fespaceH1->GetMesh()->bdr_attributes.Max();
  int nTagsMax = std::max(nVTags,nJTags);

  //Initialise the arrays with the default tags
  ess_bdr_J = BCsdTags[0];
  ess_bdr_v = BCsdTags[1];

  //Set the Tag value that opposes the default
  int OppTagValJ = ((BCsdTags[0] == 1) ? 0 : 1);
  int OppTagValV = ((BCsdTags[1] == 1) ? 0 : 1);

  //Set the tags
  if( BCsTags.Size() == 2){//Only do this if the size is correct
    for(int I=0; I<2; I++){
      if((BCsTags[I]->Size() > 0)and(BCsTags[I] != NULL)){//Makes sure array if not empty/NULL
        for(int J=0; J<BCsTags[I]->Size(); J++){
          int K = (*BCsTags[I])[J];
          if((K > nTagsMax)and(K < 0)) continue;
          if(I==0)ess_bdr_J[K] = OppTagValJ;
          if(I==1)ess_bdr_v[K] = OppTagValV;
        }
      }
    }
  }

  //Find the True Dofs
  fespaceRT->GetEssentialTrueDofs(ess_bdr_J, ess_tdof_J);
  fespaceH1->GetEssentialTrueDofs(ess_bdr_v, ess_tdof_v);

  //J-Field and v-Field BC-values Defaults
  for(int J=0; J<2; J++) 
    for(int I=0; I<nTagsMax; I++) 
      DirchVal[J][I] = BCdVals[J];


  //J-Field and v-Field BC-values Unique
  if(BCVals.Size() == BCsTags.Size()){//Check whether Meta Array values correspond
    if( BCsTags.Size() == 2){//Only do this if the size is correct
      for(int I=0; I<2; I++){//Loop over all fields
        if((BCsTags[I]->Size() > 0)and(BCsTags[I] != NULL)){//Makes sure array if not empty/NULL
          for(int J=0; J<BCsTags[I]->Size(); J++){
            int K = (*BCsTags[I])[J];
            if((K > nTagsMax)or(K < 0)) continue;
            DirchVal[I][K] =  (*BCVals[I])[J];
          }
        }
      }
    }
  }
  //End of function
};


//Set the BCs by setting the solution vectors and
void DarcyEMProblem::SetFieldBCs(){
  //Set the boundary Solution function for
  //the current field
  int nJTags = fespaceRT->GetMesh()->bdr_attributes.Max();
  int nVTags = fespaceH1->GetMesh()->bdr_attributes.Max();
  int nTagsMax = max(nVTags,nJTags);

  x_vec.GetBlock(1) = 0.2;
  tx_vec.GetBlock(1) = 0.2;

  //Set the J-Field BCs by looping over the
  //active boundaries
  for(int P=0; P <2; P++){ //Over all field-blocks
    for(int I=0; I<nTagsMax; I++){
      int K=0;
      if((P==0)and(I < nJTags)) K = ess_bdr_J[I];
      if((P==1)and(I < nVTags)) K = ess_bdr_v[I];
      if(K == 1){ //Checks if boundary is active
        Array<int> ess_tdof, ess_bdr_tmp(nTagsMax);
        ess_bdr_tmp = 0;
        ess_bdr_tmp[I] = 1;
        if(P==0) fespaceRT->GetEssentialTrueDofs(ess_bdr_tmp, ess_tdof);
        if(P==1) fespaceH1->GetEssentialTrueDofs(ess_bdr_tmp, ess_tdof);
	  
        x_vec.GetBlock(P).SetSubVector( ess_tdof, DirchVal[P][I] );
        b_vec.GetBlock(P).SetSubVector( ess_tdof, DirchVal[P][I] );
        tx_vec.GetBlock(P).SetSubVector( ess_tdof, DirchVal[P][I] );
        tb_vec.GetBlock(P).SetSubVector( ess_tdof, DirchVal[P][I] );
      }
    }
  }
};


//Builds a preconditioner needed to
//accelerate the Darcy problem solver (Optional)
void DarcyEMProblem::BuildPreconditioner()
{
  if(darcyEMOp == NULL) throw("Error Problem Operator not built");

  //Construct the a Schurr Complement
  //Gauss-Seidel block Preconditioner
  matM = static_cast<HypreParMatrix*>( opM1.Ptr() );
  matB = static_cast<HypreParMatrix*>( opB.Ptr() );
  matC = static_cast<HypreParMatrix*>( opC.Ptr() );
  (*matB) *= -1.0; 
  (*matC) *= -1.0;

  Md = new HypreParVector(MPI_COMM_WORLD, matM->GetGlobalNumRows(),matM->GetRowStarts());
  matM->GetDiag(*Md);
  MinvBt = matB->Transpose();
  MinvBt->InvScaleRows(*Md);
  matS = ParMult(matC, MinvBt);
  if(DIM==3) invM = new HypreADS(*matM, fespaceRT);
  if(DIM!=3) invM = new HypreDiagScale(*matM);
  invS = new HypreBoomerAMG(*matS);

  invS->SetInterpolation(6);
  invS->SetCoarsening(8);
  invS->SetRelaxType(6);
  invS->SetCycleNumSweeps(1,1);
  invS->SetCycleType(2);

  darcyEMPr = new BlockDiagonalPreconditioner(block_trueOffsets);
  darcyEMPr->SetDiagonalBlock(0,invM);
  darcyEMPr->SetDiagonalBlock(1,invS);
}

//Sets the linear/non-linear solver
//for the Darcy problem
void DarcyEMProblem::Set_Solver(bool verbosity){
  int maxIter(2500);
  double rtol(1.e-7);
  double atol(1.e-10);
  solver = new MINRESSolver(MPI_COMM_WORLD);
  solver->SetAbsTol(atol);
  solver->SetRelTol(rtol);
  solver->SetMaxIter(maxIter);
  solver->SetPrintLevel(verbosity);
  if(darcyEMOp != NULL) solver->SetOperator(*darcyEMOp);
  if(darcyEMPr != NULL) solver->SetPreconditioner(*darcyEMPr);
};

 //Solves the system of equations
//for the Darcy problem
void DarcyEMProblem::Solve(bool verbosity){
  if(darcyEMOp != NULL){
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
void DarcyEMProblem::SetFields(){
  FieldNames.push_back("Current");
  Fields.push_back(new ParGridFunction);
  Fields[0]->MakeRef(fespaceRT, x_vec.GetBlock(0), 0);
  Fields[0]->Distribute(&(tx_vec.GetBlock(0)));

  FieldNames.push_back("Potential");
  Fields.push_back(new ParGridFunction);
  Fields[1]->MakeRef(fespaceH1, x_vec.GetBlock(1), 0);
  Fields[1]->Distribute(&(tx_vec.GetBlock(1)));
};


DarcyEMProblem::~DarcyEMProblem(){
   if(JForm     != NULL) delete JForm;
   if(VForm     != NULL) delete VForm;
   if(darcyEMOp != NULL) delete darcyEMOp;
   if(darcyEMPr != NULL) delete darcyEMPr;
   if(JJ_Form   != NULL) delete JJ_Form;
   if(JJForm    != NULL) delete JJForm;
   if(JVForm    != NULL) delete JVForm;
   if(VJForm    != NULL) delete VJForm;
};
#endif
