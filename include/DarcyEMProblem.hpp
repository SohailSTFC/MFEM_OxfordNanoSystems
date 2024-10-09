#ifndef DARCYEMPROBLEM_HPP
#define DARCYEMPROBLEM_HPP 

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

//
//
//  The problem class
//
//
class DarcyEMProblem
{
  protected:
    //Finite element spaces
    ParFiniteElementSpace *fespaceRT;    //Raviart Thomas elements
    ParFiniteElementSpace *fespaceMixed; //Mixed finite element space

    //Block Matrix structure offsets
    Array<int> block_offsets;     // number of variables + 1 (2-variables J and v)
    Array<int> block_trueOffsets; // number of variables + 1 (2-variables J and v)

    //The permiativity of space
    real_t sigma;

    //The linear and Bilinear forms (residual and jacobian)
    ParBilinearForm JJForm;
    ParLinearForm   JForm;
    ParLinearForm   VForm;

    //Form operators
   BlockOperator *darcyOp = new BlockOperator(block_trueOffsets);
   Array<int> empty_tdof_list;  // empty
   OperatorPtr opM, opB;
   TransposeOperator *Bt = NULL;


  public:
	//The constructor
    DarcyEMProblem(ParFiniteElementSpace *f1, ParFiniteElementSpace *f2, real_t sig);

    //Compute the RHS of the Linear system
	void Mult(const Vector &r, Vector &y);
};


//
//
//  Implementation of the problem class
//
//

//The constructor
DarcyEMProblem::DarcyEMProblem(ParFiniteElementSpace *f1RT, ParFiniteElementSpace *f2M, real_t sig)
: fespaceRT(f1RT), fespaceMixed(f2M), JJForm(fespaceRT), sigma(sig)
{




   // 8. Define the two BlockStructure of the problem.  block_offsets is used
   //    for Vector based on dof (like ParGridFunction or ParLinearForm),
   //    block_trueOffstes is used for Vector based on trueDof (HypreParVector
   //    for the rhs and solution of the linear system).  The offsets computed
   //    here are local to the processor.
/*   Array<int> block_offsets(3); // number of variables + 1
   block_offsets[0] = 0;
   block_offsets[1] = R_space->GetVSize();
   block_offsets[2] = W_space->GetVSize();
   block_offsets.PartialSum();

   Array<int> block_trueOffsets(3); // number of variables + 1
   block_trueOffsets[0] = 0;
   block_trueOffsets[1] = R_space->TrueVSize();
   block_trueOffsets[2] = W_space->TrueVSize();
   block_trueOffsets.PartialSum();

   block_offsets
   block_trueOffsets*/
};


void DarcyEMProblem::Mult(const Vector &r, Vector &y)
{
	

};
#endif