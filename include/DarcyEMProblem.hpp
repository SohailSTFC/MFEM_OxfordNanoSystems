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

    //Boundary conditions
    Array<int> ess_tdof_list;

    //The permiativity of space
    real_t sigma;

    //The linear and Bilinear forms (residual and jacobian)
    ParBilinearForm *Jacobian;
    ParLinearForm   *AVForm;

  public:
	//The 
    DarcyEMProblem(ParFiniteElementSpace *f1, ParFiniteElementSpace *f2, real_t sig);
	
};


//
//
//  Implementation of the problem class
//
//

//The constructor
DarcyEMProblem::DarcyEMProblem(ParFiniteElementSpace *f1RT, ParFiniteElementSpace *f2M, real_t sig)
: fespaceRT(f1RT), fespaceMixed(f2M), sigma(sig)
{
	
};
#endif