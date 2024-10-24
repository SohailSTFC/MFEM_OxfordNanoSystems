#ifndef BOUNDARYANDINITIALSOLUTION_HPP
#define BOUNDARYANDINITIALSOLUTION_HPP 

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Define the analytical solution and forcing terms / boundary conditions
void uFun_ex(const Vector & x, Vector & u);
real_t pFun_ex(const Vector & x);
void fFun(const Vector & x, Vector & f);
real_t gFun(const Vector & x);
real_t f_natural(const Vector & x);

void uFun_ex(const Vector & x, Vector & u)
{
   real_t xi(x(0));
   real_t yi(x(1));
   real_t zi(0.0);
   if (x.Size() == 3)
   {
      zi = x(2);
   }

   u(0) = - exp(xi)*sin(yi)*cos(zi);
   u(1) = - exp(xi)*cos(yi)*cos(zi);

   if (x.Size() == 3)
   {
      u(2) = exp(xi)*sin(yi)*sin(zi);
   }
}

// Change if needed
real_t pFun_ex(const Vector & x)
{
   real_t xi(x(0));
   real_t yi(x(1));
   real_t zi(0.0);

   if (x.Size() == 3)
   {
      zi = x(2);
   }

   return exp(xi)*sin(yi)*cos(zi);
}

void fFun(const Vector & x, Vector & f)
{
   f = 0.0;
}

real_t gFun(const Vector & x)
{
   if (x.Size() == 3)
   {
      return -pFun_ex(x);
   }
   else
   {
      return 0;
   }
}

real_t f_natural(const Vector & x)
{
   return (-pFun_ex(x));
}

//
// Boundary conditions
//

void current_dirichlet_bc(const Vector & x, Vector & u){
  u = 0.00;	  
}
#endif