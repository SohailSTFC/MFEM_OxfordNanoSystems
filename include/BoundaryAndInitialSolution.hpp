#ifndef BOUNDARYANDINITIALSOLUTION_HPP
#define BOUNDARYANDINITIALSOLUTION_HPP 

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Define the analytical solution and forcing terms / boundary conditions
double pFun_ex(const Vector & x);
void fFun(const Vector & x, Vector & f);
double gFun(const Vector & x);
double f_natural(const Vector & x);


// Change if needed
double pFun_ex(const Vector & x)
{
   double xi(x(0));
   double yi(x(1));
   double zi(0.0);

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

double gFun(const Vector & x)
{
   if (x.Size() == 3)
   {
      return -pFun_ex(x);
   }
   else
   {
      return 0.0;
   }
}

double f_natural(const Vector & x)
{
   return (-pFun_ex(x));
}
#endif