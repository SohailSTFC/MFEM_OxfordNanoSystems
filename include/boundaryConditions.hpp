#ifndef BOUNDARYCONDITIONS_HPP
#define BOUNDARYCONDITIONS_HPP

#include "mfem.hpp"
#include "../../mfem-4.6/general/forall.hpp"

using namespace mfem;
using namespace std;

//
// Applies the value at the boundary from a reference Vector/GridFunction
// to the solution vector for use in non-Hmogenous Dirchelet BCs
//
void applyDirchValues(const Vector &k, Vector &y, Array<int> dofs)
{
  if(dofs.Size() > 0){ //Only apply if there are constrained DOF's
    const bool use_dev = dofs.UseDevice() || k.UseDevice() || y.UseDevice();
    const int n = dofs.Size();
    // Use read+write access for X - we only modify some of its entries
    auto d_X = y.ReadWrite(use_dev);
    auto d_y = k.Read(use_dev);
    auto d_dofs = dofs.Read(use_dev);
    mfem::forall_switch(use_dev, n, [=] MFEM_HOST_DEVICE (int i)
    {
      const int dof_i = d_dofs[i];
      if (dof_i >= 0)   d_X[dof_i]    =  d_y[dof_i];
      if (!(dof_i >= 0))d_X[-1-dof_i] = -d_y[-1-dof_i];
    });
  }
};


//
// Applies the value-elimination at the boundary from a
// to the residual or solution vector for use in non-Homogenous
// and homogenous Dirchelet BCs
//
void applyDirchElimination(Vector &y, Array<int> dofs)
{
  if(dofs.Size() > 0){ //Only apply if there are constrained DOF's
    const bool use_dev = dofs.UseDevice() || y.UseDevice();
    const int n = dofs.Size();
    // Use read+write access for X - we only modify some of its entries
    auto d_X = y.ReadWrite(use_dev);
    auto d_dofs = dofs.Read(use_dev);
    mfem::forall_switch(use_dev, n, [=] MFEM_HOST_DEVICE (int i)
    {
      const int dof_i = d_dofs[i];
      if (dof_i >= 0)   d_X[dof_i]    =  0.0;
      if (!(dof_i >= 0))d_X[-1-dof_i] =  0.0;
    });
  }
};


//
// Sets the entries of one vector from another vector
//
void setValues(const Vector &k, Vector &y){ y = k; };


#endif
