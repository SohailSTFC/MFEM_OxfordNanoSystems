#ifndef STOKESPROBLEMFACADECLASS_HPP
#define STOKESPROBLEMFACADECLASS_HPP 

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

#include "BoundaryAndInitialSolution.hpp"
#include "boundaryConditions.hpp"
#include "referenceTypes.hpp"

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


class StokesProbFacade
{
  private:
    int dim;
    ParFiniteElementSpace *H1_Uspace=NULL; //H1 elements of Order p+1
    ParFiniteElementSpace *H1_Pspace=NULL; //H1 elements of Order p

    Array<Coefficient*>       ess_bdr_PCoeffs;
    Array<VectorCoefficient*> ess_bdr_UCoeffs;
    Array<int> U_EssBCTags_markers, P_EssBCTags_markers;
  public:

  //Construct the Stoke problem Facade
  StokesProbFacade(ParFiniteElementSpace *h1_Uspace
                 , ParFiniteElementSpace *h1_Pspace
                 , int DIM);

  //Set the Finite element spaces
  void SetFESpaces(I_map<ParFiniteElementSpace*> *feSpaces
                 , vector<std::string> *FieldNames
                 , vector<int> *sv_flag);

  //Set the Boundary conditions
  void SetBCs(I_map<Array<int>*> *essBCtag_markers
            , I_map<Array<Coefficient*>*> *ess_bdr_sCoeffs
            , I_map<Array<VectorCoefficient*>*> *ess_bdr_vCoeffs);

  //Set the block bilinear forms
  void SetBlockBilinearForms(const I_map<ParFiniteElementSpace*> *feSpaces
                           , IJ_map<ParMixedBilinearForm*> *MBLForms
                           , IJ_map<int> *TrFlag);
};

//
// Construct the problem
// Facade used to construct
// the equation system
//
StokesProbFacade::StokesProbFacade(ParFiniteElementSpace *h1_Uspace
                                 , ParFiniteElementSpace *h1_Pspace
			                     , int DIM):
               H1_Uspace(h1_Uspace), H1_Pspace(h1_Pspace), dim(DIM){};

//
// Sets the FESpaces, the Field
// names, and whether the space
// is a scalar or vector space
//
void StokesProbFacade::SetFESpaces(I_map<ParFiniteElementSpace*> *feSpaces
                                 , vector<std::string> *FieldNames
                                 , vector<int> *sv_flag)
{
  (*feSpaces)[0] = new ParFiniteElementSpace(*H1_Uspace);
  (*feSpaces)[1] = new ParFiniteElementSpace(*H1_Pspace);
  FieldNames->push_back("Velocity");
  FieldNames->push_back("Pressure");
  sv_flag->push_back(1);
  sv_flag->push_back(0);
};


//
// Sets the Boundary condition 
// and coefficients for each of
// the tagged surfaces
//
void StokesProbFacade::SetBCs(I_map<Array<int>*> *essBCtag_markers
                            , I_map<Array<Coefficient*>*> *ess_bdr_sCoeffs
                            , I_map<Array<VectorCoefficient*>*> *ess_bdr_vCoeffs)
{
    int nUTags = H1_Uspace->GetMesh()->bdr_attributes.Max();;
    int nPTags = H1_Pspace->GetMesh()->bdr_attributes.Max();;
    U_EssBCTags_markers = Array<int>(nUTags);
	P_EssBCTags_markers = Array<int>(nPTags);

    //Set the tagged surfaces
    U_EssBCTags_markers = 1;
    P_EssBCTags_markers = 0;

    U_EssBCTags_markers[7] = 0;
    U_EssBCTags_markers[5] = 0;
    P_EssBCTags_markers[4] = 1;

    //Construct the coeeficients
    Vector ZeroV(dim);
	ZeroV = 0.0;
    ess_bdr_UCoeffs = Array<VectorCoefficient*>(nUTags);
    ess_bdr_PCoeffs = Array<Coefficient*>(nPTags);

    //Initially set Zero Coeffs
    for(int I=0; I<nUTags; I++) ess_bdr_UCoeffs[I] = new VectorConstantCoefficient(ZeroV);
    for(int I=0; I<nPTags; I++) ess_bdr_PCoeffs[I] = new ConstantCoefficient(0.00);

    //Delete and reset the Coeffs
    delete ess_bdr_UCoeffs[4];
    ess_bdr_UCoeffs[4] = new VectorFunctionCoefficient(dim, InletFunc);

    //Make copies and send them to the
    //equation system
    (*ess_bdr_sCoeffs)[1] = new Array<Coefficient*>(ess_bdr_PCoeffs);
    (*ess_bdr_vCoeffs)[0] = new Array<VectorCoefficient*>(ess_bdr_UCoeffs);
};


//
// Sets the block bilinear forms
// for the Stokes Problem on H1-H1
// mixed space of two different Orders
//	
void StokesProbFacade::SetBlockBilinearForms(const I_map<ParFiniteElementSpace*> *feSpaces
                                           , IJ_map<ParMixedBilinearForm*> *MBLForms
                                           , IJ_map<int> *TrFlag)
{
  (*MBLForms)[{0,0}] = new ParMixedBilinearForm((*feSpaces)[0], (*feSpaces)[0]);
  (*TrFlag)[{0,0}]   = 0;

  (*MBLForms)[{0,1}] = new ParMixedBilinearForm((*feSpaces)[0], (*feSpaces)[1]); 
  (*TrFlag)[{0,1}]   = 0;

  (*MBLForms)[{1,0}] = new ParMixedBilinearForm((*feSpaces)[0], (*feSpaces)[1]);
  (*TrFlag)[{0,1}]   = 1;

  ConstantCoefficient One(1.00);
  ConstantCoefficient Mu(mu);

  //Mixed Bilinear forms (for solution)
  (*MBLForms)[{0,0}] ->AddDomainIntegrator(new VectorDiffusionIntegrator(mu));
  (*MBLForms)[{0,0}] ->Assemble();

  (*MBLForms)[{0,1}]->AddDomainIntegrator(new VectorDivergenceIntegrator(One));
  (*MBLForms)[{0,1}]->Assemble();

  (*MBLForms)[{1,0}]->AddDomainIntegrator(new VectorDivergenceIntegrator(One));
  (*MBLForms)[{1,0}]->Assemble();
};
#endif