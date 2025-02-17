
// Define a nonlinear integrator that computes (f(u), v) and its linearized
// operator, (u df(u), v).
//
// Note that the action (f(u), v) can be computed using DomainLFIntegrator
// and the Jacobian matrix linearized operator can be computed using
// MassIntegrator with the appropriate coefficients.
class NonlinearMassIntegrator : public NonlinearFormIntegrator
{
   FiniteElementSpace &fes;
   GridFunction gf;
   Array<int> dofs;

public:
   NonlinearMassIntegrator(FiniteElementSpace &fes_) : fes(fes_), gf(&fes) { }

   virtual void AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      const Vector &elfun, Vector &elvect)
   {
      fes.GetElementDofs(Tr.ElementNo, dofs);
      gf.SetSubVector(dofs, elfun);
      NonlinearGridFunctionCoefficient coeff(gf, f);
      DomainLFIntegrator integ(coeff);
      integ.AssembleRHSElementVect(el, Tr, elvect);
   }

   virtual void AssembleElementGrad(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    const Vector &elfun, DenseMatrix &elmat)
   {
      fes.GetElementDofs(Tr.ElementNo, dofs);
      gf.SetSubVector(dofs, elfun);
      NonlinearGridFunctionCoefficient coeff(gf, df);
      MassIntegrator integ(coeff);
      integ.AssembleElementMatrix(el, Tr, elmat);
   }
};
