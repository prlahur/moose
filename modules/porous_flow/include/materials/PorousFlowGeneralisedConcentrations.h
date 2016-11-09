/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWGENERALISEDCONCENTRATIONS_H
#define POROUSFLOWGENERALISEDCONCENTRATIONS_H

#include "PorousFlowMaterialVectorBase.h"

//Forward Declarations
class PorousFlowGeneralisedConcentrations;

template<>
InputParameters validParams<PorousFlowGeneralisedConcentrations>();

/**
 * Material designed to form a std::vector<std::vector>
 * of mass fractions from the individual mass fraction variables
 */
class PorousFlowGeneralisedConcentrations : public PorousFlowMaterialVectorBase
{
public:
  PorousFlowGeneralisedConcentrations(const InputParameters & parameters);

protected:
  /// Mass fraction matrix
  MaterialProperty<std::vector<std::vector<Real> > > & _gen_conc;

  /// Old value of mass fraction matrix
  MaterialProperty<std::vector<std::vector<Real> > > & _gen_conc_old;

  /// Gradient of the mass fraction matrix at the quad points
  MaterialProperty<std::vector<std::vector<RealGradient> > > & _grad_gen_conc;

  /// Derivative of the mass fraction matrix with respect to the porous flow variables
  MaterialProperty<std::vector<std::vector<std::vector<Real> > > > & _dgen_conc_dvar;

  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  /**
   * Builds the mass-fraction variable matrix at the quad point
   * @param qp the quad point
   */
  void build_gen_conc(unsigned int qp);

  /**
   * Number of mass-fraction variables provided by the user
   * This needs to be _num_phases*_num_primary
   */
  const unsigned int _num_passed_gc_vars;

  /// the variable number of the mass-fraction variables
  std::vector<unsigned int> _gc_vars_num;

  /// the mass-fraction variables
  std::vector<const VariableValue *> _gc_vars;

  /// the gradient of the mass-fraction variables
  std::vector<const VariableGradient *> _grad_gc_vars;
};

#endif //POROUSFLOWGENERALISEDCONCENTRATIONS_H
