/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWCHEMISTRYSINK
#define POROUSFLOWCHEMISTRYSINK

#include "Kernel.h"
#include "PorousFlowDictator.h"

class PorousFlowChemistrySink;

template<>
InputParameters validParams<PorousFlowChemistrySink>();

/**
 * Sink of chemicals.  THIS IS CURRENTLY VERY SPECIALISED!!
 */
class PorousFlowChemistrySink : public Kernel
{
public:
  PorousFlowChemistrySink(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Generalised concentrations in each phase
  const MaterialProperty<std::vector<std::vector<Real> > > & _gen_concentration;

  /// Derivative of the generalised concentrations wrt PorousFlow variables
  const MaterialProperty<std::vector<std::vector<std::vector<Real> > > > & _dgen_concentration_dvar;

  /// PorousFlow Dictator UserObject
  const PorousFlowDictator & _dictator;

  /// The phase for this Kernel
  const unsigned int _phase;

  /// The primary species for this Kernel
  const unsigned int _primary_species;

  /// The coefficient for the sink
  const Real _coeff;

  /// The equilibrium concentration
  const Real _ceq;
};

#endif // POROUSFLOWCHEMISTRYSINK
