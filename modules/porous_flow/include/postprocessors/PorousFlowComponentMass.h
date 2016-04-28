/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWCOMPONENTMASS_H
#define POROUSFLOWCOMPONENTMASS_H

#include "ElementIntegralVariablePostprocessor.h"
#include "PorousFlowDictator.h"

// Forward Declarations
class PorousFlowComponentMass;

template<>
InputParameters validParams<PorousFlowComponentMass>();

/**
 * Postprocessor produces the mass of a given fluid component in a region
 */
class PorousFlowComponentMass: public ElementIntegralPostprocessor
{
public:
  PorousFlowComponentMass(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  /// Holds info on the PorousFlow variables
  const PorousFlowDictator & _dictator_UO;
  /// The component number that this Postprocessor applies to
  std::vector<unsigned int> _component_index;
  /// The (optional) phase indices that this Postprocessor applies to. Default is all phases
  std::vector<unsigned int> _phase_index;
  const MaterialProperty<Real> & _porosity;
  const MaterialProperty<std::vector<Real> > & _fluid_density;
  const MaterialProperty<std::vector<Real> > & _fluid_saturation;
  const MaterialProperty<std::vector<std::vector<Real> > > & _mass_fraction;
  /// Saturation threshold - only fluid mass at saturations below this are calculated
  const Real _saturation_threshold;
};

#endif //POROUSFLOWCOMPONENTMASS_H
