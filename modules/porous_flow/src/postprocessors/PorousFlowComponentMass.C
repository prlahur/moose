/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowComponentMass.h"

template<>
InputParameters validParams<PorousFlowComponentMass>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredParam<std::vector<unsigned int> >("component_index", "The index(s) corresponding to the component for this Postprocessor");
  params.addRequiredParam<UserObjectName>("PorousFlowDictator_UO", "The UserObject that holds the list of PorousFlow variable names.");
  params.addParam<std::vector<unsigned int> >("phase_index", "The index(s) of the fluid phase that this Postprocessor applies to. Default is all phases");
  params.addRangeCheckedParam<Real>("saturation_threshold", 1.0, "saturation_threshold >= 0 & saturation_threshold <= 1", "The saturation threshold below which the mass is calculated for a specific phase. Default is 1.0. Note: only one phase_index can be entered");
  params.addClassDescription("Calculates the mass of a fluid component");
  return params;
}

PorousFlowComponentMass::PorousFlowComponentMass(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),

  _dictator_UO(getUserObject<PorousFlowDictator>("PorousFlowDictator_UO")),
  _component_index(getParam<std::vector<unsigned int> >("component_index")),
  _phase_index(getParam<std::vector<unsigned int> >("phase_index")),
  _porosity(getMaterialProperty<Real>("PorousFlow_porosity_nodal")),
  _fluid_density(getMaterialProperty<std::vector<Real> >("PorousFlow_fluid_phase_density")),
  _fluid_saturation(getMaterialProperty<std::vector<Real> >("PorousFlow_saturation_nodal")),
  _mass_fraction(getMaterialProperty<std::vector<std::vector<Real> > >("PorousFlow_mass_frac")),
  _saturation_threshold(getParam<Real>("saturation_threshold"))
{
  const unsigned int num_phases = _dictator_UO.numPhases();
  const unsigned int num_components = _dictator_UO.numComponents();

  /// Check that the number of components entered is not greater than the maximum number of components
  if (_component_index.size() > num_components)
    mooseError("The Dictator does not like it when you enter " << _component_index.size() << " components in the Postprocessor " << _name << " when only " << num_components << " are allowed. Try again");

  /// Check that the largest component_index is not greater than the number of components. Note that negative values
  /// throw an input parser error so we don't have to consider them
  if (_component_index.size() > 0)
  {
    unsigned int max_comp_num = * std::max_element(_component_index.begin(), _component_index.end());
    if (max_comp_num > num_components - 1)
      mooseError("The Dictator proclaims that the number of components in this simulation is " << num_components << " whereas you have used the Postprocessor PorousFlowComponentMass with component = " << max_comp_num << ".  The Dictator does not take such mistakes lightly");
  }

  /// Check that the number of phases entered is not more than the maximum possible phases
  if (_phase_index.size() > _dictator_UO.numPhases())
    mooseError("The Dictator decrees that the number of phases in this simulation is " << num_phases << " but you have entered " << _phase_index.size() << " phases in the Postprocessor " << _name);

  /// Also check that the phase indices entered are not greater than the number of phases to avoid a segfault. Note
  /// that the input parser takes care of negative inputs so we don't need to guard against them
  if (_phase_index.size() > 0)
  {
    unsigned int max_phase_num = * std::max_element(_phase_index.begin(), _phase_index.end());
    if (max_phase_num > num_phases - 1)
      mooseError("The Dictator proclaims that the phase index " << max_phase_num << " in the Postprocessor " << _name << " is greater than the largest phase index possible, which is " << num_phases - 1);
  }

  /// Using saturation_threshold only makes sense for a specific phase_index
  if (_saturation_threshold < 1 && _phase_index.size() != 1)
    mooseError("A single phase_index must be entered when prescribing a saturation_threshold in the Postprocessor " << _name);

  /// If phase_index is empty, we want to sum over all phases
  if (_phase_index.size() == 0)
    for (unsigned int i = 0; i < num_phases; ++i)
      _phase_index.push_back(i);
}

Real
PorousFlowComponentMass::computeQpIntegral()
{
  Real mass = 0.0;
  unsigned int ph, cp;
  for (unsigned int i = 0; i < _phase_index.size(); ++i)
  {
    ph = _phase_index[i];
    if (_fluid_saturation[_qp][ph] <= _saturation_threshold)
      for (unsigned int j = 0; j < _component_index.size(); ++j)
      {
        cp = _component_index[j];
        mass += _fluid_density[_qp][ph] * _fluid_saturation[_qp][ph] * _mass_fraction[_qp][ph][cp];
      }
  }

  return _porosity[_qp] * mass;
}
