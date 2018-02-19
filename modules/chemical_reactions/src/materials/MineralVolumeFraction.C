/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MineralVolumeFraction.h"
#include "Assembly.h"

template <>
InputParameters
validParams<MineralVolumeFraction>()
{
  InputParameters params = validParams<Material>();
  params.addCoupledVar("mineral_species", "The coupled mineral species");
  params.addRequiredParam<std::vector<Real>>("molar_volume",
                                             "Molar volume of coupled mineral species");
  params.addClassDescription("Total volume fraction of coupled mineral species");
  return params;
}

MineralVolumeFraction::MineralVolumeFraction(const InputParameters & parameters)
  : Material(parameters),
    _num_coupled_species(coupledComponents("mineral_species")),
    _molar_volume(getParam<std::vector<Real>>("molar_volume")),
    _mineral_volume_frac(declareProperty<Real>("mineral_volume_frac")),
    _current_elem_volume(_assembly.elemVolume())
{
  _coupled_species_conc.resize(_num_coupled_species);

  for (unsigned int i = 0; i < _num_coupled_species; ++i)
    _coupled_species_conc[i] = &coupledValue("mineral_species", i);

  // Make sure that each coupled species has a corresponding molar volume
  if (_molar_volume.size() != _num_coupled_species)
    mooseError(
        "The number of entries in molar_volume is not equal to the number of coupled_species");
}

void
MineralVolumeFraction::computeQpProperties()
{
  Real species_vf = 0.0;

  for (unsigned int i = 0; i < _num_coupled_species; ++i)
    species_vf += (*_coupled_species_conc[i])[_qp] * _molar_volume[i] / _current_elem_volume;

  // Volume fraction cannot be less than 0 or greater than 1
  species_vf = species_vf < 0.0 ? 0 : (species_vf > 1.0 ? 1 : species_vf);

  _mineral_volume_frac[_qp] = species_vf;
}
