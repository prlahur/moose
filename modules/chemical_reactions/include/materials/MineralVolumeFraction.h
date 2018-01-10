/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef MINERALVOLUMEFRACTION_H
#define MINERALVOLUMEFRACTION_H

#include "Material.h"

class MineralVolumeFraction;

template <>
InputParameters validParams<MineralVolumeFraction>();

/**
 * The volume fraction of each coupled mineral species calculated
 * using the molar concentration and molar volume.
 */
class MineralVolumeFraction : public Material
{
public:
  MineralVolumeFraction(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Vector of coupled mineral species concentration
  std::vector<const VariableValue *> _coupled_species_conc;
  /// Number of coupled mineral species
  const unsigned int _num_coupled_species;
  /// Molar volume of coupled mineral species
  const std::vector<Real> _molar_volume;
  /// Mineral volume fraction
  MaterialProperty<Real> & _mineral_volume_frac;
  /// Volume of the current element
  const Real & _current_elem_volume;
};

#endif // MINERALVOLUMEFRACTION_H
