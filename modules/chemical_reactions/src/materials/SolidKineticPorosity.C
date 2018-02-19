/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "SolidKineticPorosity.h"

template <>
InputParameters
validParams<SolidKineticPorosity>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("base_porosity", "Base porosity with no mineral species");
  params.addClassDescription("Porosity including contribution due to solid mineral species");
  return params;
}

SolidKineticPorosity::SolidKineticPorosity(const InputParameters & parameters)
  : Material(parameters),
    _base_porosity(getParam<Real>("base_porosity")),
    _mineral_volume_frac(getMaterialProperty<Real>("mineral_volume_frac")),
    _porosity(declareProperty<Real>("porosity"))
{
}

void
SolidKineticPorosity::initQpStatefulProperties()
{
  computeQpProperties();
}

void
SolidKineticPorosity::computeQpProperties()
{
  Real porosity = _base_porosity - _mineral_volume_frac[_qp];

  // Porosity cannot be less than zero or greater than 1
  porosity = porosity < 0.0 ? 0 : (porosity > 1.0 ? 1 : porosity);

  _porosity[_qp] = porosity;
}
