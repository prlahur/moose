/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef SOLIDKINETICPOROSITY_H
#define SOLIDKINETICPOROSITY_H

#include "Material.h"

class SolidKineticPorosity;

template <>
InputParameters validParams<SolidKineticPorosity>();

/**
 * The porosity of the porous medium including the contribution due to
 * solid mineral species.
 *
 * In this case, porosity is defined as one minus the sum of the total unreactive
 * rock volume fraction and the total reactive mineral volume fraction. Total reactive
 * mineral volume fraction is calculated in the MineralVolumeFraction material.
 *
 * The initial porosity is input using the base_porosity input parameter, which represents
 * the porosity due to the unreactive rock (supplied as 1 - unreactive volume fraction)
 */
class SolidKineticPorosity : public Material
{
public:
  SolidKineticPorosity(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Base porosity (without mineral contribution)
  const Real _base_porosity;
  /// Total mineral volume fraction
  const MaterialProperty<Real> & _mineral_volume_frac;
  /// Porosity
  MaterialProperty<Real> & _porosity;
};

#endif // SOLIDKINETICPOROSITY_H
