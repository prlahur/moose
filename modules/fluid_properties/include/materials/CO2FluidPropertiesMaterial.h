/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef CO2FLUIDPROPERTIESMATERIAL_H
#define CO2FLUIDPROPERTIESMATERIAL_H

#include "CO2FluidProperties.h"
#include "Material.h"

class CO2FluidPropertiesMaterial;

template<>
InputParameters validParams<CO2FluidPropertiesMaterial>();

/**
 * Computes values of pressure and its derivatives using (u, v) formulation
 */
class CO2FluidPropertiesMaterial : public Material
{
public:
  CO2FluidPropertiesMaterial(const InputParameters & parameters);
  virtual ~CO2FluidPropertiesMaterial();

protected:
  virtual void computeQpProperties();

  /// Pressure (Pa)
  const VariableValue & _pressure;
  /// Temperature (K)
  const VariableValue & _temperature;
  MaterialProperty<Real> & _rho;
  MaterialProperty<Real> & _mu;

  /// CO2 Fluid properties
  const CO2FluidProperties & _fp;
};

#endif /* CO2FLUIDPROPERTIESMATERIAL_H */
