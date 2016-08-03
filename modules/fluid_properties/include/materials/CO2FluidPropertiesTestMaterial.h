/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef CO2FLUIDPROPERTIESTESTMATERIAL_H
#define CO2FLUIDPROPERTIESTESTMATERIAL_H

#include "CO2FluidProperties.h"
#include "Material.h"

class CO2FluidPropertiesTestMaterial;

template<>
InputParameters validParams<CO2FluidPropertiesTestMaterial>();

/**
 * Computes values of pressure and its derivatives using (u, v) formulation
 */
class CO2FluidPropertiesTestMaterial : public Material
{
public:
  CO2FluidPropertiesTestMaterial(const InputParameters & parameters);
  virtual ~CO2FluidPropertiesTestMaterial();

protected:
  virtual void computeQpProperties();

  /// Pressure (Pa)
  const VariableValue & _pressure;
  /// Temperature (K)
  const VariableValue & _temperature;
  /// Sublimation pressure (Pa)
  MaterialProperty<Real> & _psub;
  /// Melting pressure (Pa)
  MaterialProperty<Real> & _pmelt;
  /// Vapour pressure (Pa)
  MaterialProperty<Real> & _pvap;
  /// Saturated vapour density (kg/m^3)
  MaterialProperty<Real> & _rhovap;
  /// Saturated liquid density (kg/m^3)
  MaterialProperty<Real> & _rhosat;
  /// Partial density (kg/m^3)
  MaterialProperty<Real> & _rhopartial;

  /// CO2 Fluid properties UserObject
  const CO2FluidProperties & _fp;
  /// Bool flag to calculate sublimation pressure
  const bool _sublimation;
  /// Bool flag to calculate melting pressure
  const bool _melting;
  /// Bool flag to calculate vapour pressure
  const bool _vapour;
};

#endif /* CO2FLUIDPROPERTIESMATERIAL_H */
