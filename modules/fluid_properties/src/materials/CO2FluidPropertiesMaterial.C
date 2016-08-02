/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#include "CO2FluidPropertiesMaterial.h"

template<>
InputParameters validParams<CO2FluidPropertiesMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("pressure", "Pressure (Pa)");
  params.addRequiredCoupledVar("temperature", "Temperature");
  params.addRequiredParam<UserObjectName>("fp", "The name of the user object for fluid properties");
  return params;
}

CO2FluidPropertiesMaterial::CO2FluidPropertiesMaterial(const InputParameters & parameters) :
    Material(parameters),
    _pressure(coupledValue("pressure")),
    _temperature(coupledValue("temperature")),

    _rho(declareProperty<Real>("density")),
    _mu(declareProperty<Real>("viscosity")),

    _fp(getUserObject<CO2FluidProperties>("fp"))
{
}

CO2FluidPropertiesMaterial::~CO2FluidPropertiesMaterial()
{
}

void
CO2FluidPropertiesMaterial::computeQpProperties()
{
  _rho[_qp] = _fp.rho(_pressure[_qp], _temperature[_qp]);
  _mu[_qp] = _fp.mu(_rho[_qp], _temperature[_qp]);
}
