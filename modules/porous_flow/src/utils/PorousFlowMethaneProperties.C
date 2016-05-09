/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowMethaneProperties.h"
#include "PorousFlowIdealGasProperties.h"

namespace PorousFlowMethaneProperties
{

std::string
fluidName()
{
  return "Methane";
}

Real
molarMass()
{
  return _Mch4;
}

Real
criticalPressure()
{
  return _p_critical;
}

Real
criticalTemperature()
{
  return _t_critical;
}

Real
density(Real pressure, Real temperature)

{
  return PorousFlowIdealGasProperties::density(pressure, temperature, _Mch4);
}

Real
viscosity(Real temperature)
{
  Real a[6] = {2.968267e-1, 3.711201e-2, 1.218298e-5, -7.02426e-8, 7.543269e-11,
    -2.7237166e-14};

  Real viscosity;
  Real tk = temperature + _t_c2k;
  Real tk2 = tk * tk;
  Real tk3 = tk2 * tk;
  Real tk4 = tk3 * tk;
  Real tk5 = tk4 * tk;

  viscosity = a[0] + a[1] * tk + a[2] * tk2 + a[3] * tk3 + a[4] * tk4 + a[5] * tk5;

  return viscosity * 1.e-6;
}

Real
dDensity_dP(Real temperature)
{
  return PorousFlowIdealGasProperties::dDensity_dP(temperature, _Mch4);
}

Real
dDensity_dT(Real pressure, Real temperature)
{
  return PorousFlowIdealGasProperties::dDensity_dT(pressure, temperature, _Mch4);
}

Real
dViscosity_dT(Real temperature)
{
  Real a[6] = {2.968267e-1, 3.711201e-2, 1.218298e-5, -7.02426e-8, 7.543269e-11,
    -2.7237166e-14};

  Real dviscosity;
  Real tk = temperature + _t_c2k;
  Real tk2 = tk * tk;
  Real tk3 = tk2 * tk;
  Real tk4 = tk3 * tk;

  dviscosity = a[1] + 2.0 * a[2] * tk + 3.0 * a[3] * tk2 + 4.0 * a[4] * tk3 + 5.0 * a[5] * tk4;

  return dviscosity * 1.e-6;
}

std::vector<Real>
henryConstants()
{
  std::vector<Real> ch4henry;
  ch4henry.push_back(-10.44708);
  ch4henry.push_back(4.66491);
  ch4henry.push_back(12.12986);

  return ch4henry;
}
}
