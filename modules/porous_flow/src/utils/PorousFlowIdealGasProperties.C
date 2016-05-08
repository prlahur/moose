/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowIdealGasProperties.h"

namespace PorousFlowIdealGasProperties
{

std::string
fluidName()
{
  return "Ideal";
}

Real
density(Real pressure, Real temperature, Real molar_mass)
{
  return pressure * molar_mass / (_R * (temperature + _t_c2k));
}


Real
dDensity_dP(Real temperature, Real molar_mass)
{
  return molar_mass / (_R * (temperature + _t_c2k));
}

Real
dDensity_dT(Real pressure, Real temperature, Real molar_mass)
{
  Real tk = temperature + _t_c2k;
  return - pressure * molar_mass / (_R * tk * tk);
}
}
