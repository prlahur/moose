/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowConstantBulkProperties.h"

namespace PorousFlowConstantBulkProperties
{
std::string
fluidName()
{
  return "Constant bulk";
}

Real
density(Real pressure, Real density_p0, Real bulk_modulus)
{
  return density_p0 * std::exp(pressure / bulk_modulus);
}

Real
dDensity_dP(Real pressure, Real density_p0, Real bulk_modulus)
{
  return density(pressure, density_p0, bulk_modulus) / bulk_modulus;
}
}
