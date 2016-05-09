/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWCONSTANTBULKPROPERTIES_H
#define POROUSFLOWCONSTANTBULKPROPERTIES_H

#include "MooseTypes.h"

namespace PorousFlowConstantBulkProperties
{
  /**
   * Fluid name
   * @return fluid Name
   */
  std::string fluidName();

  /**
   * Density as a function of  pressure and temperature assuming an
   * ideal gas.
   *
   * @param pressure pressure (Pa)
   * @param density_p0 density at zero pressure (kg/m^3)
   * @param bulk_modulus bulk modulus (Pa)
   * @return density (kg/m^3)
   */
  Real density(Real pressure,  Real density_p0, Real bulk_modulus);

  /**
   * Derivative of the density wrt pressure
   *
   * @param pressure pressure (Pa)
   * @param density_p0 density at zero pressure (kg/m^3)
   * @param bulk_modulus bulk modulus (Pa)
   * @return derivative of density (kg/m^3) wrt pressure
   */
  Real dDensity_dP(Real pressure, Real density_p0, Real bulk_modulus);
}

#endif // POROUSFLOWCONSTANTBULKPROPERTIES_H
