/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWMETHANEPROPERTIES_H
#define POROUSFLOWMETHANEPROPERTIES_H

#include "MooseTypes.h"
#include "MooseError.h"

namespace PorousFlowMethaneProperties
{
  /**
   * Fluid name
   * @return fluid Name
   */
  std::string fluidName();

  /**
   * Methane molar mass.
   * @return molar mass (kg/mol)
   */
  Real molarMass();

  /**
   * CH4 gas density as a function of  pressure and temperature assuming an
   * ideal gas
   *
   * @param pressure gas pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return density (kg/m^3)
   */
  Real density(Real pressure, Real temperature);

  /**
   * CH4 gas viscosity as a function of temperature.
   * From Irvine Jr, T. F. and Liley, P. E. (1984) Steam and Gas Tables with
   * Computer Equations.
   *
   * @param temperature fluid temperature (C)
   * @return viscosity (Pa.s)
   */
  Real viscosity(Real temperature);

  /**
   * Derivative of the density of gaseous CH4 as a function of
   * pressure.
   *
   * @param temperature gas temperature (C)
   * @return derivative of CH4 density (kg/m^3) with respect to pressure
   */
  Real dDensity_dP(Real temperature);

  /**
   * Derivative of the density of gaseous CH4 as a function of
   * temperature.
   *
   * @param pressure gas pressure (Pa)
   * @param temperature gas temperature (C)
   * @return derivative of CH4 density (kg/m^3) with respect to temperature
   */
  Real dDensity_dT(Real pressure, Real temperature);

  /**
   * Derivative of the viscosity of gaseous CH4 as a function of temperature
   *
   * @param temperature gas temperature (C)
   * @return derivative of CH4 viscosity wrt temperature
   */
  Real dViscosity_dT(Real temperature);

  /**
   * Henry's law constant coefficients for dissolution of CH4 into water.
   * From Guidelines on the Henry's constant and vapour
   * liquid distribution constant for gases in H20 and D20 at high
   * temperatures, IAPWS (2004).
   *
   * @return constants for Henry's constant (-)
   */
  std::vector<Real> henryConstants();

  /// Conversion of temperature from Celcius to Kelvin
  const Real _t_c2k = 273.15;
  /// Molar mass of pure CH4
  const Real _Mch4 = 16.0425e-3;
}

#endif // POROUSFLOWMETHANEPROPERTIES_H
