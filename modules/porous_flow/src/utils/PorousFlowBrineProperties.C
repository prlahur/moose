/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowBrineProperties.h"
#include "PorousFlowWaterProperties.h"

namespace PorousFlowBrineProperties
{
  /// Reference constants used in to calculate thermophysical properties of water.
  const Real _Mh2o = PorousFlowWaterProperties::molarMass();
  const Real _t_c2k = PorousFlowWaterProperties::_t_c2k;

std::string
fluidName()
{
  return "Brine";
}

Real
molarMass()
{
  return _Mnacl;
}

Real
molarMassH2O()
{
  return _Mh2o;
}

Real
density(Real pressure, Real temperature, Real xnacl)
{
  Real n1, n2, n11, n12, n1x1, n20, n21, n22, n23, n2x1, Tv;
  Real water_density;

  /// The correlation requires the pressure in bar, not Pa.
  Real pbar = pressure * 1.0e-5;
  Real pbar2 = pbar * pbar;
  Real pbar3 = pbar2 * pbar;

  /// The correlation requires mole fraction. First calculate the average molar mass
  /// from the mass fraction, then compute the mole fraction
  Real Mbrine = 1.0 / (xnacl / _Mnacl + (1.0 - xnacl) / _Mh2o);
  Real Xnacl = xnacl / _Mnacl * Mbrine;

  n11 = -54.2958 - 45.7623 * std::exp(-9.44785e-4 * pbar);
  n21 = -2.6142 - 2.39092e-4 * pbar;
  n22 = 0.0356828 + 4.37235e-6 * pbar + 2.0566e-9 * pbar2;
  n1x1 = 330.47 + 0.942876 * std::sqrt(pbar) + 0.0817193 * pbar - 2.47556e-8 * pbar2
         + 3.45052e-10 * pbar3;
  n2x1 = -0.0370751 + 0.00237723 * std::sqrt(pbar) + 5.42049e-5 * pbar
         + 5.84709e-9 * pbar2 - 5.99373e-13 * pbar3;
  n12 = - n1x1 - n11;
  n20 = 1.0 - n21 * std::sqrt(n22);
  n23 = n2x1 - n20 -n21 * std::sqrt(1.0 + n22);

  /// The temperature Tv (C) where the brine has the same molar volume as pure water
  n1 = n1x1 + n11 * (1.0 - Xnacl) + n12 * (1.0 - Xnacl) * (1.0 - Xnacl);
  n2 = n20 + n21 * std::sqrt(Xnacl + n22) + n23 * Xnacl;
  Tv = n1 + n2 * temperature;

  /// The brine density is then given by the density of water at temperature Tv.
  water_density = PorousFlowWaterProperties::density(pressure, Tv);

  return water_density * Mbrine / _Mh2o;
}

Real
viscosity(Real pressure, Real temperature, Real xnacl)
{
  /// Correlation requires molar concentration (mol/kg)
  Real mol = xnacl / ((1.0 - xnacl) * _Mnacl);
  Real mol2 = mol * mol;
  Real mol3 = mol2 * mol;

  Real a = 1.0 + 0.0816 * mol + 0.0122 * mol2 + 0.128e-3 * mol3 + 0.629e-3 * temperature
           * (1.0 - std::exp(-0.7 * mol));

  /// The brine viscosity is then given by a multiplied by the viscosity of pure water
  Real water_density = PorousFlowWaterProperties::density(pressure, temperature);
  return a * PorousFlowWaterProperties::viscosity(temperature, water_density);
}

Real
haliteDensity(Real pressure, Real temperature)
{
  /// Correlation needs pressure in bar
  Real pbar = pressure * 1.e-5;
  /// Halite density at 0 Pa
  Real density_p0 = 2.17043e3 - 2.4599e-1 * temperature - 9.5797e-5 * temperature * temperature;
  /// Halite density as a function of pressure
  Real l = 5.727e-3 + 2.715e-3 * std::exp(temperature / 733.4);

  return density_p0 + l * pbar;
}

Real
haliteSolubility(Real temperature)
{
  Real solubility = (26.18 + 7.2e-3 * temperature + 1.06e-4 * temperature * temperature)/100.;

  return solubility;
}

Real
pSat(Real temperature, Real xnacl)
{
  /// Temperature in K
  Real tk = temperature + _t_c2k;

  /// Correlation requires molar concentration (mol/kg)
  Real mol = xnacl / ((1.0 - xnacl) * _Mnacl);
  Real mol2 = mol * mol;
  Real mol3 = mol2 * mol;
  Real mol4 = mol3 * mol;
  Real mol5 = mol4 * mol;

  Real a = 1.0 + 5.93582e-6 * mol - 5.19386e-5 * mol2 + 1.23156e-5 * mol3;
  Real b = 1.1542e-6 * mol + 1.41254e-7 * mol2 - 1.92476e-8 * mol3 - 1.70717e-9 * mol4
           + 1.0539e-10 * mol5;

  /**
   * The temperature of pure water at the same pressure as the brine is now calculated.
   * Note that this correlation uses temperature in K, not C
   */
  Real th20 = std::exp(std::log(tk) / (a + b * tk)) - _t_c2k;

  /**
   * The brine vapour pressure is then found by evaluating the saturation pressure for pure water
   * using this effective temperature
   */
  return PorousFlowWaterProperties::pSat(th20);
}

Real
dDensity_dP(Real pressure, Real temperature, Real /* xnacl */)
{
  //TODO: implement brine density derivative wrt P
  return PorousFlowWaterProperties::dDensity_dP(pressure, temperature);
}

Real
dViscosity_dDensity(Real pressure, Real temperature, Real /* density */, Real xnacl)
{
  /// Correlation requires molar concentration (mol/kg)
  Real mol = xnacl / ((1.0 - xnacl) * _Mnacl);
  Real mol2 = mol * mol;
  Real mol3 = mol2 * mol;

  Real a = 1.0 + 0.0816 * mol + 0.0122 * mol2 + 0.128e-3 * mol3 + 0.629e-3 * temperature
           * (1.0 - std::exp(-0.7 * mol));

  /**
   * The derivative of brine viscosity wrt density is then given by a multiplied by the
   * derivative of viscosity of pure water wrt density
   */
  Real water_density = PorousFlowWaterProperties::density(pressure, temperature);
  return a * PorousFlowWaterProperties::dViscosity_dDensity(temperature, water_density);
}
}
