/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWBRINEPROPERTIES_H
#define POROUSFLOWBRINEPROPERTIES_H

#include "MooseTypes.h"
#include "MooseError.h"

namespace PorousFlowBrineProperties
{
  /**
   * Fluid name
   * @return fluid Name
   */
  std::string fluidName();

  /**
   * NaCl molar mass.
   * @return molar mass (kg/mol)
   */
  Real molarMass();

  /**
   * Water molar mass.
   * @return molar mass of water (kg/mol)
   */
  Real molarMassH2O();

  /**
   * Density of brine.
   * From Driesner, The system H2o-NaCl. Part II: Correlations for molar volume,
   * enthalpy, and isobaric heat capacity from 0 to 1000 C, 1 to 500 bar, and 0
   * to 1 Xnacl, Geochimica et Cosmochimica Acta 71, 4902-4919 (2007).
   *
   * @param pressure brine pressure (Pa)
   * @param temperature brine temperature (C)
   * @param xnacl salt mass fraction (-)
   * @return brine density (kg/m^3)
   */
  Real density(Real pressure, Real temperature, Real xnacl);

  /**
   * Viscosity of brine.
   * From Phillips et al, A technical databook for geothermal energy utilization,
   * LbL-12810 (1981).
   *
   * @param pressure brine pressure (Pa)
   * @param temperature brine temperature (C)
   * @param xnacl salt mass fraction (-)
   * @return brine viscosity (Pa.s)
   */
  Real viscosity(Real pressure, Real temperature, Real xnacl);

  /**
   * Density of halite (solid NaCl)
   * From Driesner, The system H2o-NaCl. Part II: Correlations for molar volume,
   * enthalpy, and isobaric heat capacity from 0 to 1000 C, 1 to 500 bar, and 0
   * to 1 Xnacl, Geochimica et Cosmochimica Acta 71, 4902-4919 (2007).
   *
   * @param pressure pressure (Pa)
   * @param temperature halite temperature (C)
   * @return density (kg/m^3)
   */
  Real haliteDensity(Real pressure, Real temperature);

  /**
   * Halite solubility
   * Originally from Potter et al., A new method for determining the solubility
   * of salts in aqueous solutions at elevated temperatures, J. Res. U.S. Geol.
   * Surv., 5, 389-395 (1977). Equation describing halite solubility is repeated
   * in Chou, Phase relations in the system NaCI-KCI-H20. III: Solubilities of
   * halite in vapor-saturated liquids above 445°C and redetermination of phase
   * equilibrium properties in the system NaCI-HzO to 1000°C and 1500 bars,
   * Geochimica et Cosmochimica Acta 51, 1965-1975 (1987).
   *
   * @param temperature temperature (C)
   * @return halite solubility (kg/kg)
   *
   * This correlation is valid for 0 <= T << 424.5 C
   */
  Real haliteSolubility(Real temperature);

  /**
   * Brine vapour pressure
   * From Haas, Physical properties of the coexisting phases and thermochemical
   * properties of the H20 component in boiling NaCl solutions, Geological Survey
   * Bulletin, 1421-A (1976).
   *
   * @param temperature brine temperature (C)
   * @param xnacl salt mass fraction (-)
   * @return brine vapour pressure (Pa)
   */
  Real pSat(Real temperature, Real xnacl);

  /**
   * Derivative of brine density with respect to presure.
   *
   * @param pressure brine pressure (Pa)
   * @param temperature brine temperature (C)
   * @param xnacl salt mass fraction (-)
   * @return derivative of brine density wrt pressure (kg/m^3/Pa)
   */
  Real dDensity_dP(Real pressure, Real temperature, Real xnacl);

  /**
   * Derivative of brine viscosity with respect to density
   *
   * @param pressure brine pressure (Pa)
   * @param temperature brine temperature (C)
   * @param density brine density (kg/m^3)
   * @param xnacl salt mass fraction (-)
   * @return derivative of brine viscosity wrt density
   */
  Real dViscosity_dDensity(Real pressure, Real temperature, Real density, Real xnacl);

  /// Molar mass of NaCl (kg/mol)
  const Real _Mnacl = 58.443e-3;
}

#endif // POROUSFLOWBRINEPROPERTIES_H
