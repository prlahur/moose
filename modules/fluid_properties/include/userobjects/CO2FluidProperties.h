/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef CO2FLUIDPROPERTIES_H
#define CO2FLUIDPROPERTIES_H

#include "SinglePhaseFluidPropertiesPT.h"

class CO2FluidProperties;

template<>
InputParameters validParams<CO2FluidProperties>();

/**
 * CO2 fluid properties
 * From "A New Equation of State for Carbon Dioxide Covering the Fluid Region
 * from the Triple-Point Temperature to 1100K at Pressures up to 800 MPa",
 * R. Span and W. Wagner, J. Phys. Chem. Ref. Data, volume 25, no. 6, 1996,
 * pp. 1509-1596
 */
class CO2FluidProperties : public SinglePhaseFluidPropertiesPT

{
public:
  CO2FluidProperties(const InputParameters & parameters);
  virtual ~CO2FluidProperties();

  /**
   * CO2 gas density as a function of pressure and temperature
   * (assuming an ideal gas).
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (K)
   * @return density (kg/m^3)
   */
  virtual Real rho(Real pressure, Real temperature) const override;

  /**
   * CO2 gas density as a function of pressure and temperature, and
   * derivatives wrt pressure and temperature (assuming an ideal gas).
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (K)
   * @param[out] rho density (kg/m^3)
   * @param[out] drho_dp derivative of density wrt pressure
   * @param[out] drho_dT derivative of density wrt temperature
   */
  virtual void rho_dpT(Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const override;

  /**
   * CO2 viscosity as a function of  density and temperature.
   * From Fenghour et al., The viscosity of carbon dioxide, J. Phys. Chem. Ref.
   * Data, 27, 31-44 (1998)
   * Note: critical enhancement not included
   * Valid for 217 K < T < 1000K and rho < 1400 kg/m^3
   *
   * @param density fluid density (kg/m^3)
   * @param temperature fluid temperature (K)
   * @return viscosity (Pa.s)
   */
  virtual Real mu(Real density, Real temperature) const override;

  /**
   * CO2 viscosity as a function of  density and temperature.
   * From Fenghour et al., The viscosity of carbon dioxide, J. Phys. Chem. Ref.
   * Data, 27, 31-44 (1998)
   * Note: critical enhancement not included
   * Valid for 217 K < T < 1000K and rho < 1400 kg/m^3
   *
   * @param density fluid density (kg/m^3)
   * @param temperature fluid temperature (K)
   * @param[out] mu viscosity (Pa.s)
   * @param[out] dmu_drho derivative of viscosity wrt density
   * @param[out] dmu_dT derivative of viscosity wrt temperature
   */
  virtual void mu_drhoT(Real density, Real temperature, Real & mu, Real & dmu_drho, Real & dmu_dT) const override;

  /**
   * CO2 molar mass
   * @return molar mass (kg/mol)
   */
  Real molarMass() const;

  /**
   * CO2 critical pressure
   * @return critical pressure (Pa)
   */
  Real criticalPressure() const;

  /**
   * CO2 critical temperature
   * @return critical temperature (K)
   */
  Real criticalTemperature() const;

  /**
   * CO2 critical density
   * @return critical density (kg/m^3)
   */
  Real criticalDensity() const;

  /**
   * CO2 triple point pressure
   * @return triple point pressure (Pa)
   */
  Real triplePointPressure() const;

  /**
   * CO2 triple point temperature
   * @return triple point temperature (K)
   */
  Real triplePointTemperature() const;

  /**
   * Span and Wagner Equation Of State for CO2 (reference at top of header)
   *
   * @param density CO2 density (kg/m^3)
   * @param temperature CO2 temperature (K)
   * @param[out] pressure (Pa)
   * @param[out] enthalpy (kJ/kg)
   * @param[out] internal energy (kJ/kg)
   * @param[out] isochoric heat capacity (kJ/kg)
   * @param all if false, only calculate presure, otherwise calculate all properties
   */
  void eosSW(Real density, Real temperature, Real & pressure, Real & enthalpy, Real & internal_energy, Real & cv, bool all) const;

  /**
   * Melting pressure as a function of temperature. Used to delineate solid phase
   * from the liquid phase. Valid for temperatures greater than the triple point
   * temperature.
   *
   * Eq. 3.10, from Span and Wagner (reference at top of header)
   *
   * @param temperature CO2 temperature (K)
   * @return melting pressure (Pa)
   */
  Real meltingPressure(Real temperature) const;

  /**
   * Sublimation pressure as a function of temperature. Used to delineate solid phase
   * from the gas phase. Valid for temperatures less than the triple point
   * temperature.
   *
   * Eq. 3.12, from Span and Wagner (reference at top of header)
   *
   * @param temperature CO2 temperature (K)
   * @return sublimation pressure (Pa)
   */
  Real sublimationPressure(Real temperature) const;

  /**
   * Vapour pressure as a function of temperature. Used to delineate liquid phase
   * from the gas phase. Valid for temperatures between the triple point temperature
   * and the critical temperature.
   *
   * Eq. 3.13, from Span and Wagner (reference at top of header)
   *
   * @param temperature CO2 temperature (K)
   * @return vapour pressure (Pa)
   */
  Real vapourPressure(Real temperature) const;

  /**
   * Saturated liquid density of CO2. Valid for temperatures between the triple
   * point temperature and critical temperature.
   *
   * Eq. 3.14, from Span and Wagner (reference at top of header)
   *
   * @param temperature CO2 temperature (K)
   * @return saturated liquid density (kg/m^3)
   */
  Real saturatedLiquidDensity(Real temperature) const;

  /**
   * Saturated vapour density of CO2. Valid for temperatures between the triple
   * point temperature and critical temperature.
   *
   * Eq. 3.15, from Span and Wagner (reference at top of header)
   *
   * @param temperature CO2 temperature (K)
   * @return saturated vapour density (kg/m^3)
   */
  Real saturatedVapourDensity(Real temperature) const;

  /**
   * Calculates pressure as a function of density and temperature using the
   * Span and Wagner EOS.
   *
   * @param density CO2 density (kg/m^3)
   * @param temperature CO2 temperature (K)
   * @return CO2 pressure (Pa)
   */
  Real pressureEOS(Real density, Real temperature) const;

  /**
   * Calculate density, enthalpy and internal energy of CO2 for a given pressure and
   * temperature using the Span and Wagner EOS (see documentation for eosSW()).
   * Uses Brent's method to determine density, then uses that density to calculate
   * enthalpy and internal energy.
   *
   * Valid for 0 K < temperature < 1100 K, P > 0 Pa
   *
   * @param pressure CO2 pressure (Pa)
   * @param temperature CO2 temperature (K)
   * @param[out] density (kg/m^3)
   * @param[out] enthalpy (kJ/kg)
   * @param[out] internal energy (kJ/kg)
   * @param[out] isochoric heat capacity (kJ/kg)
   */
  void eosSWProperties(Real pressure, Real temperature, Real & density, Real & enthalpy, Real & internal_energy, Real & cv) const;

  /**
   * Henry's law constant coefficients for dissolution of CO2 into water.
   * From Guidelines on the Henry's constant and vapour
   * liquid distribution constant for gases in H20 and D20 at high
   * temperatures, IAPWS (2004).
   *
   * @return constants for Henry's constant (-)
   */
  std::vector<Real> henryConstants() const;

  /**
   * Partial density of dissolved CO2. From Garcia, Density of aqueous
   * solutions of CO2, LBNL-49023 (2001).
   *
   * @param temperature fluid temperature (K)
   * @return partial molar density (kg/m^3)
   */
  Real partialDensity(Real temperature) const;

  /// Internal energy from pressure and temperature (kJ/kg)
  virtual Real e(Real pressure, Real temperature) const override;
  /// Internal energy and its derivatives wrt pressure and temperature
  virtual void e_dpT(Real pressure, Real temperature, Real & e, Real & de_dp, Real & de_dT) const override;
  /// Density and internal energy from pressure and temperature and derivatives wrt pressure and temperature
  virtual void rho_e_dpT(Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT, Real & e, Real & de_dp, Real & de_dT) const override;
  /// Speed of sound (m/s)
  virtual Real c(Real pressure, Real temperature) const override;
  /// Isobaric specific heat capacity (kJ/kg/K)
  virtual Real cp(Real pressure, Real temperature) const override;
  /// Isochoric specific heat capacity (kJ/kg/K)
  virtual Real cv(Real pressure, Real temperature) const override;
  /// Thermal conductivity (W/m/K)
  virtual Real k(Real pressure, Real temperature) const override;
  /// Specific entropy (kJ/kg/K)
  virtual Real s(Real pressure, Real temperature) const override;
  /// Specific enthalpy (kJ/kg)
  virtual Real h(Real p, Real T) const override;
  /// Specific enthalpy and its derivatives
  virtual void h_dpT(Real pressure, Real temperature, Real & h, Real & dh_dp, Real & dh_dT) const override;
  /// Thermal expansion coefficient (-)
  virtual Real beta(Real pressure, Real temperature) const override;

protected:
/**
 * Constants used in the correlations. From
 * A New Equation of State for Carbon Dioxide Covering the Fluid Region
 * from the Triple-Point Temperature to 1100K at Pressures up to 800 MPa",
 * R. Span and W. Wagner, J. Phys. Chem. Ref. Data, volume 25, no. 6, 1996,
 * pp. 1509-1596
 */

/// Molar mass of CO2 (kg/mol)
const Real _Mco2 = 44.0e-3;
/// Critical pressure (Pa)
const Real _critical_pressure = 7.3773e6;
/// Critical temperature (K)
const Real _critical_temperature = 304.1282;
/// Critical density (kg/m^3)
const Real _critical_density = 467.6;
/// Triple point pressure (Pa)
const Real _triple_point_pressure = 0.51795e6;
/// Triple point temperature (K)
const Real _triple_point_temperature = 216.592;
/// Universal gas constant (J/mol/K)
const Real _R = 8.31451;
};

#endif /* CO2FLUIDPROPERTIES_H */
