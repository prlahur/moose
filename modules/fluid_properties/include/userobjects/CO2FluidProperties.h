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
   * Note: this is calculated iteratively
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
  Real pressure(Real density, Real temperature) const;

  /**
   * Internal function to calculate pressure as a function of density and
   * temperature using the Span and Wagner EOS. This function is called by
   * pressure(density, temperature)
   *
   * @param density CO2 density (kg/m^3)
   * @param temperature CO2 temperature (K)
   * @return CO2 pressure (Pa)
   */
  Real pressureSW(Real density, Real temperature) const;

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
   * Helmholtz free energy for CO2 from the Span and Wagner EOS
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return phi Helmholtz free energy
   */
  Real phiSW(Real delta, Real tau) const;

  /**
   * Derivative of Helmholtz free energy wrt delta
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return derivative of Helmholtz free energy wrt delta
   */
  Real dphiSW_dd(Real delta, Real tau) const;

  /**
   * Derivative of Helmholtz free energy wrt tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return derivative of Helmholtz free energy wrt tau
   */
  Real dphiSW_dt(Real delta, Real tau) const;

  /**
   * Second derivative of Helmholtz free energy wrt delta
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of Helmholtz free energy wrt delta
   */
  Real d2phiSW_dd2(Real delta, Real tau) const;

  /**
   * Second derivative of Helmholtz free energy wrt tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of Helmholtz free energy wrt tau
   */
  Real d2phiSW_dt2(Real delta, Real tau) const;

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

/// Coefficients for the ideal gas component of the Helmholtz free energy
std::vector<Real> _a0 {1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678};
std::vector<Real> _theta0 {3.15163, 6.11190, 6.77708, 11.32384, 27.08792};

/// Coefficients for the residual component of the Helmholtz free energy
std::vector<Real> _n1 {0.38856823203161, 2.9385475942740, -5.5867188534934,
  -0.76753199592477, 0.31729005580416, 0.54803315897767, 0.12279411220335};
std::vector<Real> _d1 {1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0};
std::vector<Real> _t1 {0.0, 0.75, 1.0, 2.0, 0.75, 2.0, 0.75};
std::vector<Real> _n2 {2.1658961543220, 1.5841735109724, -0.23132705405503,
   0.058116916431436, -0.55369137205382, 0.48946615909422, -0.024275739843501,
   0.062494790501678, -0.12175860225246, -0.37055685270086, -0.016775879700426,
  -0.11960736637987, -0.045619362508778, 0.035612789270346,-0.0074427727132052,
  -0.0017395704902432, -0.021810121289527, 0.024332166559236,-0.037440133423463,
   0.14338715756878, -0.13491969083286, -0.023151225053480, 0.012363125492901,
   0.0021058321972940, -0.00033958519026368, 0.0055993651771592, -0.00030335118055646};
std::vector<Real> _d2 {1.0, 2.0, 4.0, 5.0, 5.0, 5.0, 6.0, 6.0, 6.0, 1.0, 1.0, 4.0,
  4.0, 4.0, 7.0, 8.0, 2.0, 3.0, 3.0, 5.0, 5.0, 6.0, 7.0, 8.0, 10.0, 4.0, 8.0};
std::vector<Real> _t2 {1.5, 1.5, 2.5, 0.0, 1.5, 2.0, 0.0, 1.0, 2.0, 3.0, 6.0, 3.0,
  6.0, 8.0, 6.0, 0.0, 7.0, 12.0, 16.0, 22.0, 24.0, 16.0, 24.0, 8.0, 2.0, 28.0, 14.0};
std::vector<Real> _c2 {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0,
  2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 6.0};
std::vector<Real> _n3 {-213.65488688320, 26641.569149272, -24027.212204557,
  -283.41603423999, 212.47284400179};
std::vector<Real> _d3 {2.0, 2.0, 2.0, 3.0, 3.0};
std::vector<Real> _t3 {1.0, 0.0, 1.0, 3.0, 3.0};
std::vector<Real> _alpha3 {25.0, 25.0, 25.0, 15.0, 20.0};
std::vector<Real> _beta3 {325.0, 300.0, 300.0, 275.0, 275.0};
std::vector<Real> _gamma3 {1.16, 1.19, 1.19, 1.25, 1.25};
std::vector<Real> _eps3 {1.0, 1.0, 1.0, 1.0, 1.0};
std::vector<Real> _n4 {-0.66642276540751, 0.72608632349897, 0.055068668612842};
std::vector<Real> _a4 {3.5, 3.5, 3.5};
std::vector<Real> _b4 {0.875, 0.925, 0.875};
std::vector<Real> _beta4 {0.3, 0.3, 0.3};
std::vector<Real> _A4 {0.7, 0.7, 0.7};
std::vector<Real> _B4 {0.3, 0.3, 1.0};
std::vector<Real> _C4 {10.0, 10.0, 12.5};
std::vector<Real> _D4 {275.0, 275.0, 275.0};
};

#endif /* CO2FLUIDPROPERTIES_H */
