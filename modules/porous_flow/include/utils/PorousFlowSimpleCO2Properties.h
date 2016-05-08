/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWSIMPLECO2PROPERTIES_H
#define POROUSFLOWSIMPLECO2PROPERTIES_H

#include "MooseTypes.h"
#include "MooseError.h"

namespace PorousFlowSimpleCO2Properties
{
  /**
   * Fluid name
   * @return fluid Name
   */
  std::string fluidName();

  /**
   * CO2 molar mass.
   * @return molar mass (kg/mol)
   */
  Real molarMass();

  /**
   * CO2 critical pressure.
   * @return critical pressure (Pa)
   */
  Real criticalPressure();

  /**
   * CO2 critical temperature.
   * @return critical temperature (C)
   */
  Real criticalTemperature();

  /**
   * CO2 density as a function of pressure and temperature.
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return density (kg/m^3)
   */
  Real density(Real pressure, Real temperature);

  /**
   * CO2 viscosity as a function of  pressure and temperature.
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @param density phase density (kg/m^3)
   * @return viscosity (Pa.s)
   */
  Real viscosity(Real pressure, Real temperature, Real density);

  /**
   * Henry's law constant coefficients for dissolution of CO2 into water.
   * From Guidelines on the Henry's constant and vapour
   * liquid distribution constant for gases in H20 and D20 at high
   * temperatures, IAPWS (2004).
   *
   * @return constants for Henry's constant (-)
   */
  std::vector<Real> henryConstants();

  /**
   * Partial density of dissolved CO2. From Garcia, Density of aqueous
   * solutions of CO2, LBNL-49023 (2001).
   *
   * @param temperature fluid temperature (C)
   * @return partial molar density (kg/m^3)
   */
  Real partialDensity(Real temperature);

  /**
   * Density of supercritical CO2. From Ouyang, New correlations for predicting
   * the density and viscosity of supercritical Carbon Dioxide under conditions
   * expected in Carbon Capture and Sequestration operations, The Open Petroleum
   * Engineering Journal, 4, 13-21 (2011)
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return density (kg/m^3)
   */
  Real supercriticalDensity(Real pressure, Real temperature);

  /**
   * CO2 gas density as a function of  pressure and temperature.
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return density (kg/m^3)
   */
  Real gasDensity(Real pressure, Real temperature);

  /**
   * CO2 gas viscosity as a function of temperature and density.
   * From Fenghour et al., The viscosity of carbon dioxide, J. Phys. Chem. Ref.
   * Data, 27, 31-44 (1998).
   * Note: used up to the critical point. The supercritical viscosity is calculated
   * in supercriticalViscosity(). As a result, the critical enhancement is not
   * implemented.
   *
   * @param temperature fluid temperature (C)
   * @param density fluid density (kg/m^3)
   * @return viscosity (Pa.s)
   */
  Real gasViscosity(Real temperature, Real density);

  /**
   * Viscosity of supercritical CO2. From Ouyang, New correlations for predicting
   * the density and viscosity of supercritical Carbon Dioxide under conditions
   * expected in Carbon Capture and Sequestration operations, The Open Petroleum
   * Engineering Journal, 4, 13-21 (2011)
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return viscosity (Pa.s)
   */
  Real supercriticalViscosity(Real pressure, Real temperature);

  /**
   * Derivative of the density of CO2 as a function of
   * pressure.
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return derivative of CO2 density (kg/m^3) with respect to pressure (Pa)
   */
  Real dDensity_dP(Real pressure, Real temperature);

  /**
   * Derivative of the density of CO2 as a function of
   * temperature.
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return derivative of CO2 density (kg/m^3) with respect to temperature (C)
   */
  Real dDensity_dT(Real pressure, Real temperature);

  /**
   * Derivative of the density of gaseous CO2 as a function of
   * pressure.
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return derivative of CO2 density (kg/m^3) with respect to pressure (Pa)
   */
  Real dGasDensity_dP(Real pressure, Real temperature);

  /**
   * Derivative of the density of gaseous CO2 as a function of
   * temperature.
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return derivative of CO2 density (kg/m^3) with respect to temperature (C)
   */
  Real dGasDensity_dT(Real pressure, Real temperature);

  /**
   * Derivative of the density of supercritical CO2 as a function of
   * pressure.
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return derivative of CO2 density (kg/m^3) with respect to pressure (Pa)
   */
  Real dSupercriticalDensity_dP(Real pressure, Real temperature);

  /**
   * Derivative of the density of supercritical CO2 as a function of
   * temperature.
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return derivative of CO2 density (kg/m^3) with respect to temperature (C)
   */
  Real dSupercriticalDensity_dT(Real pressure, Real temperature);

  /**
   * The derivative of CO2 viscosity with respect to density
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return derivative of CO2 viscosity wrt density
   */
  Real dViscosity_dDensity(Real pressure, Real temperature, Real density);

  /**
   * Derivative of the viscosity of gaseous CO2 as a function of density
   *
   * @param temperature fluid temperature (C)
   * @param density fluid density (kg/m^3)
   * @return derivative of CO2 viscosity wrt density
   */
  Real dGasViscosity_dDensity(Real temperature, Real density);

  /**
   * Derivative of the viscosity of supercritical CO2 as a function of density
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return derivative of CO2 viscosity wrt density
   */
  Real dSupercriticalViscosity_dDensity(Real pressure, Real temperature);

  /**
   * Derivative of the viscosity of supercritical CO2 as a function of pressure
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (C)
   * @return derivative of CO2 viscosity wrt pressure
   */
  Real dSupercriticalViscosity_dP(Real pressure, Real temperature);

  /**
   * Constants used in the correlations
   */

  /// Conversion of temperature from Celcius to Kelvin
  const Real _t_c2k = 273.15;
  /// Molar mass of CO2
  const Real _Mco2 = 44.0e-3;
  /// Critical pressure (Pa)
  const Real _p_critical = 7.3773e6;
  /// Critical temperature (C)
  const Real _t_critical = 30.9782;
  // Conversion from pressure in Pa to psia
  const Real _pa2psia = 1.45037738007e-4;
  /// Arrays of coefficients used in the calculation of supercritical density
  const Real _b_scd[5][5] = {{-2.148322085348e5, 1.168116599408e4, -2.302236659392e2,
    1.967428940167, -6.184842764145e-3}, {4.757146002428e2, -2.619250287624e1,
    5.215134206837e-1, -4.494511089838e-3, 1.423058795982e-5}, {-3.713900186613e-1,
    2.072488876536e-2, -4.169082831078e-4, 3.622975674137e-6, -1.155050860329e-8},
    {1.228907393482e-4, -6.930063746226e-6, 1.406317206628e-7, -1.230995287169e-9,
    3.948417428040e-12}, {-1.466408011784e-8, 8.338008651366e-10, -1.704242447194e-11,
    1.500878861807e-13, -4.838826574173e-16}};
  const Real _c_scd[5][5] = {{6.897382693936e2, 2.730479206931, -2.254102364542e-2,
    -4.651196146917e-3, 3.439702234956e-5}, {2.213692462613e-1, -6.547268255814e-3,
    5.982258882656e-5, 2.274997412526e-6, -1.888361337660e-8}, {-5.118724890479e-5,
    2.019697017603e-6, -2.311332097185e-8, -4.079557404679e-10, 3.893599641874e-12},
    {5.517971126745e-9, -2.415814703211e-10, 3.121603486524e-12, 3.171271084870e-14,
    -3.560785550401e-16}, {-2.184152941323e-13, 1.010703706059e-14, -1.406620681883e-16,
      -8.957731136447e-19, 1.215810469539e-20}};
  /// Arrays of coefficients used in the calculation of supercritical viscosity
  const Real _b_scv[5][5] = {{-1.958098980443E+01, 1.123243298270E+00, -2.320378874100E-02,
    2.067060943050E-04, -6.740205984528E-07}, {4.187280585109E-02, -2.425666731623E-03,
    5.051177210444E-05, -4.527585394282E-07, 1.483580144144E-09}, {-3.164424775231E-05,
    1.853493293079E-06, -3.892243662924E-08, 3.511599795831E-10, -1.156613338683E-12},
    {1.018084854204E-08, -6.013995738056E-10, 1.271924622771E-11, -1.154170663233E-13,
    3.819260251596E-16}, {-1.185834697489E-12, 7.052301533772E-14, -1.500321307714E-15,
    1.368104294236E-17, -4.545472651918E-20}};
  const Real _c_scv[5][5] = {{1.856798626054E-02, 3.083186834281E-03, -1.004022090988E-04,
    8.331453343531E-07, -1.824126204417E-09}, {6.519276827948E-05, -3.174897980949E-06,
    7.524167185714E-08, -6.141534284471E-10, 1.463896995503E-12}, {-1.310632653461E-08,
    7.702474418324E-10, -1.830098887313E-11, 1.530419648245E-13, -3.852361658746E-16},
    {1.335772487425E-12, -8.113168443709E-14, 1.921794651400E-15, -1.632868926659E-17,
    4.257160059035E-20}, {-5.047795395464E-17, 3.115707980951E-18, -7.370406590957E-20,
    6.333570782917E-22, -1.691344581198E-24}};
}

#endif // POROUSFLOWSIMPLECO2PROPERTIES_H
