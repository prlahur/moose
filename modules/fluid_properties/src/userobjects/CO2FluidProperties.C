/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CO2FluidProperties.h"
#include "BrentsMethod.h"

template<>
InputParameters validParams<CO2FluidProperties>()
{
  InputParameters params = validParams<SinglePhaseFluidProperties>();
  params.addClassDescription("Fluid properties for methane (CH4)");
  return params;
}

CO2FluidProperties::CO2FluidProperties(const InputParameters & parameters) :
    SinglePhaseFluidProperties(parameters)
{
}

CO2FluidProperties::~CO2FluidProperties()
{
}

Real
CO2FluidProperties::molarMass()
{
  return _Mco2; //kg/mol
}

Real
CO2FluidProperties::criticalPressure()
{
  return _critical_pressure;
}

Real
CO2FluidProperties::criticalTemperature()
{
  return _critical_temperature;
}

Real
CO2FluidProperties::criticalDensity()
{
  return _critical_density;
}

Real
CO2FluidProperties::meltingPressure(Real temperature)
{
  if (temperature < _triple_point_temperature)
    mooseError("Temperature is below the triple point temperature in PorousFlowCO2PropertiesSW::meltingPressure");

  /// Correlation requires temperature in K
  Real tk = temperature + _t_c2k;
  Real ttpk = _triple_point_temperature + _t_c2k;
  Real rt = tk / ttpk;

  Real pressure = _triple_point_pressure * (1.0 + 1955.5390 * (rt - 1.0) + 2055.4593 * std::pow(rt - 1.0, 2.0));

  return pressure;
}

Real
CO2FluidProperties::sublimationPressure(Real temperature)
{
  if (temperature > _triple_point_temperature)
    mooseError("Temperature is above the triple point temperature in PorousFlowCO2PropertiesSW::sublimationPressure");

  /// Correlation requires temperature in K
  Real tk = temperature + _t_c2k;
  Real ttpk = _triple_point_temperature + _t_c2k;
  Real rt = tk / ttpk;

  Real pressure = _triple_point_pressure * std::exp((- 14.740846 * (1.0 - rt) + 2.4327015 *
    std::pow(1.0 - rt, 1.9) - 5.3061778 * std::pow(1.0 - rt, 2.9)) / rt);

  return pressure;
}

Real
CO2FluidProperties::vapourPressure(Real temperature)
{
  if (temperature < _triple_point_temperature || temperature > _critical_temperature)
    mooseError("Temperature is out of range in PorousFlowCO2PropertiesSW::vapourPressure");

  /// Correlation requires temperature in K
  Real tk = temperature + _t_c2k;
  Real tck = _critical_temperature + _t_c2k;
  Real rt = tk / tck;

  Real pressure = _critical_pressure * std::exp((-7.0602087 * (1.0 - rt) + 1.9391218 * std::pow(1.0 - rt, 1.5)
    - 1.6463597 * std::pow(1.0 - rt, 2.0) - 3.2995634 * std::pow(1.0 - rt, 4.0)) / rt);

  return pressure;
}

Real
CO2FluidProperties::saturatedLiquidDensity(Real temperature)
{
  if (temperature < _triple_point_temperature || temperature > _critical_temperature)
    mooseError("Temperature is out of range in PorousFlowCO2PropertiesSW::saturatedLiquiDensity");

  /// Correlation requires temperature in K
  Real tk = temperature + _t_c2k;
  Real tck = _critical_temperature + _t_c2k;
  Real rt = tk / tck;

  Real density = _critical_density * std::exp(1.9245108 * std::pow(1.0 - rt, 0.34)
    - 0.62385555 * std::pow(1.0 - rt, 0.5) - 0.32731127 * std::pow(1.0 - rt, 10.0 / 6.0)
    + 0.39245142 * std::pow(1.0 - rt, 11.0 / 6.0));

  return density;
}

Real
CO2FluidProperties::saturatedVapourDensity(Real temperature)
{
  if (temperature < _triple_point_temperature || temperature > _critical_temperature)
    mooseError("Temperature is out of range in PorousFlowCO2PropertiesSW::saturatedVapourDensity");

  /// Correlation requires temperature in K
  Real tk = temperature + _t_c2k;
  Real tck = _critical_temperature + _t_c2k;
  Real rt = tk / tck;

  Real density = _critical_density * std::exp(- 1.7074879 * std::pow(1.0 - rt, 0.34)
    - 0.82274670 * std::pow(1.0 - rt, 0.5) - 4.6008549 * (1.0 - rt) - 10.111178
    * std::pow(1.0 - rt, 7.0 / 3.0) - 29.742252 * std::pow(1.0 - rt, 14.0 / 3.0));

  return density;
}

/**
 * The eosSW function calculates the thermophysical properties of CO2 using a Helmholtz
 * formulation. Each of the properties is given by the free energy phi = phi0 + phir,
 * and their derivatives wrt tau (T_c/T) and delta (rho/rho_c)
 * The relations for the properties are:
 * Pressure: p = rho R T (1 + delta * dphir_dd)
 * Internal energy: u = R T tau (dphi0_dt + dphir_dt)
 * Enthalpy: h = R T (1 + tau (dphi0_dt + dphir_dt) + delta * dphir_dd)
 * Isochoric heat capacity: cv = - R tau^2 (d2phi0_dt2 + d2phir_dt2)
 * Isobaric heat capacity: cp = R (- tau^2 (d2phi0_dt2 + d2phir_dt2) + (1 + delta dphir_dd - delta tau d2phir_ddt)^2
 *                              / (1 + 2 delta dphir_dd + delta^2 d2phir_dd2)
 * Speed of sound: s^2 = R T (1 + 2 delta dphir_dd + delta^2 d2phir_dd2)   - (1 + delta dphir_dd
 *                       - delta tau d2phir_ddt)^2 / (tau^2 (d2phi0_dt2 + d2phir_dt2))
 */
void
CO2FluidProperties::eosSW(Real density, Real temperature, Real & pressure, Real & enthalpy, Real & internal_energy, Real & cv, bool all)
{
  /// Temperature in K
  Real tk = temperature + _t_c2k;

  /// Check the validity of the inputs
  if (tk <= 0.0 || tk > 1100.0 || density <= 0.0)
   mooseError("Temperature or density out of range in PorousFlowCO2PropertiesSW::eosSW");

  /// Scale the input density and temperature (in K)
  Real delta = density / _critical_density;
  Real tau = (_critical_temperature + _t_c2k) / tk;

  /// Local variables used in the calculation
  Real Psi, sum, theta, Delta, Psi_dt, Delta_dt;

  /**
   * Next set the various coefficients used in the formulation. Since C++ arrays
   * start with element 0, in order to retain the original notation the arrays here
   * are declared with one more element, and the first element is set to zero
   */

  Real n[43] = {0.0, 0.38856823203161, 2.9385475942740, -5.5867188534934,
    -0.76753199592477, 0.31729005580416, 0.54803315897767, 0.12279411220335,
    2.1658961543220, 1.5841735109724, -0.23132705405503, 0.058116916431436,
    -0.55369137205382, 0.48946615909422, -0.024275739843501, 0.062494790501678,
    -0.12175860225246, -0.37055685270086, -0.016775879700426,-0.11960736637987,
    -0.045619362508778, 0.035612789270346,-0.0074427727132052,-0.0017395704902432,
    -0.021810121289527, 0.024332166559236,-0.037440133423463, 0.14338715756878,
    -0.13491969083286, -0.023151225053480, 0.012363125492901, 0.0021058321972940,
    -0.00033958519026368, 0.0055993651771592, -0.00030335118055646, -213.65488688320,
    26641.569149272, -24027.212204557, -283.41603423999, 212.47284400179,
    -0.66642276540751, 0.72608632349897, 0.055068668612842};

  Real d[40] = {0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 2.0, 4.0, 5.0,
    5.0, 5.0, 6.0, 6.0, 6.0, 1.0, 1.0, 4.0, 4.0, 4.0, 7.0, 8.0, 2.0, 3.0, 3.0, 5.0,
    5.0, 6.0, 7.0, 8.0, 10.0, 4.0, 8.0, 2.0, 2.0, 2.0, 3.0, 3.0};

  Real t[40] = {0.0, 0.0, 0.75, 1.0, 2.0, 0.75, 2.0, 0.75, 1.5, 1.5, 2.5,
    0.0, 1.5, 2.0, 0.0, 1.0, 2.0, 3.0, 6.0, 3.0, 6.0, 8.0, 6.0, 0.0, 7.0, 12.0,
    16.0, 22.0, 24.0, 16.0, 24.0, 8.0, 2.0, 28.0, 14.0, 1.0, 0.0, 1.0, 3.0, 3.0};

  Real c[35] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0,
    4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 6.0};

  Real a0[9] = {0.0, 8.37304456, -3.70454304, 2.5, 1.99427042, 0.62105248,
    0.41195293, 1.04028922, 0.08327678};

  Real theta0[9]= {0.0, 0.0, 0.0, 0.0, 3.15163, 6.11190, 6.77708, 11.32384,
    27.08792};

  /**
   * The following are mainly empty arrays (apart from certain values) that are used to make
   * the summations easier
   */
  Real alpha[40] = {0.0};
  Real beta[43] = {0.0};
  Real gamma[40] = {0.0};
  Real eps[40] = {0.0};
  Real a[43] = {0.0};
  Real b[43] = {0.0};
  Real A[43] = {0.0};
  Real B[43] = {0.0};
  Real C[43] = {0.0};
  Real D[43] = {0.0};

  alpha[35] = 25.0;
  alpha[36] = 25.0;
  alpha[37] = 25.0;
  alpha[38] = 15.0;
  alpha[39] = 20.0;

  beta[35] = 325.0;
  beta[36] = 300.0;
  beta[37] = 300.0;
  beta[38] = 275.0;
  beta[39] = 275.0;
  beta[40] = 0.3;
  beta[41] = 0.3;
  beta[42] = 0.3;

  gamma[35] = 1.16;
  gamma[36] = 1.19;
  gamma[37] = 1.19;
  gamma[38] = 1.25;
  gamma[39] = 1.22;

  eps[35] = 1.0;
  eps[36] = 1.0;
  eps[37] = 1.0;
  eps[38] = 1.0;
  eps[39] = 1.0;

  a[40] = 3.5;
  a[41] = 3.5;
  a[42] = 3.5;

  b[40] = 0.875;
  b[41] = 0.925;
  b[42] = 0.875;

  A[40] = 0.7;
  A[41] = 0.7;
  A[42] = 0.7;

  B[40] = 0.3;
  B[41] = 0.3;
  B[42] = 1.0;

  C[40] = 10.0;
  C[41] = 10.0;
  C[42] = 12.5;

  D[40] = 275.0;
  D[41] = 275.0;
  D[42] = 275.0;

  /**
   * When the boolean flag 'all' is false, only the pressure is calculated. This is
   * used when iteratively calculating the density for a given pressure and temperature.
   * Therefore, only the terms required to calculate the pressure are calculated when the
   * flag is false.
   *
   */

  /// Only dphir_dd is required to compute the pressure
  sum = 0.0;
  for(unsigned int i = 1; i < 8; ++i)
    sum += n[i] * d[i] * std::pow(delta, d[i] - 1.0) * std::pow(tau, t[i]);

  for (unsigned int i = 8; i < 35; ++i)
    sum += n[i] * std::exp(- std::pow(delta, c[i])) * (std::pow(delta, d[i] - 1.0)
      * std::pow(tau, t[i]) * (d[i] - c[i] * std::pow(delta, c[i])));

  for(unsigned int i = 35; i < 40; ++i)
    sum += n[i] * std::pow(delta, d[i]) * std::pow(tau, t[i])
      * std::exp(- alpha[i] * std::pow(delta - eps[i], 2.0) - beta[i] * std::pow(tau - gamma[i], 2.0))
      * (d[i] / delta - 2.0 * alpha[i] * (delta - eps[i]));

  for(unsigned int i = 40; i < 43 ; ++i)
  {
    theta = 1.0 - tau + A[i] * std::pow(std::pow(delta - 1.0, 2.0), 1.0 / (2.0 * beta[i]));
    Delta = std::pow(theta, 2.0) + B[i] * std::pow(std::pow(delta - 1.0, 2.0), a[i]);
    Psi = std::exp(- C[i] * std::pow(delta - 1.0, 2.0) - D[i] * std::pow(tau - 1.0, 2.0));
    Psi_dt = - 2.0 * C[i] * (delta - 1.0) * Psi;
    Delta_dt = (delta - 1.0) * (A[i] * theta * 2.0 / beta[i]*
      std::pow(std::pow(delta - 1.0, 2.0), 1.0 / (2.0 * beta[i]) - 1.0) + 2.0 * B[i] * a[i]
      * std::pow(std::pow(delta - 1.0, 2.0), a[i] - 1.0));

    sum += n[i] * (std::pow(Delta, b[i]) * (Psi + delta * Psi_dt) + b[i]
      * std::pow(Delta, b[i] - 1.0) * Delta_dt * delta * Psi);
  }
  Real dphir_dd = sum;

  /// pressure
  pressure = _R * tk * density * (1.0 + delta * dphir_dd) / _Mco2;

  /// Calculate all other properties if 'all' is true
  if (all)
  {
    /**
     * To calculate all of the other properties, we need to first compute
     * all of the remaining derivatives of the Helmholtz free energy wrt tau and delta
     * (noting that we have already calculated dphir_dd).
     *
     * To begin, calculate the derivatives of the ideal part of the free energy psi0
     */

    /// dphi0_dt
    sum = 0.0;
    for(unsigned int i = 4; i < 9; ++i)
      sum += a0[i] * theta0[i] * (1.0 / (1.0 - std::exp(- theta0[i] * tau)) - 1.0);

   Real dphi0_dt = a0[2] + a0[3] / tau + sum;

   /// d2phi0_dt2
   sum = 0.0;
   for(unsigned int i = 4; i < 9; ++i)
     sum += a0[i] * theta0[i] * theta0[i] * std::exp(- theta0[i] * tau) *
       std::pow(1.0 - std::exp(- theta0[i] * tau), - 2.0);

   Real d2phi0_dt2 = - a0[3] / tau / tau + sum;

   /**
    * Now calculate the remaining required derivatives of the residual part
    * of the free energy psir
    */

   /// dphir_dt
   sum = 0.0;
   for (unsigned int i = 1; i < 8; ++i)
     sum += n[i] * t[i] * std::pow(delta, d[i]) * std::pow(tau, t[i] - 1.0);

   for (unsigned int i = 8; i < 35; ++i)
     sum += n[i] * std::exp( - std::pow(delta, c[i])) * std::pow(delta, d[i])
       * t[i] * std::pow(tau, t[i] - 1.0);

   for(unsigned int i = 35; i < 40; ++i)
     sum += n[i] * std::pow(delta, d[i]) * std::pow(tau, t[i]) * std::exp(- alpha[i]
       * std::pow(delta - eps[i], 2) - beta[i] * std::pow(tau - gamma[i], 2)) * (t[i] / tau
       - 2.0 * beta[i] * (tau - gamma[i]));

   for(unsigned int i = 40; i < 43; ++i)
   {
     theta = 1.0 - tau + A[i] * std::pow(std::pow(delta - 1.0, 2.0), 1.0 / (2.0 * beta[i]));
     Delta = std::pow(theta, 2) + B[i] * std::pow(std::pow(delta - 1.0, 2.0), a[i]);
     Psi = std::exp(- C[i] * std::pow(delta - 1.0, 2.0) - D[i] * std::pow(tau - 1.0, 2.0));

     sum += n[i] * delta * Psi * (- 2.0 * theta * b[i] * std::pow(Delta, b[i] - 1.0)
       - 2.0 * std::pow(Delta, b[i]) * D[i] * (tau - 1.0));
   }
   Real dphir_dt = sum;

   /// Calculate d2phir_dt2
   sum = 0.0;
   for (unsigned int i = 1; i < 8; ++i)
     sum += n[i] * t[i] * (t[i] - 1.0) * std::pow(delta, d[i]) * std::pow(tau, t[i] - 2.0);

   for (unsigned int i = 8; i < 35; ++i)
     sum += n[i] * t[i] * (t[i] - 1.0) * std::pow(delta, d[i]) * std::exp(- std::pow(delta, c[i]))
        * std::pow(tau, t[i] - 2.0);

   for(unsigned int i = 35; i < 40; ++i)
     sum += n[i] * std::pow(delta, d[i]) * std::pow(tau, t[i]) * std::exp(- alpha[i]
       * std::pow(delta - eps[i], 2.0) - beta[i] * std::pow(tau - gamma[i], 2.0)) * (std::pow(t[i] / tau
       - 2.0 * beta[i] * (tau - gamma[i]), 2.0) - t[i] / tau / tau - 2.0 * beta[i]);

   Real dPsi_dt, dDelta_dt, d2Delta_dt2, d2Psi_dt2;
   for(unsigned int i = 40; i < 43; ++i)
   {
     theta = 1.0 - tau + A[i] * std::pow(std::pow(delta - 1.0, 2.0), 1.0 / (2.0 * beta[i]));
     Delta = std::pow(theta, 2) + B[i] * std::pow(std::pow(delta - 1.0, 2.0), a[i]);
     Psi = std::exp(- C[i] * std::pow(delta - 1.0, 2.0) - D[i] * std::pow(tau - 1.0, 2.0));
     dDelta_dt = -2.0 * theta * b[i] * std::pow(Delta, b[i] - 1.0);
     d2Delta_dt2 = 2.0 * b[i] * std::pow(Delta, b[i] - 1.0) + 4.0 * theta * theta * b[i] *
       (b[i] - 1.0) * std::pow(Delta, b[i] - 2.0);
     dPsi_dt = - 2.0 * D[i] * (tau - 1.0) * Psi;
     d2Psi_dt2 = 2.0 * D[i] * (2.0 * D[i] * (tau - 1.0) * (tau - 1.0) - 1.0) * Psi;
     sum += n[i] * delta * (Psi * d2Delta_dt2 + 2.0 * dDelta_dt * dPsi_dt + Delta * d2Psi_dt2);
   }
   Real d2phir_dt2 = sum;

   /// Calculate d2phir_dd2
   sum = 0.0;
   for (unsigned int i = 1; i < 8; ++i)
     sum += n[i] * d[i] * (d[i] - 1.0) * std::pow(delta, d[i] - 2.0) * std::pow(tau, t[i]);

   for (unsigned int i = 8; i < 35; ++i)
     sum += n[i] * std::exp(- std::pow(delta, c[i])) * (std::pow(delta, d[i] - 2.0) * std::pow(tau, t[i]) *
       ((d[i] - c[i] * std::pow(delta, c[i])) * (d[i] - 1.0 - c[i] * std::pow(delta, c[i]))) - c[i] * c[i] *
       std::pow(delta, c[i]));

   for(unsigned int i = 35; i < 40; ++i)
     sum += n[i] * std::pow(tau, t[i]) * std::exp(- alpha[i] * std::pow(delta - eps[i], 2.0) - beta[i] *
       std::pow(tau - gamma[i], 2.0)) * (- 2.0 * alpha[i] * std::pow(delta, d[i]) + 4.0 * alpha[i] *
       alpha[i] * std::pow(delta, d[i]) * std::pow(delta - eps[i], 2.0) - 4.0 * d[i] * alpha[i] *
       std::pow(delta, d[i] - 1.0) * (delta - eps[i]) + d[i] * (d[i] - 1.0) * std::pow(delta, d[i] - 2.0));

   Real dDelta_dd, d2Delta_dd2, d2Psi_dd2;
   for(unsigned int i = 40; i < 43; ++i)
   {
     theta = 1.0 - tau + A[i] * std::pow(std::pow(delta - 1.0, 2.0), 1.0 / (2.0 * beta[i]));
     Delta = std::pow(theta, 2) + B[i] * std::pow(std::pow(delta - 1.0, 2.0), a[i]);
     Psi = std::exp(- C[i] * std::pow(delta - 1.0, 2.0) - D[i] * std::pow(tau - 1.0, 2.0));
     dDelta_dd = 0; /// UP to here - need to fix following part
     d2Delta_dt2 = 2.0 * b[i] * std::pow(Delta, b[i] - 1.0) + 4.0 * theta * theta * b[i] *
       (b[i] - 1.0) * std::pow(Delta, b[i] - 2.0);
     dPsi_dt = - 2.0 * D[i] * (tau - 1.0) * Psi;
     d2Psi_dt2 = 2.0 * D[i] * (2.0 * D[i] * (tau - 1.0) * (tau - 1.0) - 1.0) * Psi;
     sum += n[i] * delta * (Psi * d2Delta_dt2 + 2.0 * dDelta_dt * dPsi_dt + Delta * d2Psi_dt2);
   }
   Real d2phir_dd2 = sum;

   /// Isobaric heat capacity: cp = R (- tau^2 (d2phi0_dt2 + d2phir_dt2) + (1 + delta dphir_dd - delta tau d2phir_ddt)^2
   ///                             / (1 + 2 delta dphir_dd + delta^2 d2phir_dd2)


  /**
   * Compute the enthalpy, internal energy, isochoric and isobaric heat capacities, and the
   * speed of sound.
   * Note: the reference state for the enthalpy is 298.15 K and 0.101325 MPa.
   * Note: divide by 1000 to get kJ
   */
   enthalpy = _R * tk / _Mco2 * (tau * (dphi0_dt + dphir_dt) + 1.0 + delta * dphir_dd) / 1000.0;
   internal_energy = _R * tk / _Mco2 * tau * (dphi0_dt + dphir_dt) / 1000.0;
   cv = - _R * tau * tau * (d2phi0_dt2 + d2phir_dt2);

  }
}

Real
CO2FluidProperties::pressureEOS(Real density, Real temperature)
{
  /// Check that the input parameters are within the region of validity
  if (temperature <= - _t_c2k || temperature > 1100.0 || density <= 0.0)
    mooseError("Parameters out of range in PorousFlowCO2PropertiesSW::pressure");

  Real pressure = 0.0;
  Real enthalpy = 0.0;
  Real internal_energy = 0.0;
  Real cv = 0.0;

  if (temperature > _triple_point_temperature && temperature < _critical_temperature)
  {
    Real gas_density = saturatedVapourDensity(temperature);
    Real liquid_density = saturatedLiquidDensity(temperature);

    if (density < gas_density || density > liquid_density)
      eosSW(density, temperature, pressure, enthalpy, internal_energy, cv, false);
    else
      pressure = vapourPressure(temperature);
    }
  else
    eosSW(density, temperature, pressure, enthalpy, internal_energy, cv, false);

  return pressure;
}

Real
CO2FluidProperties::pressureDifference(Real density, Real temperature, Real pressure)
{
  return pressureEOS(density, temperature) - pressure;
}

void
CO2FluidProperties::eosSWProperties(Real pressure, Real temperature, Real & density, Real & enthalpy, Real & internal_energy, Real & cv)
{
  /// Check that the pressure and temperature are within the valid range
  if (pressure <= 0.0)
    mooseError("Input pressure in PorousFlowCO2PropertiesSW::eosSWProperties must be greater than 0");

  if (temperature < - _t_c2k || temperature > 1100.0 + _t_c2k)
    mooseError("Input temperature in PorousFlowCO2PropertiesSW::eosSWProperties must be between (-273.15 C < temperature < 1375.15 C)");

  /// Check that the pressure and temperature are not in the solid phase region
  if(((temperature > _triple_point_temperature) && (pressure > meltingPressure(temperature)))
  || ((temperature < _triple_point_temperature) && (pressure > sublimationPressure(temperature))))
    mooseError("Input pressure and temperature in PorousFlowCO2PropertiesSW::eosSWProperties correspond to solid CO2 phase");

  /// Determine a bracketing interval for the density used in the root finding algortihm
  Real lower_density = 1.0;
  Real upper_density = 1000.0;

  BrentsMethod::bracket(pressureDifference, lower_density, upper_density, temperature, pressure);

  /// Now find the density using Brent's method
  Real eps = 1.0e-12;
  density = BrentsMethod::root(pressureDifference, lower_density, upper_density, temperature, pressure, eps);

  /// Using this density, calculate all other properties
  eosSW(density, temperature, pressure, enthalpy, internal_energy, cv, true);
}

Real
CO2FluidProperties::rho(Real pressure, Real temperature) const
{
  return 0.0; ///FIXME
}

void
CO2FluidProperties::rho_dpT(Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  Real tk = temperature + _t_c2k;
  ///FIXME
}

Real
CO2FluidProperties::mu(Real /*pressure*/, Real temperature) const
{
 return 0.; ///FIXME
}

Real
CO2FluidProperties::viscosity(Real temperature, Real density)
{
  Real tk = temperature + _t_c2k;
  Real tstar = tk / 251.196;
  Real a[5] = {0.235156, -0.491266, 5.211155e-2, 5.347906e-2, -1.537102e-2};
  Real d[5] = {0.4071119e-2, 0.7198037e-4, 0.2411697e-16, 0.2971072e-22, -0.1627888e-22};
  unsigned int j[5] = {1, 1, 4, 1, 2};
  unsigned int i[5] = {1, 2, 6, 8, 8};

  // Zero-denisty viscosity
  Real sum = 0.0;

  for (unsigned int n = 0; n < 5; ++n)
    sum += a[n] * std::pow(std::log(tstar), n);

  Real theta = std::exp(sum);
  Real mu0 = 1.00697 * std::sqrt(tk) / theta;

  Real b[5];
  for (unsigned int n = 0; n < 5; ++n)
    b[n] = d[n] /  std::pow(tstar, j[n] - 1);

  // Excess viscosity due to density
  Real mu = 0.0;

  for (unsigned int n = 0; n < 5; ++n)
    mu += b[n] * std::pow(density, i[n]);

  return (mu0 + mu) * 1e-6; // convert to Pa.s
}

Real
CO2FluidProperties::partialDensity(Real temperature)
{
  Real t2 = temperature * temperature;
  Real t3 = t2 * temperature;

  Real V = 37.51 - 9.585e-2 * temperature + 8.74e-4 * t2 - 5.044e-7 * t3;

  return 1.e6 * _Mco2 / V;
}

Real
CO2FluidProperties::dViscosity_dDensity(Real temperature, Real density)
{
  Real tk = temperature + _t_c2k;
  Real tstar = tk / 251.196;
  Real d[5] = {0.4071119e-2, 0.7198037e-4, 0.2411697e-16, 0.2971072e-22, -0.1627888e-22};
  unsigned int j[5] = {1, 1, 4, 1, 2};
  unsigned int i[5] = {1, 2, 6, 8, 8};

  Real b[5];
  for (unsigned int n = 0; n < 5; ++n)
    b[n] = d[n] /  std::pow(tstar, j[n] - 1);

  // Excess viscosity due to density
  Real dmu = 0.0;

  for (unsigned int n = 0; n < 5; ++n)
    dmu += b[n] * i[n] * std::pow(density, i[n] - 1);

  return  dmu * 1e-6; // convert to Pa.s
}

std::vector<Real>
CO2FluidProperties::henryConstants()
{
  std::vector<Real> co2henry;
  co2henry.push_back(-8.55445);
  co2henry.push_back(4.01195);
  co2henry.push_back(9.52345);

  return co2henry;
}

void
CO2FluidProperties::mu_dpT(Real pressure, Real temperature, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  ///FIXME
  const Real a[6] = {2.968267e-1, 3.711201e-2, 1.218298e-5, -7.02426e-8, 7.543269e-11,
    -2.7237166e-14};

  const Real tk = temperature + _t_c2k;
  const Real tk2 = tk * tk;
  const Real tk3 = tk2 * tk;
  const Real tk4 = tk3 * tk;

  mu = this->mu(pressure, temperature);
  dmu_dp = 0.0;
  dmu_dT = a[1] + 2.0 * a[2] * tk + 3.0 * a[3] * tk2 + 4.0 * a[4] * tk3 + 5.0 * a[5] * tk4;
}
