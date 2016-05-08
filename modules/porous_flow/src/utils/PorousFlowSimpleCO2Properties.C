/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowSimpleCO2Properties.h"

namespace PorousFlowSimpleCO2Properties
{

std::string
fluidName()
{
  return "CO2";
}

Real
molarMass()
{
  return _Mco2; //kg/mol
}

Real
criticalPressure()
{
  return _p_critical;
}

Real
criticalTemperature()
{
  return _t_critical;
}

Real
density(Real pressure, Real temperature)
{
  Real rho;

  if (pressure <= _p_critical)
    rho = gasDensity(pressure, temperature);
  else
    rho = supercriticalDensity(pressure, temperature);

  return rho;
}

Real
viscosity(Real pressure, Real temperature, Real density)
{
  Real mu;

  if (pressure <= _p_critical)
    mu = gasViscosity(temperature, density);
  else
    mu = supercriticalViscosity(pressure, temperature);

  return mu;
}

Real
dDensity_dP(Real pressure, Real temperature)
{
  Real drho;

  if (pressure <= _p_critical)
    drho = dGasDensity_dP(pressure, temperature);
  else
    drho = dSupercriticalDensity_dP(pressure, temperature);

  return drho;
}

Real
dDensity_dT(Real pressure, Real temperature)
{
  Real drho;

  if (pressure <= _p_critical)
    drho = dGasDensity_dT(pressure, temperature);
  else
    drho = dSupercriticalDensity_dT(pressure, temperature);

  return drho;
}

Real
gasViscosity(Real temperature, Real density)
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
supercriticalViscosity(Real pressure, Real temperature)
{

  Real a[5];

  Real t1 = temperature;
  Real t2 = t1 * t1;
  Real t3 = t2 * t1;
  Real t4 = t3 * t1;

  // Correlation uses pressure in psia
  Real p1 = pressure * _pa2psia;
  Real p2 = p1 * p1;
  Real p3 = p2 * p1;
  Real p4 = p3 * p1;

  if (p1 <= 3000)
    for (unsigned int i = 0; i < 5; ++i)
      a[i] = _b_scv[i][0] + _b_scv[i][1] * t1 + _b_scv[i][2] * t2 + _b_scv[i][3] * t3 + _b_scv[i][4] * t4;
  else
    for (unsigned int i = 0; i < 5; ++i)
      a[i] = _c_scv[i][0] + _c_scv[i][1] * t1 + _c_scv[i][2] * t2 + _c_scv[i][3] * t3 + _c_scv[i][4] * t4;

 Real mu = a[0] + a[1] * p1 + a[2] * p2 + a[3] * p3 + a[4] * p4;

 return mu * 1e-3; //cP to Pa.s
}

Real
partialDensity(Real temperature)
{
  Real t2 = temperature * temperature;
  Real t3 = t2 * temperature;

  Real V = 37.51 - 9.585e-2 * temperature + 8.74e-4 * t2 - 5.044e-7 * t3;

  return 1.e6 * _Mco2 / V;
}

Real
supercriticalDensity(Real pressure, Real temperature)
{
  Real a[5];

  Real t1 = temperature;
  Real t2 = t1 * t1;
  Real t3 = t2 * t1;
  Real t4 = t3 * t1;

  // Correlation uses pressure in psia
  Real pa2psia = 1.45037738007e-4;
  Real p1 = pressure * pa2psia;
  Real p2 = p1 * p1;
  Real p3 = p2 * p1;
  Real p4 = p3 * p1;

  if (p1 <= 3000)
    for (unsigned int i = 0; i < 5; ++i)
      a[i] = _b_scd[i][0] + _b_scd[i][1] * t1 + _b_scd[i][2] * t2 + _b_scd[i][3] * t3 + _b_scd[i][4] * t4;
  else
    for (unsigned int i = 0; i < 5; ++i)
      a[i] = _c_scd[i][0] + _c_scd[i][1] * t1 + _c_scd[i][2] * t2 + _c_scd[i][3] * t3 + _c_scd[i][4] * t4;

  return a[0] + a[1] * p1 + a[2] * p2 + a[3] * p3 + a[4] * p4;
}

Real
gasDensity(Real pressure, Real temperature)
{
  Real tk = temperature + _t_c2k;
  Real tc = std::pow((tk * 1.e-2), 10./3.);
  Real pc = pressure * 1.e-6;

  Real vc1 = 1.8882e-4 * tk;
  Real vc2 = - pc * (8.24e-2 + 1.249e-2 * pc) / tc;

  return pc / (vc1 + vc2);
}

Real
dGasDensity_dP(Real pressure, Real temperature)
{
  Real tk = temperature + _t_c2k;
  Real tc = std::pow((tk * 1.e-2), 10./3.);
  Real pc = pressure * 1.e-6;

  Real vc1 = 1.8882e-4 * tk;
  Real vc2 = - pc * (8.24e-2 + 1.249e-2 * pc) / tc;
  Real dvc2 = - (8.24e-2 + 2.498e-2 * pc) / tc;

  return (vc1 + vc2 - pc * dvc2) / ((vc1 + vc2) * (vc1 + vc2)) * 1e-6;
}

Real
dGasDensity_dT(Real pressure, Real temperature)
{
  Real tk = temperature + _t_c2k;
  Real tc = std::pow((tk * 1.e-2), 10./3.);
  Real pc = pressure * 1.e-6;

  Real vc1 = 1.8882e-4 * tk;
  Real vc2 = - pc * (8.24e-2 + 1.249e-2 * pc) / tc;
  Real dtc = (0.1 / 3.0) * std::pow((tk * 1.e-2), 7./3.);
  Real dvc1 = 1.8882e-4;
  Real dvc2 = - vc2 / tc;

  return - pc / (vc1 + vc2) / (vc1 + vc2) * (dvc1 + dvc2 * dtc);
}

Real
dSupercriticalDensity_dP(Real pressure, Real temperature)
{
  Real t1 = temperature;
  Real t2 = t1 * t1;
  Real t3 = t2 * t1;
  Real t4 = t3 * t1;

  // Correlation uses pressure in psia
  Real p1 = pressure * _pa2psia;
  Real p2 = p1 * p1;
  Real p3 = p2 * p1;

  Real a[5];

  if (p1 <= 3000)
    for (unsigned int i = 0; i < 5; ++i)
      a[i] = _b_scd[i][0] + _b_scd[i][1] * t1 + _b_scd[i][2] * t2 + _b_scd[i][3] * t3 + _b_scd[i][4] * t4;
  else
    for (unsigned int i = 0; i < 5; ++i)
      a[i] = _c_scd[i][0] + _c_scd[i][1] * t1 + _c_scd[i][2] * t2 + _c_scd[i][3] * t3 + _c_scd[i][4] * t4;

  return (a[1] + 2.0 * a[2] * p1 + 3.0 * a[3] * p2 + 4.0 * a[4] * p3) * _pa2psia;
}

Real
dSupercriticalDensity_dT(Real pressure, Real temperature)
{
  Real t1 = temperature;
  Real t2 = t1 * t1;
  Real t3 = t2 * t1;

  // Correlation uses pressure in psia
  Real p1 = pressure * _pa2psia;
  Real p2 = p1 * p1;
  Real p3 = p2 * p1;
  Real p4 = p3 * p1;

  Real a[5];

  if (p1 <= 3000)
    for (unsigned int i = 0; i < 5; ++i)
      a[i] = _b_scd[i][1] + 2.0 * _b_scd[i][2] * t1 + 3.0 * _b_scd[i][3] * t2 + 4.0 * _b_scd[i][4] * t3;
   else
     for (unsigned int i = 0; i < 5; ++i)
       a[i] = _c_scd[i][1] + 2.0 * _c_scd[i][2] * t1 + 3.0 * _c_scd[i][3] * t2 + 4.0 * _c_scd[i][4] * t3;

  return a[0] + a[1] * p1 + a[2] * p2 + a[3] * p3 + a[4] * p4;
}

Real
dViscosity_dDensity(Real pressure, Real temperature, Real density)
{
  Real dmu_drho;

  if (pressure <= _p_critical)
    dmu_drho = dGasViscosity_dDensity(temperature, density);
  else
    dmu_drho = dSupercriticalViscosity_dDensity(pressure, temperature);

  return dmu_drho;
}

Real
dGasViscosity_dDensity(Real temperature, Real density)
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

Real
dSupercriticalViscosity_dDensity(Real pressure, Real temperature)
{
  Real dmu_drho = 0.;
  /**
   * The correlation for supercritical CO2 gives viscosity as a function of pressure. Therefore, The
   * derivative of viscosity wrt density is given by the chain rule.
   * Note that if d(density)/d(pressure) = 0, so should d(viscosity)/d(density)
   */
  Real drho_dp = dSupercriticalDensity_dP(pressure, temperature);
  if (drho_dp != 0.0)
    dmu_drho += dSupercriticalViscosity_dP(pressure, temperature) / drho_dp;

  return dmu_drho;
}

Real
dSupercriticalViscosity_dP(Real pressure, Real temperature)
{
  Real a[5];

  Real t1 = temperature;
  Real t2 = t1 * t1;
  Real t3 = t2 * t1;
  Real t4 = t3 * t1;

  // Correlation uses pressure in psia
  Real p1 = pressure * _pa2psia;
  Real p2 = p1 * p1;
  Real p3 = p2 * p1;

  if (p1 <= 3000)
    for (unsigned int i = 0; i < 5; ++i)
      a[i] = _b_scv[i][0] + _b_scv[i][1] * t1 + _b_scv[i][2] * t2 + _b_scv[i][3] * t3 + _b_scv[i][4] * t4;
  else
    for (unsigned int i = 0; i < 5; ++i)
      a[i] = _c_scv[i][0] + _c_scv[i][1] * t1 + _c_scv[i][2] * t2 + _c_scv[i][3] * t3 + _c_scv[i][4] * t4;

 Real dmu = a[1] + 2.0 * a[2] * p1 + 3.0 * a[3] * p2 + 4.0 * a[4] * p3;

 return dmu * _pa2psia * 1e-3; // cP to Pa.s
}

std::vector<Real>
henryConstants()
{
  std::vector<Real> co2henry;
  co2henry.push_back(-8.55445);
  co2henry.push_back(4.01195);
  co2henry.push_back(9.52345);

  return co2henry;
}
}
