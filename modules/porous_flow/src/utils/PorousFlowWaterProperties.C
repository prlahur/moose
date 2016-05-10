/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowWaterProperties.h"

namespace PorousFlowWaterProperties
{

std::string
fluidName()
{
  return "Water";
}

Real
molarMass()
{
  return _Mh2o;
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

unsigned int
inRegion(Real pressure, Real temperature)
{
  /**
   * Valid for 273.15 K <= T <= 1073.15 K, p <= 100 MPa
   *          1073.15 K <= T <= 2273.15 K, p <= 50 Mpa
   *
   * Check that the provided pressure and temperature is in the valid region
   */
  if (temperature > 0.0 && temperature <= 800.0)
  {
    if (pressure < 0.0 || pressure > 100.0e6)
      mooseError("Pressure "<< pressure << " is out of range in PorousFlowWaterProperties::inRegion");
  }
  else if (temperature > 800.0 && temperature <= 2000.0)
  {
    if (pressure < 0.0 || pressure > 50.0e6)
      mooseError("Pressure "<< pressure << " is out of range in PorousFlowWaterProperties::inRegion");
  }
  else
    mooseError("Temperature " << temperature << " is out of range in PorousFlowWaterProperties::inRegion");

  /// Determine the phase region that the (P, T) point lies in
  unsigned int region;

  if (temperature >= 0.0 && temperature <= 350.0)
  {
    if (pressure > pSat(temperature) && pressure <= 100.0e6)
      region = 1;
    else
      region = 2;
  }
  else if (temperature > 350.0 && temperature <= 590.0)
  {
    if (pressure <= b23p(temperature))
      region = 2;
    else
      region = 3;
  }
  else if (temperature > 590.0 && temperature <= 800.0)
    region = 2;
  else
   region = 5;

  return region;
}

Real
density(Real pressure, Real temperature)
{
  Real density;

  /// Determine which region the point is in
  unsigned int region = inRegion(pressure, temperature);

  switch(region)
  {
    case 1:
      density = densityRegion1(pressure, temperature);
      break;

    case 2:
      density = densityRegion2(pressure, temperature);
      break;

    case 3:
      density = densityRegion3(pressure, temperature);
      break;

    case 5:
      density = densityRegion5(pressure, temperature);
      break;
  }

  return density;
}

Real
viscosity(Real temperature, Real density)
{
  Real t0[4], t1[6], d1[7];

  Real mu_star = 1.e-6;

  Real tbar = (temperature + _t_c2k) / _t_critical;
  Real rhobar = density / _rho_critical;

  t0[0] = 1.0;
  t0[1] = 1.0 / tbar;
  t0[2] = t0[1] * t0[1];
  t0[3] = t0[2] * t0[1];

  t1[0] = 1.0;
  t1[1] = 1.0 / tbar - 1.0;
  t1[2] = t1[1] * t1[1];
  t1[3] = t1[2] * t1[1];
  t1[4] = t1[3] * t1[1];
  t1[5] = t1[4] * t1[1];

  d1[0] = 1.0;
  d1[1] = rhobar - 1.0;
  d1[2] = d1[1] * d1[1];
  d1[3] = d1[2] * d1[1];
  d1[4] = d1[3] * d1[1];
  d1[5] = d1[4] * d1[1];
  d1[6] = d1[5] * d1[1];

  /// Calculate mu0
  Real sum0 = 0.0;
  for (unsigned int i = 0; i < 4; ++i)
     sum0 += _h0v[i] * t0[i];

  Real mu0 = 100.0 * std::sqrt(tbar) / sum0;

  /// Now calculate mu1
  Real sum1 = 0.0;
  for (unsigned int i = 0; i < 21; ++i)
     sum1 += t1[_iv[i]] * _h1v[i] * d1[_jv[i]];

  Real mu1 = std::exp(rhobar * sum1);

  /// The water viscosity (in Pa.s) is then given by
  return mu_star * mu0 * mu1;
}

Real
pSat(Real temperature)
{
  Real tk = temperature + _t_c2k;
  Real theta, theta2, a, b, c, p;

  /**
   * Check whether the input temperature is within the region of validity of this equation.
   * Valid for 273.15 K <= t <= 647.096 K
   */
  if (tk >= _t_c2k && tk <= _t_critical)
  {
    theta = tk + _n4[8] / (tk - _n4[9]);
    theta2 = theta * theta;
    a = theta2 + _n4[0] * theta + _n4[1];
    b = _n4[2] * theta2 + _n4[3] * theta + _n4[4];
    c = _n4[5] * theta2 + _n4[6] * theta + _n4[7];
    p = std::pow(2.0 * c/(-b + std::sqrt(b * b - 4.0 * a * c)), 4.0);
  }
  else
    mooseError("PorousFlowWaterProperties::pSat: Temperature is outside range 0 C <= t <= 400.964 C");

  return p * 1.e6;
}

Real
tSat(Real pressure)
{
  Real beta, beta2, e, f, g, d, t;

  /**
   * Check whether the input pressure is within the region of validity of this equation.
   * Valid for 611.213 Pa <= p <= 22.064 MPa
   */
  if (pressure >= 611.23 && pressure <= _p_critical)
  {
    beta = std::pow(pressure / 1.e6, 0.25);
    beta2 = beta * beta;
    e = beta2 + _n4[2] * beta + _n4[5];
    f = _n4[0] * beta2 + _n4[3] * beta + _n4[6];
    g = _n4[1] * beta2 + _n4[4] * beta + _n4[7];
    d = 2.0 * g /(-f - std::sqrt(f * f - 4.0 * e * g));
    t = (_n4[9] + d - std::sqrt((_n4[9] + d) * (_n4[9] + d) - 4.0 * (_n4[8] + _n4[9] * d))) / 2.0;
  }
  else
    mooseError("PorousFlowWaterProperties::tSat: Pressure is outside range 611.213 Pa <= p <= 22.064 MPa");

  return t - _t_c2k;
}

Real
b23p(Real temperature)
{
  Real tk = temperature + _t_c2k;
  Real p23;

  /**
   * Check whether the input temperature is within the region of validity of this equation.
   * Valid for 623.15 K <= t <= 863.15 K
   */
  if (tk >= 623.15 && tk <= 863.15)
    p23 = (_n23[0] + _n23[1] * tk + _n23[2] * tk * tk) * 1.e6;
  else
    mooseError("PorousFlowWaterProperties::b23p: Temperature is outside range of 350 C <= t <= 590 C");

  return p23;
}

Real
b23t(Real pressure)
{
  Real t23;

  /**
   * Check whether the input pressure is within the region of validity of this equation.
   * Valid for 16.529 MPa <= p <= 100 MPa
   */
  if (pressure >= 16.529e6 && pressure <= 100.0e6)
    t23 = _n23[3] + std::sqrt((pressure / 1.e6 - _n23[4]) / _n23[2]) - _t_c2k;
  else
    mooseError("PorousFlowWaterProperties::b23t: Pressure is outside range 16.529 MPa <= p <= 100 MPa");

  return t23;
}

Real
densityRegion1(Real pressure, Real temperature)
{
  Real p_star1 = 16.53e6;
  Real t_star1 = 1386.0;
  Real tk = temperature + _t_c2k;

  /// Now evaluate the sums
  Real sum1 = 0.0;
  Real tau1 = t_star1 / tk;
  Real pi1 = pressure / p_star1;

  for (unsigned int i = 0; i < 34; ++i)
  {
    sum1 -= _n1[i] * _I1[i] * std::pow(7.1 - pi1, _I1[i] - 1) * std::pow(tau1 - 1.222, _J1[i]);
}

  /// The density of water in this region is then given by
  return p_star1 / (sum1 * _Rw * tk) / 1000.0;
}

Real
densityRegion2(Real pressure, Real temperature)
{
  Real p_star2 = 1.e6;
  Real t_star2 = 540.0;
  Real tk = temperature + _t_c2k;

  /// Ideal gas component of region 2 - Eq. (16)
  Real tau2 = t_star2 / tk;
  Real pi2 = pressure / p_star2;

  /// Residual component of Gibbs free energy - Eq. (17).
  Real sumr2 = 0.0;

  for (unsigned int i = 0; i < 43; ++i)
    sumr2 += _n2[i] * _I2[i] * std::pow(pi2, _I2[i] - 1) * std::pow(tau2 - 0.5, _J2[i]);

  /// The density in Region 2 is then given by
  return p_star2 / (_Rw * tk * (1.0 / pi2 + sumr2)) / 1000.0;
}

Real
dDensity_dP(Real pressure, Real temperature)
{
  Real ddensity = 0.0;

  /// Determine which region the point is in
  unsigned int region = inRegion(pressure, temperature);

  switch(region)
  {
    case 1:
      ddensity = dDensityRegion1_dP(pressure, temperature);
      break;

    case 2:
      ddensity = dDensityRegion2_dP(pressure, temperature);
      break;

    case 3:
      ddensity = dDensityRegion3_dP(pressure, temperature);
      break;

    case 5:
      ddensity = dDensityRegion5_dP(pressure, temperature);
      break;
  }

  return ddensity;
}

Real
dDensityRegion1_dP(Real pressure, Real temperature)
{
  Real p_star1 = 16.53e6;
  Real t_star1 = 1386.0;
  Real tk = temperature + _t_c2k;

  /// Now evaluate the sums
  Real sum1 = 0.0;
  Real sum2 = 0.0;
  Real tau1 = t_star1 / tk;
  Real pi1 = pressure / p_star1;

  for (unsigned int i = 0; i < 34; ++i)
  {
    sum1 -= _n1[i] * _I1[i] * std::pow(7.1 - pi1, _I1[i] - 1) * std::pow(tau1 - 1.222, _J1[i]);
    sum2 += _n1[i] * _I1[i] * (_I1[i] - 1) * std::pow(7.1 - pi1, _I1[i] - 2) * std::pow(tau1 - 1.222, _J1[i]);
  }

  /// The derivative of the density of water with respect to pressure in this region is then given by
  return - sum2 / (sum1 * sum1 * _Rw * tk) / 1000.0;
}

Real
dDensityRegion2_dP(Real pressure, Real temperature)
{
  Real p_star2 = 1.e6;
  Real t_star2 = 540.0;
  Real tk = temperature + _t_c2k;

  /// Ideal gas component of region 2 - Eq. (16)
  Real tau2 = t_star2 / tk;
  Real pi2 = pressure / p_star2;

  /// Residual component of Gibbs free energy - Eq. (17).
  Real sumr2 = 0.0;
  Real sumdr2 = 0.0;

  for (unsigned int i = 0; i < 43; ++i)
  {
    sumr2 += _n2[i] * _I2[i] * std::pow(pi2, _I2[i] - 1) * std::pow(tau2 - 0.5, _J2[i]);
    sumdr2 += _n2[i] * _I2[i] * (_I2[i] - 1) * std::pow(pi2, _I2[i] - 2) * std::pow(tau2 - 0.5, _J2[i]);
  }

  /// The derivative of the density in Region 2 with respect to pressure is then given by
  return - (- 1.0 / (pi2 * pi2) + sumdr2) / (_Rw * tk * (1.0 / pi2 + sumr2) * (1.0 / pi2 + sumr2)) / 1000.0;
}

Real
dViscosity_dDensity(Real temperature, Real density)
{
  Real t1[6], d1[7];

  Real tbar = (temperature + _t_c2k) / _t_critical;
  Real rhobar = density / _rho_critical;

  t1[0] = 1.0;
  t1[1] = 1.0 / tbar - 1.0;
  t1[2] = t1[1] * t1[1];
  t1[3] = t1[2] * t1[1];
  t1[4] = t1[3] * t1[1];
  t1[5] = t1[4] * t1[1];

  d1[0] = 1.0;
  d1[1] = rhobar - 1.0;
  d1[2] = d1[1] * d1[1];
  d1[3] = d1[2] * d1[1];
  d1[4] = d1[3] * d1[1];
  d1[5] = d1[4] * d1[1];
  d1[6] = d1[5] * d1[1];

  /// Prefactor to derivative of viscosity
  Real sum1 = 0.0;
  Real sum2 = 0.0;
  for (unsigned int i = 0; i < 21; ++i)
  {
    sum1 += t1[_iv[i]] * _h1v[i] * d1[_jv[i]];
    sum2 += t1[_iv[i]] * _jv[i] * _h1v[i] * d1[_jv[i]] / (rhobar - 1.0);
  }

  /// The derivative of viscosity wrt density is then
  return viscosity(temperature, density) * (sum1 + rhobar * sum2) / _rho_critical;
}

Real
dDensity_dT(Real /* pressure */, Real /* temperature */)
{
  return 0.; // FIXME: not implemented yet
}

Real
tempXY(Real pressure, std::string xy)
{
  Real pi = pressure / 1.0e6;

  const int I[11][5] = {{0, 1, 2, -1, -2}, {0, 1, 2, 3, 4}, {0, 1, 2, 3, 4}, {0, 1, 2, 3, 4}, {0, 1, 2, 3, 4},
                       {0, 1, 2, 3, 4}, {0, 1, 2, -1, -2}, {0, 1, 2, 3, 0}, {0, 1, 2, 3, 4}, {0, 1, 2, 3, 4}, {0, 1, 2, -1, -2}};

  const Real n[11][5] = {{0.154793642129415e4, -0.187661219490113e3, 0.213144632222113e2, -0.191887498864292e4, 0.918419702359447e3},
                        {0.585276966696349e3, 0.278233532206915e1, -0.127283549295878e-1, 0.159090746562729e-3, 0.0},
                        {-0.249284240900418e5, 0.428143584791546e4, -0.269029173140130e3, 0.751608051114157e1, -0.787105249910383e-1},
                        {0.584814781649163e3, -0.616179320924617, 0.260763050899562, -0.587071076864459e-2, 0.515308185433082e-4},
                        {0.617229772068439e3, -0.770600270141675e1, 0.697072596851896, -0.157391839848015e-1, 0.137897492684194e-3},
                        {0.535339483742384e3, 0.761978122720128e1, -0.158365725441648, 0.192871054508108e-2, 0.0},
                        {0.969461372400213e3, -0.332500170441278e3, 0.642859598466067e2, 0.773845935768222e3, -0.152313732937084e4},
                        {0.565603648239126e3, 0.529062258221222e1, -0.102020639611016, 0.122240301070145e-2, 0.0},
                        {0.584561202520006e3, -0.102961025163669e1, 0.243293362700452, -0.294905044740799e-2, 0.0},
                        {0.528199646263062e3, 0.890579602135307e1, -0.222814134903755, 0.286791682263697e-2},
                        {0.728052609145380e1, 0.973505869861952e2, 0.147370491183191e2, 0.329196213998375e3, 0.873371668682417e3}};

  /// Choose the constants based on the string xy
  unsigned int row;

  if (xy == "ab")
    row = 0;
  else if (xy == "cd")
    row = 1;
  else if (xy == "gh")
    row = 2;
  else if (xy == "ij")
    row = 3;
  else if (xy == "jk")
    row = 4;
  else if (xy == "mn")
    row = 5;
  else if (xy == "op")
    row = 6;
  else if (xy == "qu")
    row = 7;
  else if (xy == "rx")
    row = 8;
  else if (xy == "uv")
    row = 9;
  else if (xy == "wx")
    row = 10;
  else if (xy == "ef")
    row = 0; // not used
  else
    mooseError("PorousFlowWaterProperties::tempXY: Invalid boundary specified");

  Real sum = 0.0;

  if (xy == "ab" || xy == "op" || xy == "wx")
    for (unsigned int i = 0; i < 5; ++i)
      sum += n[row][i] * std::pow(std::log(pi), I[row][i]);
  else if (xy == "ef")
      sum += 3.727888004 * (pi - _p_critical / 1.0e6) + _t_critical;
  else
    for (unsigned int i = 0; i < 5; ++i)
      sum += n[row][i] * std::pow(pi, I[row][i]);

  return sum - _t_c2k;
}

unsigned int
subregion3(Real pressure, Real temperature)
{
  Real pMPa = pressure / 1.0e6;
  const Real P3cd  = 19.00881189173929;
  unsigned int subregion = 0;

  if (pMPa > 40.0 && pMPa <= 100.0)
  {
    if (temperature <= tempXY(pressure, "ab"))
      subregion = 0;
    else if (temperature > tempXY(pressure, "ab"))
      subregion = 1;
    else
      mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 1");
  }
  else if (pMPa > 25.0 && pMPa <= 40.0)
  {
    if (temperature <= tempXY(pressure, "cd"))
      subregion = 2;
    else if (temperature > tempXY(pressure, "cd") && temperature <= tempXY(pressure, "ab"))
      subregion = 3;
    else if (temperature > tempXY(pressure, "ab") && temperature <= tempXY(pressure, "ef"))
      subregion = 4;
    else if (temperature > tempXY(pressure, "ef"))
      subregion = 5;
    else
      mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 2");
  }
  else if (pMPa > 23.5 && pMPa <= 25.0)
  {
    if (temperature <= tempXY(pressure, "cd"))
      subregion = 2;
    else if (temperature > tempXY(pressure, "cd") && temperature <= tempXY(pressure, "gh"))
      subregion = 6;
    else if (temperature > tempXY(pressure, "gh") && temperature <= tempXY(pressure, "ef"))
      subregion = 7;
    else if (temperature > tempXY(pressure, "ef") && temperature <= tempXY(pressure, "ij"))
      subregion = 8;
    else if (temperature > tempXY(pressure, "ij") && temperature <= tempXY(pressure, "jk"))
      subregion = 9;
    else if (temperature > tempXY(pressure, "jk"))
      subregion = 10;
    else
      mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 3");
  }
  else if (pMPa > 23.0 && pMPa <= 23.5)
  {
    if (temperature <= tempXY(pressure, "cd"))
      subregion = 2;
    else if (temperature > tempXY(pressure, "cd") && temperature <= tempXY(pressure, "gh"))
      subregion = 11;
    else if (temperature > tempXY(pressure, "gh") && temperature <= tempXY(pressure, "ef"))
      subregion = 7;
    else if (temperature > tempXY(pressure, "ef") && temperature <= tempXY(pressure, "ij"))
      subregion = 8;
    else if (temperature > tempXY(pressure, "ij") && temperature <= tempXY(pressure, "jk"))
      subregion = 9;
    else if (temperature > tempXY(pressure, "jk"))
      subregion = 10;
    else
      mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 4");
  }
  else if (pMPa > 22.5 && pMPa <= 23.0)
  {
    if (temperature <= tempXY(pressure, "cd"))
      subregion = 2;
    else if (temperature > tempXY(pressure, "cd") && temperature <= tempXY(pressure, "gh"))
      subregion = 11;
    else if (temperature > tempXY(pressure, "gh") && temperature <= tempXY(pressure, "mn"))
      subregion = 12;
    else if (temperature > tempXY(pressure, "mn") && temperature <= tempXY(pressure, "ef"))
      subregion = 13;
    else if (temperature > tempXY(pressure, "ef") && temperature <= tempXY(pressure, "op"))
      subregion = 14;
    else if (temperature > tempXY(pressure, "op") && temperature <= tempXY(pressure, "ij"))
      subregion = 15;
    else if (temperature > tempXY(pressure, "ij") && temperature <= tempXY(pressure, "jk"))
      subregion = 9;
    else if (temperature > tempXY(pressure, "jk"))
      subregion = 10;
    else
      mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 5");
  }
  else if (pMPa > pSat(370.0) * 1.0e-6 && pMPa <= 22.5) /// pSat(370.0) = 21.04 MPa
  {
    if (temperature <= tempXY(pressure, "cd"))
      subregion = 2;
    else if (temperature > tempXY(pressure, "cd") && temperature <= tempXY(pressure, "qu"))
      subregion = 16;
    else if (temperature > tempXY(pressure, "qu") && temperature <= tempXY(pressure, "rx"))
    {
      if (pMPa > 22.11 && pMPa <= 22.5)
      {
        if (temperature <= tempXY(pressure, "uv"))
          subregion = 20;
        else if (temperature > tempXY(pressure, "uv") && temperature <= tempXY(pressure, "ef"))
          subregion = 21;
        else if (temperature > tempXY(pressure, "ef") && temperature <= tempXY(pressure, "wx"))
          subregion = 22;
        else if (temperature > tempXY(pressure, "wx") && temperature <= tempXY(pressure, "rx"))
          subregion = 23;
        else
          mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 6");
      }
      else if (pMPa > 22.064 && pMPa <= 22.11)
      {
        if (temperature <= tempXY(pressure, "uv"))
          subregion = 20;
        else if (temperature > tempXY(pressure, "uv") && temperature <= tempXY(pressure, "ef"))
          subregion = 24;
        else if (temperature > tempXY(pressure, "ef") && temperature <= tempXY(pressure, "wx"))
          subregion = 25;
        else  if (temperature > tempXY(pressure, "wx") && temperature <= tempXY(pressure, "rx"))
          subregion = 23;
        else
          mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 7");
      }
      else if (temperature <= tSat(pressure))
      {
        if (pMPa > 21.93161551 && pMPa <= 22.064)
          if (temperature > tempXY(pressure, "qu") && temperature <= tempXY(pressure, "uv"))
            subregion = 20;
          else
            subregion = 24;
        else if (pMPa > pSat(643.15 - _t_c2k) * 1.0e-6 && pMPa <= 21.93161551)
          subregion = 20;
        else
          mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 8");
      }
      else if (temperature > tSat(pressure))
      {
        if (pMPa > 21.90096265 && pMPa <= 22.064)
          if (temperature <= tempXY(pressure, "wx"))
            subregion = 25;
          else
            subregion = 23;
        else
          subregion = 23;
      }
    }
    else if (temperature > tempXY(pressure, "rx") && temperature <= tempXY(pressure, "jk"))
      subregion = 17;
    else
      subregion = 10;
  }
  else if (pMPa > 20.5 && pMPa <= pSat(370.0) * 1.0e-6) /// pSat(370.0) = 21.04 MPa
  {
    if (temperature <= tempXY(pressure, "cd"))
      subregion = 2;
    else if (temperature > tempXY(pressure, "cd") && temperature <= tSat(pressure))
      subregion = 16;
    else if (temperature > tSat(pressure) && temperature <= tempXY(pressure, "jk"))
      subregion = 17;
    else if (temperature > tempXY(pressure, "jk"))
      subregion = 10;
    else
      mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 9");
  }
  else if (pMPa > P3cd && pMPa <= 20.5) /// P3cd  = 19.00881189173929
  {
    if (temperature <= tempXY(pressure, "cd"))
      subregion = 2;
    else if (temperature > tempXY(pressure, "cd") && temperature <= tSat(pressure))
      subregion = 18;
    else
      subregion = 19;
  }
  else if (pMPa > pSat(350.0) * 1.0e-6 && pMPa <= P3cd)
  {
    if (temperature < tSat(pressure))
      subregion = 2;
    else
      subregion = 19;
  }
  else if (pMPa > 22.11 && pMPa <= 22.5)
  {
    if (temperature > tempXY(pressure, "qu") && temperature <= tempXY(pressure, "uv"))
      subregion = 20;
    else if (temperature > tempXY(pressure, "uv") && temperature <= tempXY(pressure, "ef"))
      subregion = 21;
    else if (temperature > tempXY(pressure, "ef") && temperature <= tempXY(pressure, "wx"))
      subregion = 22;
    else if (temperature > tempXY(pressure, "wx") && temperature <= tempXY(pressure, "rx"))
      subregion = 23;
    else
      mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 10");
  }
  else if (pMPa > 22.064 && pMPa <= 22.11)
  {
    if (temperature > tempXY(pressure, "qu") && temperature <= tempXY(pressure, "uv"))
      subregion = 20;
    else if (temperature > tempXY(pressure, "uv") && temperature <= tempXY(pressure, "ef"))
      subregion = 24;
    else if (temperature > tempXY(pressure, "ef") && temperature <= tempXY(pressure, "wx"))
      subregion = 25;
    else if (temperature > tempXY(pressure, "wx") && temperature <= tempXY(pressure, "rx"))
      subregion = 23;
    else
      mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 11");
  }
  else
    mooseError("PorousFlowWaterProperties::subregion3. Shouldn't have got here. Error location 12");

  return subregion;
}

Real
densityRegion3(Real pressure, Real temperature)
{
  /**
   * Region 3 is subdivided into 26 subregions, each with a given backwards equation
   * to directly calculate density from pressure and temperature without the need for
   * expensive iterations.
   */

  /// The subregion that the point is in:
  unsigned int subregion = subregion3(pressure, temperature);

  Real vstar, pi, theta, a, b, c, d, e;
  unsigned int N;

  vstar = _par3[subregion][0];
  pi = pressure / _par3[subregion][1] / 1.0e6;
  theta = (temperature + _t_c2k) / _par3[subregion][2];
  a = _par3[subregion][3];
  b = _par3[subregion][4];
  c = _par3[subregion][5];
  d = _par3[subregion][6];
  e = _par3[subregion][7];
  N = _par3N[subregion];

  Real sum = 0.0;
  Real volume = 0.0;

  switch(subregion)
  {
    case 0:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3a, _J3a, _n3a);
      break;

    case 1:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3b, _J3b, _n3b);
      break;

    case 2:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3c, _J3c, _n3c);
      break;

    case 3:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3d, _J3d, _n3d);
     break;

    case 4:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3e, _J3e, _n3e);
      break;

    case 5:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3f, _J3f, _n3f);
      break;

    case 6:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3g, _J3g, _n3g);
      break;

    case 7:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3h, _J3h, _n3h);
      break;

    case 8:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3i, _J3i, _n3i);
      break;

    case 9:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3j, _J3j, _n3j);
      break;

    case 10:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3k, _J3k, _n3k);
      break;

    case 11:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3l, _J3l, _n3l);
      break;

    case 12:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3m, _J3m, _n3m);
      break;

    case 13: /// Note: this is the only different case
      for (unsigned int i = 0; i < N; ++i)
        sum += _n3n[i] * std::pow(pi - a, _I3n[i]) * std::pow(theta - b, _J3n[i]);

      volume = vstar * std::exp(sum);
      break;

    case 14:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3o, _J3o, _n3o);
      break;

    case 15:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3p, _J3p, _n3p);
      break;

    case 16:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3q, _J3q, _n3q);
      break;

    case 17:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3r, _J3r, _n3r);
      break;

    case 18:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3s, _J3s, _n3s);
      break;

    case 19:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3t, _J3t, _n3t);
      break;

    case 20:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3u, _J3u, _n3u);
      break;

    case 21:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3v, _J3v, _n3v);
      break;

    case 22:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3w, _J3w, _n3w);
      break;

    case 23:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3x, _J3x, _n3x);
      break;

    case 24:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3y, _J3y, _n3y);
      break;

    case 25:
      volume = vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3z, _J3z, _n3z);
      break;
  }

  /// Density is the inverse of volume
  return 1.0 / volume;
}

Real
dDensityRegion3_dP(Real pressure, Real temperature)
{
  /**
   * Region 3 is subdivided into 26 subregions, each with a given backwards equation
   * to directly calculate density from pressure and temperature without the need for
   * expensive iterations.
   */

   /// TODO: calculates density again - pass in density as function to avoid recompilation

  /// The subregion that the point is in:
  unsigned int subregion = subregion3(pressure, temperature);

  Real vstar, pi, theta, a, b, c, d, e;
  unsigned int N;

  vstar = _par3[subregion][0];
  pi = pressure / _par3[subregion][1] / 1.0e6;
  theta = (temperature + _t_c2k) / _par3[subregion][2];
  a = _par3[subregion][3];
  b = _par3[subregion][4];
  c = _par3[subregion][5];
  d = _par3[subregion][6];
  e = _par3[subregion][7];
  N = _par3N[subregion];

  Real sum = 0.0;
  Real density = 0.0;
  Real dvolume = 0.0;
  Real dsum = 0.0;

  switch(subregion)
  {
    case 0:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3a, _J3a, _n3a));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3a, _J3a, _n3a);
      break;

    case 1:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3b, _J3b, _n3b));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3b, _J3b, _n3b);
      break;

    case 2:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3c, _J3c, _n3c));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3c, _J3c, _n3c);
      break;

    case 3:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3d, _J3d, _n3d));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3d, _J3d, _n3d);
     break;

    case 4:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3e, _J3e, _n3e));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3e, _J3e, _n3e);
      break;

    case 5:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3f, _J3f, _n3f));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3f, _J3f, _n3f);
      break;

    case 6:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3g, _J3g, _n3g));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3g, _J3g, _n3g);
      break;

    case 7:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3h, _J3h, _n3h));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3h, _J3h, _n3h);
      break;

    case 8:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3i, _J3i, _n3i));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3i, _J3i, _n3i);
      break;

    case 9:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3j, _J3j, _n3j));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3j, _J3j, _n3j);
      break;

    case 10:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3k, _J3k, _n3k));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3k, _J3k, _n3k);
      break;

    case 11:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3l, _J3l, _n3l));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3l, _J3l, _n3l);
      break;

    case 12:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3m, _J3m, _n3m));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3m, _J3m, _n3m);
      break;

    case 13: /// Note: this is the only different case
      for (unsigned int i = 0; i < N; ++i)
      {
        sum += _n3n[i] * std::pow(pi - a, _I3n[i]) * std::pow(theta - b, _J3n[i]);
        dsum += _n3n[i] * _I3n[i] * std::pow(pi - a, _I3n[i] - 1) * std::pow(theta - b, _J3n[i]);
      }
      density = 1.0 / (vstar * std::exp(sum));
      dvolume = vstar * dsum * std::exp(sum);
      break;

    case 14:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3o, _J3o, _n3o));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3o, _J3o, _n3o);
      break;

    case 15:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3p, _J3p, _n3p));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3p, _J3p, _n3p);
      break;

    case 16:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3q, _J3q, _n3q));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3q, _J3q, _n3q);
      break;

    case 17:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3r, _J3r, _n3r));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3r, _J3r, _n3r);
      break;

    case 18:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3s, _J3s, _n3s));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3s, _J3s, _n3s);
      break;

    case 19:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3t, _J3t, _n3t));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3t, _J3t, _n3t);
      break;

    case 20:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3u, _J3u, _n3u));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3u, _J3u, _n3u);
      break;

    case 21:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3v, _J3v, _n3v));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3v, _J3v, _n3v);
      break;

    case 22:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3w, _J3w, _n3w));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3w, _J3w, _n3w);
      break;

    case 23:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3x, _J3x, _n3x));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3x, _J3x, _n3x);
      break;

    case 24:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3y, _J3y, _n3y));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3y, _J3y, _n3y);
      break;

    case 25:
      density = 1.0 / (vstar * subregionVolume(pi, theta, a, b, c, d, e, N, _I3z, _J3z, _n3z));
      dvolume = vstar * dSubregionVolume_dP(pi, theta, a, b, c, d, e, N, _I3z, _J3z, _n3z);
      break;
  }

  return - density * density * dvolume / _par3[subregion][1] / 1.0e6 ;
}

Real
subregionVolume(Real pi, Real theta, Real a, Real b, Real c, Real d, Real e, unsigned int N, const int I[], const int J[], const Real n[])
{
  Real sum = 0.0;

  for (unsigned int i = 0; i < N; ++i)
    sum += n[i] * std::pow(std::pow(pi - a, c), I[i]) * std::pow(std::pow(theta - b, d), J[i]);

  return std::pow(sum, e);
}

Real
dSubregionVolume_dP(Real pi, Real theta, Real a, Real b, Real c, Real d, Real e, unsigned int N, const int I[], const int J[], const Real n[])
{
  Real sum = 0.0;
  Real dsum = 0.0;

  for (unsigned int i = 0; i < N; ++i)
  {
    sum += n[i] * std::pow(std::pow(pi - a, c), I[i]) * std::pow(std::pow(theta - b, d), J[i]);
    dsum += n[i] * I[i] * c * std::pow(pi - a, c - 1.0) * std::pow(std::pow(pi - a, c), I[i] - 1) * std::pow(std::pow(theta - b, d), J[i]);
  }

  return e * std::pow(sum, e - 1.0) * dsum;
}

Real
densityRegion5(Real pressure, Real temperature)
{
  Real tk = temperature + _t_c2k;
  Real t_star5 = 1000.0;
  Real p_star5 = 1.0e6;

  Real pi = pressure / p_star5;
  Real tau = t_star5 / tk;

  /// The ideal gas part
  Real part0 = 1.0 / pi;

  /// The residual part
  Real part1 = 0.0;
  for (unsigned int i = 0; i < 6; ++i)
    part1 += _n5r[i] * _I5r[i] * std::pow(pi, _I5r[i] - 1) * std::pow(tau, _J5r[i]);

  /// The density is then
  return p_star5 / (_Rw * tk * (part0 + part1)) / 1000.0;
}

Real
dDensityRegion5_dP(Real pressure, Real temperature)
{
  Real tk = temperature + _t_c2k;
  Real t_star5 = 1000.0;
  Real p_star5 = 1.0e6;

  Real pi = pressure / p_star5;
  Real tau = t_star5 / tk;

  /// The derivative of the ideal gas part wrt pressure
  Real part0 = 1.0 / pi;
  Real dpart0 = - 1.0 / pi / pi;

  /// The derivative of the residual part wrt pressure
  Real part1 = 0.0;
  Real dpart1 = 0.0;
  for (unsigned int i = 0; i < 6; ++i)
  {
    part1 += _n5r[i] * _I5r[i] * std::pow(pi, _I5r[i] - 1) * std::pow(tau, _J5r[i]);
    dpart1 += _n5r[i] * _I5r[i] * (_I5r[i] - 1) * std::pow(pi, _I5r[i] - 2) * std::pow(tau, _J5r[i]);
  }

  /// The derivative of density wrt pressure is then
  return - (dpart0 + dpart1) / (_Rw * tk * (part0 + part1) * (part0 + part1)) / 1000.0;
}
}
