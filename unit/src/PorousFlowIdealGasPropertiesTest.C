/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#include "PorousFlowIdealGasPropertiesTest.h"

CPPUNIT_TEST_SUITE_REGISTRATION( PorousFlowIdealGasPropertiesTest );

PorousFlowIdealGasPropertiesTest::PorousFlowIdealGasPropertiesTest() :
    _peps(1.0E-1),
    _teps(1.0E-8)
{
}

void
PorousFlowIdealGasPropertiesTest::density()
{
  /// Input variables
  Real pressure = 1.0e6;
  Real temperature = 150.0;
  Real molar_mass = 18.015e-3;
  Real t_c2k = PorousFlowIdealGasProperties::_t_c2k;
  Real R = PorousFlowIdealGasProperties::_R;
  /// Computed density and derivatives
  Real density = pressure * molar_mass / (R * (temperature + t_c2k));
  Real ddensity_dp = density / pressure;
  Real ddensity_dt = - density * density / pressure / molar_mass * R;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(density, PorousFlowIdealGasProperties::density(pressure, temperature, molar_mass), 1.0E-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ddensity_dp, PorousFlowIdealGasProperties::dDensity_dP(temperature, molar_mass), 1.0E-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ddensity_dt, PorousFlowIdealGasProperties::dDensity_dT(pressure, temperature, molar_mass), 1.0E-12);
}
