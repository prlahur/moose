/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowSimpleCO2PropertiesTest.h"

CPPUNIT_TEST_SUITE_REGISTRATION( PorousFlowSimpleCO2PropertiesTest );

PorousFlowSimpleCO2PropertiesTest::PorousFlowSimpleCO2PropertiesTest() :
    _eps(1.0e-2),
    _teps(1.0e-4)
{
}

void
PorousFlowSimpleCO2PropertiesTest::density()
{
  /// Test the density
  Real density = 29.774;
  Real delta = 0.05 * density;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(density, PorousFlowSimpleCO2Properties::density(2.0e6, 100.0), delta);
  density = 99.632;
  delta = 0.05 * density;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(density, PorousFlowSimpleCO2Properties::density(6.0e6, 100.0), delta);
  density = 141.27;
  delta = 0.05 * density;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(density, PorousFlowSimpleCO2Properties::density(8.0e6, 100.0), delta);
  density = 480.53;
  delta = 0.05 * density;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(density, PorousFlowSimpleCO2Properties::density(20.0e6, 100.0), delta);
  /// Test derivatives wrt pressure
  Real fd;
  fd = (PorousFlowSimpleCO2Properties::density(2.0e6 + _eps, 100.0) - PorousFlowSimpleCO2Properties::density(2.0e6, 100.0)) / _eps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dDensity_dP(2.0e6, 100.0), 1.0e-8);
  fd = (PorousFlowSimpleCO2Properties::density(6.0e6 + _eps, 100.0) - PorousFlowSimpleCO2Properties::density(6.0e6, 100.0)) / _eps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dDensity_dP(6.0e6, 100.0), 1.0e-8);
  fd = (PorousFlowSimpleCO2Properties::density(8.0e6 + _eps, 100.0) - PorousFlowSimpleCO2Properties::density(8.0e6, 100.0)) / _eps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dDensity_dP(8.0e6, 100.0), 1.0e-8);
  fd = (PorousFlowSimpleCO2Properties::density(20.0e6 + _eps, 100.0) - PorousFlowSimpleCO2Properties::density(20.0e6, 100.0)) / _eps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dDensity_dP(20.0e6, 100.0), 1.0e-8);
  /// Test derivatives wrt temperature
  fd = (PorousFlowSimpleCO2Properties::density(2.0e6, 100.0 + _teps) - PorousFlowSimpleCO2Properties::density(2.0e6, 100.0)) / _teps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dDensity_dT(2.0e6, 100.0), 1.0e-3);
  fd = (PorousFlowSimpleCO2Properties::density(6.0e6, 100.0 + _teps) - PorousFlowSimpleCO2Properties::density(6.0e6, 100.0)) / _teps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dDensity_dT(6.0e6, 100.0), 1.0e-3);
  fd = (PorousFlowSimpleCO2Properties::density(8.0e6, 100.0 + _teps) - PorousFlowSimpleCO2Properties::density(8.0e6, 100.0)) / _teps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dDensity_dT(8.0e6, 100.0), 1.0e-3);
  fd = (PorousFlowSimpleCO2Properties::density(20.0e6, 100.0 + _teps) - PorousFlowSimpleCO2Properties::density(20.0e6, 100.0)) / _teps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dDensity_dT(20.0e6, 100.0), 1.0e-3);
}

void
PorousFlowSimpleCO2PropertiesTest::viscosity()
{
  /// Test viscosity
  Real density = PorousFlowSimpleCO2Properties::density(2.0e6, 100.0);
  Real viscosity = 18.654e-6;
  Real delta = 0.05 * viscosity;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(viscosity, PorousFlowSimpleCO2Properties::viscosity(2.0e6, 100.0, density), delta);
  density = PorousFlowSimpleCO2Properties::density(6.0e6, 100.0);
  viscosity = 19.589e-6;
  delta = 0.05 * viscosity;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(viscosity, PorousFlowSimpleCO2Properties::viscosity(6.0e6, 100.0, density), delta);
  density = PorousFlowSimpleCO2Properties::density(8.0e6, 100.0);
  viscosity = 20.481e-6;
  delta = 0.05 * viscosity;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(viscosity, PorousFlowSimpleCO2Properties::viscosity(8.0e6, 100.0, density), delta);
  density = PorousFlowSimpleCO2Properties::density(20.0e6, 100.0);
  viscosity = 37.190e-6;
  delta = 0.05 * viscosity;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(viscosity, PorousFlowSimpleCO2Properties::viscosity(20.0e6, 100.0, density), delta);
  /// Test derivatives wrt density (by varying pressure)
  Real fd;
  Real density2;
  density = PorousFlowSimpleCO2Properties::density(2.0e6, 100.0);
  density2 = PorousFlowSimpleCO2Properties::density(2.0e6 + _eps, 100.0);
  fd = (PorousFlowSimpleCO2Properties::viscosity(2.0e6, 100.0, density2) - PorousFlowSimpleCO2Properties::viscosity(2.0e6, 100.0, density)) / (density2 - density);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dViscosity_dDensity(2.0e6, 100.0, density), 1.0e-12);
  density = PorousFlowSimpleCO2Properties::density(6.0e6, 100.0);
  density2 = PorousFlowSimpleCO2Properties::density(6.0e6 + _eps, 100.0);
  fd = (PorousFlowSimpleCO2Properties::viscosity(6.0e6, 100.0, density2) - PorousFlowSimpleCO2Properties::viscosity(6.0e6, 100.0, density)) / (density2 - density);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dViscosity_dDensity(6.0e6, 100.0, density), 1.0e-12);
  density = PorousFlowSimpleCO2Properties::density(8.0e6, 100.0);
  density2 = PorousFlowSimpleCO2Properties::density(8.0e6 + _eps, 100.0);
  fd = (PorousFlowSimpleCO2Properties::viscosity(8.0e6 + _eps, 100.0, density2) - PorousFlowSimpleCO2Properties::viscosity(8.0e6, 100.0, density)) / (density2 - density);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dViscosity_dDensity(8.0e6, 100.0, density), 1.0e-12);
  density = PorousFlowSimpleCO2Properties::density(20.0e6, 100.0);
  density2 = PorousFlowSimpleCO2Properties::density(20.0e6 + _eps, 100.0);
  fd = (PorousFlowSimpleCO2Properties::viscosity(20.0e6 + _eps, 100.0, density2) - PorousFlowSimpleCO2Properties::viscosity(20.0e6, 100.0, density)) / (density2 - density);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowSimpleCO2Properties::dViscosity_dDensity(20.0e6, 100.0, density), 1.0e-12);
}

void
PorousFlowSimpleCO2PropertiesTest::partialDensity()
{
  /// Note that the experimental data must be scaled appropriately
  Real data = PorousFlowSimpleCO2Properties::_Mco2 * 1.0e6 / 50.0;
  Real delta = 0.02 * data;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(data, PorousFlowSimpleCO2Properties::partialDensity(200.0), delta);
  data = PorousFlowSimpleCO2Properties::_Mco2 * 1.0e6 / 74.1;
  delta = 0.02 * data;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(data, PorousFlowSimpleCO2Properties::partialDensity(300.0), delta);
}
