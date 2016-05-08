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
#include "PorousFlowMethanePropertiesTest.h"

CPPUNIT_TEST_SUITE_REGISTRATION( PorousFlowMethanePropertiesTest );

PorousFlowMethanePropertiesTest::PorousFlowMethanePropertiesTest() :
    _eps(1.0E-8)
{
}

void
PorousFlowMethanePropertiesTest::density()
{
  Real density = 11.466;
  Real delta = 0.02 * density;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(density, PorousFlowMethaneProperties::density(2.0e6, 70.0), delta);
  density = 27.761;
  delta = 0.02 * density;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(density, PorousFlowMethaneProperties::density(6.0e6, 150.0), delta);
  density = 75.844;
  delta = 0.1 * density;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(density, PorousFlowMethaneProperties::density(25.0e6, 320.0), delta);
}

void
PorousFlowMethanePropertiesTest::viscosity()
{
  Real viscosity = 12.835e-6;
  Real delta = 0.05 * viscosity;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(viscosity, PorousFlowMethaneProperties::viscosity(70), delta);
  viscosity = 15.003e-6;
  delta = 0.05 * viscosity;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(viscosity, PorousFlowMethaneProperties::viscosity(150), delta);
  /// Now test derivatives
  Real fd;
  fd = (PorousFlowMethaneProperties::viscosity(70.0 + _eps) - PorousFlowMethaneProperties::viscosity(70.0)) / _eps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowMethaneProperties::dViscosity_dT(70.0), 1.0E-10);
  fd = (PorousFlowMethaneProperties::viscosity(150.0 + _eps) - PorousFlowMethaneProperties::viscosity(150.0)) / _eps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowMethaneProperties::dViscosity_dT(150.0), 1.0E-10);
}
