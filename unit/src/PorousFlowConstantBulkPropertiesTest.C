/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowConstantBulkPropertiesTest.h"

CPPUNIT_TEST_SUITE_REGISTRATION( PorousFlowConstantBulkPropertiesTest );

PorousFlowConstantBulkPropertiesTest::PorousFlowConstantBulkPropertiesTest()
{
}

void
PorousFlowConstantBulkPropertiesTest::density()
{
  /// Input variables
  Real pressure = 1.0e6;
  Real density_p0 = 300.0;
  Real bulk_modulus = 2.0e7;
  /// Computed density and derivatives
  Real density = density_p0 * std::exp(pressure / bulk_modulus);
  Real ddensity_dp = density / bulk_modulus;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(density, PorousFlowConstantBulkProperties::density(pressure, density_p0, bulk_modulus), 1.0E-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ddensity_dp, PorousFlowConstantBulkProperties::dDensity_dP(pressure, density_p0, bulk_modulus), 1.0E-12);
}
