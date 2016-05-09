/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowBrinePropertiesTest.h"

CPPUNIT_TEST_SUITE_REGISTRATION( PorousFlowBrinePropertiesTest );

PorousFlowBrinePropertiesTest::PorousFlowBrinePropertiesTest() :
    _eps(1.0E-8),
    _Mnacl(PorousFlowBrineProperties::molarMass()),
    _Mh2o(PorousFlowBrineProperties::molarMassH2O()),
    _t_c2k(273.15)
{
}

void
PorousFlowBrinePropertiesTest::density()
{
  /**
   * Note: the verification data uses mole fraction rather than mass fraction, so first
   * the mass fraction must be calculated
   */
  Real Xmass = mol2Mass(0.1081911);
  Real density = PorousFlowBrineProperties::density(50.0e6, 100.0, Xmass);
  Real delta = 0.001 * density;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.186431e+03, density, delta);
  Xmass = mol2Mass(1.402430e-01);
  density = PorousFlowBrineProperties::density(100.0e6, 250.0, Xmass);
  delta = 0.001 * density;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.154736e+03, density, delta);
  Xmass = mol2Mass(1.402536e-01);
  density = PorousFlowBrineProperties::density(6.351518e+06, 300.0, Xmass);
  delta = 0.001 * density;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.063782e+03, density, delta);
}

void
PorousFlowBrinePropertiesTest::pSat()
{
  /**
   * Note that brine vapour pressure is given in torr (multiply by 101325 / 760 to convert to Pa),
   * and molal concentration (mol/kg) is used. Experimental data for T = 200C.
   */
  Real Xmol = 7.973;
  Real Xmass = molality2Mass(Xmol);
  Real psat = PorousFlowBrineProperties::pSat(200.0, Xmass);
  Real delta = 0.01 * psat;
  Real exp_data = (11644.93 - 3272.6) * 101325.0 / 760.0;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(exp_data, psat, delta);
  Xmol = 3.886;
  Xmass = molality2Mass(Xmol);
  psat = PorousFlowBrineProperties::pSat(200.0, Xmass);
  delta = 0.01 * psat;
  exp_data = (11644.93 - 1607.7) * 101325.0 / 760.0;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(exp_data, psat, delta);
}

void
PorousFlowBrinePropertiesTest::viscosity()
{
  /**
   * Note: the verification data uses molal concentration rather than mass fraction, so first
   * the mass fraction must be calculated.
   */
  Real Xmol = 0.0;
  Real Xmass = molality2Mass(Xmol);
  Real viscosity = PorousFlowBrineProperties::viscosity(10.0e6, 50.0, Xmass);
  Real delta = 0.01 * viscosity;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5460e-3, viscosity, delta);
  Xmol = 5.0;
  Xmass = molality2Mass(Xmol);
  viscosity = PorousFlowBrineProperties::viscosity(10.0e6, 50.0, Xmass);
  delta = 0.01 * viscosity;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.9642e-3, viscosity, delta);
  Xmol = 3.0;
  Xmass = molality2Mass(Xmol);
  viscosity = PorousFlowBrineProperties::viscosity(20.0e6, 100.0, Xmass);
  delta = 0.01 * viscosity;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4056e-3, viscosity, delta);
}

void
PorousFlowBrinePropertiesTest::solubility()
{
  Real solubility = PorousFlowBrineProperties::haliteSolubility(386.5);
  Real delta = 0.02 * solubility;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.442, solubility, delta);
  solubility = PorousFlowBrineProperties::haliteSolubility(630.0);
  delta = 0.02 * solubility;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7185, solubility, delta);
}

Real
PorousFlowBrinePropertiesTest::mol2Mass(Real xmol)
{
  return xmol * _Mnacl /(xmol * _Mnacl + (1.0 - xmol) * _Mh2o);
}

Real
PorousFlowBrinePropertiesTest::molality2Mass(Real xmol)
{
  return 1.0 / (1.0 + 1.0 / (xmol * _Mnacl));
}
