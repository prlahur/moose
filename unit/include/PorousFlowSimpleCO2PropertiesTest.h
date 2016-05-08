/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWSIMPLECO2PROPERTIESTEST_H
#define POROUSFLOWSIMPLECO2PROPERTIESTEST_H

//CPPUnit includes
#include "GuardedHelperMacros.h"

// Moose includes
#include "PorousFlowSimpleCO2Properties.h"

class PorousFlowSimpleCO2PropertiesTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( PorousFlowSimpleCO2PropertiesTest );

  /**
   * Verify calculation of CO2 density and derivatives wrt pressure and temperature.
   * Compare results with data from NIST webbook (www.nist.gov).
   * Results should agree to within 5%
   */
  CPPUNIT_TEST( density );

  /**
   * Verify calculation of CO2 viscosity and derivatives wrt density and
   * temperature. Compare results with data from NIST webbook (www.nist.gov).
   * Results should agree to within 5%
   */
  CPPUNIT_TEST( viscosity );

  /**
   * Verify calculation of CO2 partial density. Compare results with data from
   * Hnedkovsky et al, Volumes of aqueous solutions of CH4, CO2, H2S, and NH3
   * at temperatures from 298.15 K to 705 K and pressures to 35 MPa, J. Chem.
   * Thermodynamics 28, 125–142 (1996)
   */
  CPPUNIT_TEST( partialDensity );

  CPPUNIT_TEST_SUITE_END();

public:
  PorousFlowSimpleCO2PropertiesTest();

  void density();
  void viscosity();
  void partialDensity();

 private:
  Real _eps;
  Real _teps;
};

#endif  // POROUSFLOWSIMPLECO2PROPERTIESTEST_H
