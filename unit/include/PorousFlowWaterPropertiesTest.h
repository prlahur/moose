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

#ifndef POROUSFLOWWATERPROPERTIESTEST_H
#define POROUSFLOWWATERPROPERTIESTEST_H

//CPPUnit includes
#include "GuardedHelperMacros.h"

// Moose includes
#include "PorousFlowWaterProperties.h"

class PorousFlowWaterPropertiesTest : public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE( PorousFlowWaterPropertiesTest );

  /**
   * Verify that the correct region is provided for a given pressure and
   * temperature. Also verify that an error is thrown if pressure and temperature
   * are outside the range of validity
   */
  CPPUNIT_TEST( inRegion );

  /**
   * Verify calculation of water properties in region 1 using the verification values
   * given in Table 5, From Revised Release on the IAPWS Industrial
   * Formulation 1997 for the Thermodynamic Properties of Water
   * and Steam, IAPWS 2007
   */
  CPPUNIT_TEST( region1 );

  /**
   * Verify calculation of water properties in region 2 using the verification values
   * given in Table 15, From Revised Release on the IAPWS Industrial
   * Formulation 1997 for the Thermodynamic Properties of Water
   * and Steam, IAPWS 2007
   */
  CPPUNIT_TEST( region2 );

  /**
   * Verify calculation of water properties in for the boundary between regions 2 and 3
   * using the verification point (P,T) = (16.5291643 MPa, 623.15 K).
   * From Revised Release on the IAPWS Industrial
   * Formulation 1997 for the Thermodynamic Properties of Water
   * and Steam, IAPWS 2007
   */
  CPPUNIT_TEST( b23 );

  /**
   * Verify calculation of water properties in region 4 (saturation line)
   * using the verification values given in Table 35.
   * From Revised Release on the IAPWS Industrial
   * Formulation 1997 for the Thermodynamic Properties of Water
   * and Steam, IAPWS 2007
   */
  CPPUNIT_TEST( pSat );

  /**
   * Verify calculation of water properties in region 4 (saturation line)
   * using the verification values given in Table 36.
   * From Revised Release on the IAPWS Industrial
   * Formulation 1997 for the Thermodynamic Properties of Water
   * and Steam, IAPWS 2007
   */
  CPPUNIT_TEST( tSat );

  /**
   * Verify calculation of the subregion boundaries in region 3 using the
   * verification values in Table 3 and Table 11.
   * From Revised Supplementary Release on Backward Equations for
   * Specific Volume as a Function of Pressure and Temperature v(p,T)
   * for Region 3 of the IAPWS Industrial Formulation 1997 for the
   * Thermodynamic Properties of Water and Steam
   */
  CPPUNIT_TEST( tempXY );

  /**
   * Verify calculation of the subregion densities in region 3 using the
   * verification values in Table 5 and Table 13.
   * From Revised Supplementary Release on Backward Equations for
   * Specific Volume as a Function of Pressure and Temperature v(p,T)
   * for Region 3 of the IAPWS Industrial Formulation 1997 for the
   * Thermodynamic Properties of Water and Steam
   */
  CPPUNIT_TEST( region3 );

  /**
   * Verify calculation of the water properties in region 5 using the
   * verification values in Table 42 and Table 13.
   * From Revised Release on the IAPWS Industrial
   * Formulation 1997 for the Thermodynamic Properties of Water
   * and Steam, IAPWS 2007
   */
  CPPUNIT_TEST( region5 );

  /**
   * Verify calculation of the water properties in all regions using
   * a sample of date for each region
   */
  CPPUNIT_TEST( density );

  /**
   * Verify viscosity calculation using verification values in Table 4.
   * From Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary
   * Water Substance
   */
  CPPUNIT_TEST( viscosity );

  CPPUNIT_TEST_SUITE_END();

public:
  PorousFlowWaterPropertiesTest();

  void inRegion();
  void region1();
  void region2();
  void b23();
  void pSat();
  void tSat();
  void tempXY();
  void region3();
  void region5();
  void density();
  void viscosity();

 private:
  Real _peps;
  Real _teps;
  Real _deps;
  Real _t_c2k;
};

#endif  // POROUSFLOWWATERPROPERTIESTEST_H
