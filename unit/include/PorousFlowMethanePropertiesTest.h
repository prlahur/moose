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

#ifndef POROUSFLOWMETHANEPROPERTIESTEST_H
#define POROUSFLOWMETHANEPROPERTIESTEST_H

//CPPUnit includes
#include "GuardedHelperMacros.h"

// Moose includes
#include "PorousFlowMethaneProperties.h"

class PorousFlowMethanePropertiesTest : public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE( PorousFlowMethanePropertiesTest );

  /**
   * Verify calculation of methane density (derivatives have been verified in the
   * ideal gas unit tests). Compare results with data from NIST webbook (www.nist.gov)
   */
  CPPUNIT_TEST( density );

  /**
   * Verify calculation of methane viscosity and derivatives wrt density and
   * temperature. Data from Table A VII.17, from Irvine Jr, T. F. and
   *  Liley, P. E. (1984) Steam and Gas Tables with Computer Equations. 
   */
  CPPUNIT_TEST( viscosity );

  CPPUNIT_TEST_SUITE_END();

public:
  PorousFlowMethanePropertiesTest();

  void density();
  void viscosity();

 private:
  Real _eps;
  Real _t_c2k;
};

#endif  // POROUSFLOWMETHANEPROPERTIESTEST_H
