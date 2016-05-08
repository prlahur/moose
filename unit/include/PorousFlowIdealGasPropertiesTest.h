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

#ifndef POROUSFLOWIDEALGASPROPERTIESTEST_H
#define POROUSFLOWIDEALGASPROPERTIESTEST_H

//CPPUnit includes
#include "GuardedHelperMacros.h"

// Moose includes
#include "PorousFlowIdealGasProperties.h"

class PorousFlowIdealGasPropertiesTest : public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE( PorousFlowIdealGasPropertiesTest );

  /**
   * Verify calculation of ideal gas density and derivatives wrt pressure and
   * temperature
   */
  CPPUNIT_TEST( density );

  CPPUNIT_TEST_SUITE_END();

public:
  PorousFlowIdealGasPropertiesTest();

  void density();

 private:
  Real _peps;
  Real _teps;
};

#endif  // POROUSFLOWIDEALGASPROPERTIESTEST_H
