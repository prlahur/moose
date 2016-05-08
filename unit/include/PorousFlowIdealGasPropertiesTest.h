/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
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
};

#endif  // POROUSFLOWIDEALGASPROPERTIESTEST_H
