/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWCONSTANTBULKPROPERTIESTEST_H
#define POROUSFLOWCONSTANTBULKPROPERTIESTEST_H

//CPPUnit includes
#include "GuardedHelperMacros.h"

// Moose includes
#include "PorousFlowConstantBulkProperties.h"

class PorousFlowConstantBulkPropertiesTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( PorousFlowConstantBulkPropertiesTest );

  /**
   * Verify calculation of ideal gas density and derivatives wrt pressure
   */
  CPPUNIT_TEST( density );

  CPPUNIT_TEST_SUITE_END();

public:
  PorousFlowConstantBulkPropertiesTest();

  void density();
};

#endif  // POROUSFLOWCONSTANTBULKPROPERTIESTEST_H
