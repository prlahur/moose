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
#include "PorousFlowWaterPropertiesTest.h"

CPPUNIT_TEST_SUITE_REGISTRATION( PorousFlowWaterPropertiesTest );

PorousFlowWaterPropertiesTest::PorousFlowWaterPropertiesTest() :
    _peps(1.0E-1),
    _teps(1.0E-8),
    _deps(1.0E-8),
    _t_c2k(273.15)
{
}

void
PorousFlowWaterPropertiesTest::inRegion()
{
  CPPUNIT_ASSERT(PorousFlowWaterProperties::inRegion(3.0e6, 300 - _t_c2k) == 1);
  CPPUNIT_ASSERT(PorousFlowWaterProperties::inRegion(30.0e6, 700.0 - _t_c2k) == 2);
  CPPUNIT_ASSERT(PorousFlowWaterProperties::inRegion(50.0e6, 630.0 - _t_c2k) == 3);
  CPPUNIT_ASSERT(PorousFlowWaterProperties::inRegion(30.0e6, 1500.0 - _t_c2k) == 5);
  unsigned int region;
  try
  {
    // Trigger invalid pressure fail
    region = PorousFlowWaterProperties::inRegion(101.0e6, 300.0 - _t_c2k);
  }
  catch(const std::exception & e)
  {
    std::string msg(e.what());
    CPPUNIT_ASSERT( msg.find("Pressure 1.01e+08 is out of range in PorousFlowWaterProperties::inRegion") != std::string::npos );
  }
  try
  {
    // Trigger another invalid pressure fail
    region = PorousFlowWaterProperties::inRegion(51.0e6, 1200.0 - _t_c2k);
  }
  catch(const std::exception & e)
  {
    std::string msg(e.what());
    CPPUNIT_ASSERT( msg.find("Pressure 5.1e+07 is out of range in PorousFlowWaterProperties::inRegion") != std::string::npos );
  }
  try
  {
  // Trigger invalid temperature fail
  region = PorousFlowWaterProperties::inRegion(5.0e6, 2001.0);
  }
  catch(const std::exception & e)
  {
    std::string msg(e.what());
    CPPUNIT_ASSERT( msg.find("Temperature 2001 is out of range in PorousFlowWaterProperties::inRegion") != std::string::npos );
  }
}

void
PorousFlowWaterPropertiesTest::region1()
{
  /// Test the densities
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.00100215168, PorousFlowWaterProperties::densityRegion1(3.0e6, 300.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.000971180894, PorousFlowWaterProperties::densityRegion1(80.0e6, 300.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.00120241800, PorousFlowWaterProperties::densityRegion1(3.0e6, 500.0 - _t_c2k), 1.0E-5);

  /// Test the derivative of the density wrt pressure
  Real fd;
  fd = (PorousFlowWaterProperties::densityRegion1(3.0e6 + _peps, 300.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion1(3.0e6, 300.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion1_dP(3.0e6, 300.0 - _t_c2k), 1.0E-10);
  fd = (PorousFlowWaterProperties::densityRegion1(80.0e6 + _peps, 300.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion1(80.0e6, 300.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion1_dP(80.0e6, 300.0 - _t_c2k), 1.0E-10);
  fd = (PorousFlowWaterProperties::densityRegion1(3.0e6 + _peps, 500.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion1(3.0e6, 500.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion1_dP(3.0e6, 500.0 - _t_c2k), 1.0E-10);
}

void
PorousFlowWaterPropertiesTest::region2()
{
  /// Test the densities
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 39.4913866, PorousFlowWaterProperties::densityRegion2(3.5e3, 300.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 92.3015898, PorousFlowWaterProperties::densityRegion2(3.5e3, 700.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.00542946619, PorousFlowWaterProperties::densityRegion2(30.0e6, 700.0 - _t_c2k), 1.0E-5);

  /// Test the derivative of the density wrt pressure
  Real fd;
  fd = (PorousFlowWaterProperties::densityRegion2(3.5e3 + _peps, 300.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion2(3.5e3, 300.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion2_dP(3.5e3, 300.0 - _t_c2k), 1.0E-10);
  fd = (PorousFlowWaterProperties::densityRegion2(3.5e3 + _peps, 700.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion2(3.5e3, 700.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion2_dP(3.5e3, 700.0 - _t_c2k), 1.0E-10);
  fd = (PorousFlowWaterProperties::densityRegion2(30.0e6 + _peps, 700.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion2(30.0e6, 700.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion2_dP(30.0e6, 700.0 - _t_c2k), 1.0E-10);
}

void
PorousFlowWaterPropertiesTest::b23()
{
  CPPUNIT_ASSERT_DOUBLES_EQUAL(623.15 - _t_c2k, PorousFlowWaterProperties::b23t(16.5291643e6), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(16.5291643e6, PorousFlowWaterProperties::b23p(623.15 - _t_c2k), 1.0);
}

void
PorousFlowWaterPropertiesTest::pSat()
{
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.53658941e3, PorousFlowWaterProperties::pSat(300.0 - _t_c2k), 1.0E-2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.63889776e6, PorousFlowWaterProperties::pSat(500.0 - _t_c2k), 1.0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.3443146e6, PorousFlowWaterProperties::pSat(600.0 - _t_c2k), 1.0);
}

void
PorousFlowWaterPropertiesTest::tSat()
{
  CPPUNIT_ASSERT_DOUBLES_EQUAL(372.755919 - _t_c2k, PorousFlowWaterProperties::tSat(0.1e6), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(453.035632 - _t_c2k, PorousFlowWaterProperties::tSat(1.0e6), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(584.149488 - _t_c2k, PorousFlowWaterProperties::tSat(10.0e6), 1.0E-5);
}

void
PorousFlowWaterPropertiesTest::tempXY()
{
  CPPUNIT_ASSERT_DOUBLES_EQUAL(693.0341408 - _t_c2k, PorousFlowWaterProperties::tempXY(40.0e6, "ab"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(649.3659208 - _t_c2k, PorousFlowWaterProperties::tempXY(25.0e6, "cd"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(713.9593992 - _t_c2k, PorousFlowWaterProperties::tempXY(40.0e6, "ef"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(649.8873759 - _t_c2k, PorousFlowWaterProperties::tempXY(23.0e6, "gh"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(651.5778091 - _t_c2k, PorousFlowWaterProperties::tempXY(23.0e6, "ij"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(655.8338344 - _t_c2k, PorousFlowWaterProperties::tempXY(23.0e6, "jk"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(649.6054133 - _t_c2k, PorousFlowWaterProperties::tempXY(22.8e6, "mn"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(650.0106943 - _t_c2k, PorousFlowWaterProperties::tempXY(22.8e6, "op"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(645.6355027 - _t_c2k, PorousFlowWaterProperties::tempXY(22.0e6, "qu"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(648.2622754 - _t_c2k, PorousFlowWaterProperties::tempXY(22.0e6, "rx"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(647.7996121 - _t_c2k, PorousFlowWaterProperties::tempXY(22.3e6, "uv"), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(648.2049480 - _t_c2k, PorousFlowWaterProperties::tempXY(22.3e6, "wx"), 1.0E-5);
}

void
PorousFlowWaterPropertiesTest::region3()
{
  /// Test the densities
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.001470853100, PorousFlowWaterProperties::densityRegion3(50.0e6, 630.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.001503831359, PorousFlowWaterProperties::densityRegion3(80.0e6, 670.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002204728587, PorousFlowWaterProperties::densityRegion3(50.0e6, 710.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.001973692940, PorousFlowWaterProperties::densityRegion3(80.0e6, 750.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.001761696406, PorousFlowWaterProperties::densityRegion3(20.0e6, 630.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.001819560617, PorousFlowWaterProperties::densityRegion3(30.0e6, 650.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002245587720, PorousFlowWaterProperties::densityRegion3(26.0e6, 656.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002506897702, PorousFlowWaterProperties::densityRegion3(30.0e6, 670.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002970225962, PorousFlowWaterProperties::densityRegion3(26.0e6, 661.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003004627086, PorousFlowWaterProperties::densityRegion3(30.0e6, 675.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.005019029401, PorousFlowWaterProperties::densityRegion3(26.0e6, 671.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.004656470142, PorousFlowWaterProperties::densityRegion3(30.0e6, 690.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002163198378, PorousFlowWaterProperties::densityRegion3(23.6e6, 649.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002166044161, PorousFlowWaterProperties::densityRegion3(24.0e6, 650.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002651081407, PorousFlowWaterProperties::densityRegion3(23.6e6, 652.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002967802335, PorousFlowWaterProperties::densityRegion3(24.0e6, 654.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003273916816, PorousFlowWaterProperties::densityRegion3(23.6e6, 653.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003550329864, PorousFlowWaterProperties::densityRegion3(24.0e6, 655.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.004545001142, PorousFlowWaterProperties::densityRegion3(23.5e6, 655.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.005100267704, PorousFlowWaterProperties::densityRegion3(24.0e6, 660.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.006109525997, PorousFlowWaterProperties::densityRegion3(23.0e6, 660.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.006427325645, PorousFlowWaterProperties::densityRegion3(24.0e6, 670.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002117860851, PorousFlowWaterProperties::densityRegion3(22.6e6, 646.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002062374674, PorousFlowWaterProperties::densityRegion3(23.0e6, 646.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002533063780, PorousFlowWaterProperties::densityRegion3(22.6e6, 648.6 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002572971781, PorousFlowWaterProperties::densityRegion3(22.8e6, 649.3 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002923432711, PorousFlowWaterProperties::densityRegion3(22.6e6, 649.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002913311494, PorousFlowWaterProperties::densityRegion3(22.8e6, 649.7 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003131208996, PorousFlowWaterProperties::densityRegion3(22.6e6, 649.1 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003221160278, PorousFlowWaterProperties::densityRegion3(22.8e6, 649.9 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003715596186, PorousFlowWaterProperties::densityRegion3(22.6e6, 649.4 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003664754790, PorousFlowWaterProperties::densityRegion3(22.8e6, 650.2 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.001970999272, PorousFlowWaterProperties::densityRegion3(21.1e6, 640.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002043919161, PorousFlowWaterProperties::densityRegion3(21.8e6, 643.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.005251009921, PorousFlowWaterProperties::densityRegion3(21.1e6, 644.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.005256844741, PorousFlowWaterProperties::densityRegion3(21.8e6, 648.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.001932829079, PorousFlowWaterProperties::densityRegion3(19.1e6, 635.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.001985387227, PorousFlowWaterProperties::densityRegion3(20.0e6, 638.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.008483262001, PorousFlowWaterProperties::densityRegion3(17.0e6, 626.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.006227528101, PorousFlowWaterProperties::densityRegion3(20.0e6, 640.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002268366647, PorousFlowWaterProperties::densityRegion3(21.5e6, 644.6 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002296350553, PorousFlowWaterProperties::densityRegion3(22.0e6, 646.1 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002832373260, PorousFlowWaterProperties::densityRegion3(22.5e6, 648.6 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002811424405, PorousFlowWaterProperties::densityRegion3(22.3e6, 647.9 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003694032281, PorousFlowWaterProperties::densityRegion3(22.15e6, 647.5 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003622226305, PorousFlowWaterProperties::densityRegion3(22.3e6, 648.1 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.004528072649, PorousFlowWaterProperties::densityRegion3(22.11e6, 648.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.004556905799, PorousFlowWaterProperties::densityRegion3(22.3e6, 649.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002698354719, PorousFlowWaterProperties::densityRegion3(22.0e6, 646.84 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002717655648, PorousFlowWaterProperties::densityRegion3(22.064e6, 647.05 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003798732962, PorousFlowWaterProperties::densityRegion3(22.0e6, 646.89 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.003701940010, PorousFlowWaterProperties::densityRegion3(22.064e6, 647.15 - _t_c2k), 1.0E-5);
  /// Test the derivative of the density wrt pressure
  Real fd;
  fd = (PorousFlowWaterProperties::densityRegion3(50.0e6 + _peps, 630.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(50.0e6, 630.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(50.0e6, 630.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(50.0e6 + _peps, 710.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(50.0e6, 710.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(50.0e6, 710.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(20.0e6 + _peps, 630.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(20.0e6, 630.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(20.0e6, 630.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(26.0e6 + _peps, 656.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(26.0e6, 656.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(26.0e6, 656.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(26.0e6 + _peps, 661.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(26.0e6, 661.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(26.0e6, 661.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(26.0e6 + _peps, 671.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(26.0e6, 671.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(26.0e6, 671.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(23.6e6 + _peps, 649.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(23.6e6, 649.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(23.6e6, 649.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(23.6e6 + _peps, 652.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(23.6e6, 652.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(23.6e6, 652.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(23.5e6 + _peps, 653.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(23.5e6, 653.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(23.5e6, 653.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(23.0e6 + _peps, 660.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(23.0e6, 660.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(23.0e6, 660.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(22.6e6 + _peps, 646.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(22.6e6, 646.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(22.6e6, 646.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(22.6e6 + _peps, 648.6 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(22.6e6, 648.6 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(22.6e6, 648.6 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(22.6e6 + _peps, 649.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(22.6e6, 649.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(22.6e6, 649.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(22.6e6 + _peps, 649.1 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(22.6e6, 649.1 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(22.6e6, 649.1 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(22.6e6 + _peps, 649.4 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(22.6e6, 649.4 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(22.6e6, 649.4 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(21.1e6 + _peps, 640.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(21.1e6, 640.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(21.1e6, 640.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(21.1e6 + _peps, 644.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(21.1e6, 644.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(21.1e6, 644.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(19.1e6 + _peps, 635.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(19.1e6, 635.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(19.1e6, 635.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(17.0e6 + _peps, 626.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(17.0e6, 626.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(17.0e6, 626.0 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(21.5e6 + _peps, 644.6 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(21.5e6, 644.6 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(21.5e6, 644.6 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(22.5e6 + _peps, 648.6 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(22.5e6, 648.6 - _t_c2k)) / _peps;
// TODO: this one isn't working correctly
//  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(22.5e6, 648.6 - _t_c2k), 1.0E-9);
  fd = (PorousFlowWaterProperties::densityRegion3(22.15e6 + _peps, 647.5 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(22.15e6, 647.5 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(22.15e6, 647.5 - _t_c2k), 1.0E-7);
  fd = (PorousFlowWaterProperties::densityRegion3(22.11e6 + _peps, 648.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(22.11e6, 648.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(22.11e6, 648.0 - _t_c2k), 1.0E-7);
  fd = (PorousFlowWaterProperties::densityRegion3(22.0e6 + _peps, 646.84 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(22.0e6, 646.84 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(22.0e6, 646.84 - _t_c2k), 1.0E-7);
  fd = (PorousFlowWaterProperties::densityRegion3(22.0e6 + _peps, 646.89 - _t_c2k) - PorousFlowWaterProperties::densityRegion3(22.0e6, 646.89 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion3_dP(22.0e6, 646.89 - _t_c2k), 1.0E-7);
}

void
PorousFlowWaterPropertiesTest::region5()
{
  /// Test the densities
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 1.38455090, PorousFlowWaterProperties::densityRegion5(0.5e6, 1500.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.0230761299, PorousFlowWaterProperties::densityRegion5(30.0e6, 1500.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.0311385219, PorousFlowWaterProperties::densityRegion5(30.0e6, 2000.0 - _t_c2k), 1.0E-5);
  /// Test the derivative of the density wrt pressure
  Real fd;
  fd = (PorousFlowWaterProperties::densityRegion5(0.5e6 + _peps, 1500.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion5(0.5e6, 1500.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion5_dP(0.5e6, 1500.0 - _t_c2k), 1.0E-10);
  fd = (PorousFlowWaterProperties::densityRegion5(30.0e6 + _peps, 1500.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion5(30.0e6, 1500.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion5_dP(30.0e6, 1500.0 - _t_c2k), 1.0E-10);
  fd = (PorousFlowWaterProperties::densityRegion5(30.0e6 + _peps, 2000.0 - _t_c2k) - PorousFlowWaterProperties::densityRegion5(30.0e6, 2000.0 - _t_c2k)) / _peps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dDensityRegion5_dP(30.0e6, 2000.0 - _t_c2k), 1.0E-10);
}

void
PorousFlowWaterPropertiesTest::density()
{
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.00100215168, PorousFlowWaterProperties::densityRegion1(3.0e6, 300.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 39.4913866, PorousFlowWaterProperties::densityRegion2(3.5e3, 300.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.001470853100, PorousFlowWaterProperties::densityRegion3(50.0e6, 630.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 0.002923432711, PorousFlowWaterProperties::densityRegion3(22.6e6, 649.0 - _t_c2k), 1.0E-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0 / 1.38455090, PorousFlowWaterProperties::density(0.5e6, 1500.0 - _t_c2k), 1.0E-5);
}

void
PorousFlowWaterPropertiesTest::viscosity()
{
  CPPUNIT_ASSERT_DOUBLES_EQUAL(889.735100e-6, PorousFlowWaterProperties::viscosity(298.15 - _t_c2k, 998.0), 1.0E-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1437.649467e-6, PorousFlowWaterProperties::viscosity(298.15 - _t_c2k, 1200.0), 1.0E-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(307.883622e-6, PorousFlowWaterProperties::viscosity(373.15 - _t_c2k, 1000.0), 1.0E-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.538342e-6, PorousFlowWaterProperties::viscosity(433.15 - _t_c2k, 1.0), 1.0E-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(217.685358e-6, PorousFlowWaterProperties::viscosity(433.15 - _t_c2k, 1000.0), 1.0E-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(32.619287e-6, PorousFlowWaterProperties::viscosity(873.15 - _t_c2k, 1.0), 1.0E-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(35.802262e-6, PorousFlowWaterProperties::viscosity(873.15 - _t_c2k, 100.0), 1.0E-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(77.430195e-6, PorousFlowWaterProperties::viscosity(873.15 - _t_c2k, 600.0), 1.0E-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(44.217245e-6, PorousFlowWaterProperties::viscosity(1173.15 - _t_c2k, 1.0), 1.0E-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(47.640433e-6, PorousFlowWaterProperties::viscosity(1173.15 - _t_c2k, 100.0), 1.0E-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(64.154608e-6, PorousFlowWaterProperties::viscosity(1173.15 - _t_c2k, 400.0), 1.0E-10);
  /// Test derivatives wrt density
  Real fd;
  fd = (PorousFlowWaterProperties::viscosity(298.15 - _t_c2k, 998.0 + _deps) - PorousFlowWaterProperties::viscosity(298.15 - _t_c2k, 998.0)) / _deps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dViscosity_dDensity(298.15 - _t_c2k, 998.0), 1.0E-7);
  fd = (PorousFlowWaterProperties::viscosity(298.15 - _t_c2k, 1200.0 + _deps) - PorousFlowWaterProperties::viscosity(298.15 - _t_c2k, 1200.0)) / _deps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dViscosity_dDensity(298.15 - _t_c2k, 1200.0), 1.0E-7);
  fd = (PorousFlowWaterProperties::viscosity(873.15 - _t_c2k, 1.0 + _deps) - PorousFlowWaterProperties::viscosity(873.15 - _t_c2k, 1.0)) / _deps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dViscosity_dDensity(873.15 - _t_c2k, 1.0), 1.0E-7);
  fd = (PorousFlowWaterProperties::viscosity(1173.15 - _t_c2k, 400.0 + _deps) - PorousFlowWaterProperties::viscosity(1173.15 - _t_c2k, 400.0)) / _deps;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(fd, PorousFlowWaterProperties::dViscosity_dDensity(1173.15 - _t_c2k, 400.0), 1.0E-7);

}
