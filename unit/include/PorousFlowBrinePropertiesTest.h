/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWBRINEPROPERTIESTEST_H
#define POROUSFLOWBRINEPROPERTIESTEST_H

//CPPUnit includes
#include "GuardedHelperMacros.h"

// Moose includes
#include "PorousFlowBrineProperties.h"

class PorousFlowBrinePropertiesTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( PorousFlowBrinePropertiesTest );

  /**
   * Verify calculation of the brine density using values provided in the
   * supplementary data from Driesner, The system H2o-NaCl. Part II:
   *  Correlations for molar volume, enthalpy, and isobaric heat capacity
   * from 0 to 1000 C, 1 to 500 bar, and 0 to 1 Xnacl,
   * Geochimica et Cosmochimica Acta 71, 4902-4919 (2007).
   *
   * Note: Driesner uses an outdated version of IAPWS to calculate water properties,
   * hence there is a small difference in the results, so the results are compared to
   * a tolerance of +- 0.1%.
   */
  CPPUNIT_TEST( density );

  /**
   * Verify the saturation pressure using experimental data from
   *  Liu and Lindsay, Thermodynamics of sodium chloride solutions at high temperatures,
   *  Journal of Solution Chemistry, 1, 45-69 (1972)
   */
  CPPUNIT_TEST ( pSat );

  /**
   * Verify calculation of the brine viscosity using experimental values provided in
   * Phillips et al, Viscosity of NaCl and other solutions up to 350C and 50MPa pressures,
   * LBL-11586 (1980)
   */
  CPPUNIT_TEST( viscosity );

  /**
   * Verify calcualtion of halite solubility using experimental data from
   * Bodnar et al, Synthetic fluid inclusions in natural quartz, III. Determination of
   * phase equilibrium properties in the system H2O-NaCl to 1000C and 1500 bars,
   * Geocehmica et Cosmochemica Acta, 49, 1861-1873 (1985) and re-interpreted in
   * Chou (1987).
   *
   * Note that the solubility is in %, and that the average of the range quoted has
   * been used for each point.
   */
  CPPUNIT_TEST ( solubility );

  CPPUNIT_TEST_SUITE_END();

public:
  PorousFlowBrinePropertiesTest();

  void density();
  void pSat();
  void viscosity();
  void solubility();

  /**
   * Function to convert mole fraction to mass fraction
   */
  Real mol2Mass(Real xmol);

  /**
   * Function to convert molality to mass fraction
   */
  Real molality2Mass(Real xmol);

 private:
  Real _eps;
  Real _Mnacl;
  Real _Mh2o;
  Real _t_c2k;
};

#endif  // POROUSFLOWBRINEPROPERTIESTEST_H
