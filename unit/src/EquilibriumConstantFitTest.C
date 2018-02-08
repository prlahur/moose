//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "EquilibriumConstantFit.h"

#include <cmath>

const double tol = 1.0e-2;
const std::vector<double> x = {273.15, 298.15, 333.15, 373.15, 423.15, 473.15, 523.15, 573.15};
const std::vector<double> y = {
    14.9398, 13.9951, 13.0272, 12.2551, 11.6308, 11.2836, 11.1675, 11.3002};

TEST(EquilibriumConstantFitTest, constructor)
{
  EquilibriumConstantFit logk(x, y);
  EXPECT_EQ(logk.getSampleSize(), x.size());
}

TEST(EquilibriumConstantFitTest, sample)
{
  EquilibriumConstantFit logk(x, y);
  logk.generate();

  EXPECT_NEAR(logk.sample(x[1]), y[1], tol);
  EXPECT_NEAR(logk.sample(x[2]), y[2], tol);
  EXPECT_NEAR(logk.sample(x[3]), y[3], tol);
  EXPECT_NEAR(logk.sample(x[4]), y[4], tol);
  EXPECT_NEAR(logk.sample(x[5]), y[5], tol);
}
