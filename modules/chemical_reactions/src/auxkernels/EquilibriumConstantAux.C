//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EquilibriumConstantAux.h"

template <>
InputParameters
validParams<EquilibriumConstantAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("temperature", "The temperature of the aqueous phase");
  params.addParam<std::vector<Real>>("temperature_points",
                                     "Temperature points where log(K) data is evaluated");
  params.addParam<std::vector<Real>>("logk_points",
                                     "log(K) data evaluated at each value of temperature_points");
  params.addClassDescription(
      "Equilibrium constant for a given equilibrium species (in form log10(K))");
  return params;
}

EquilibriumConstantAux::EquilibriumConstantAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _temperature(coupledValue("temperature")),
    _temperature_points(getParam<std::vector<Real>>("temperature_points")),
    _logk_points(getParam<std::vector<Real>>("logk_points"))
{
  // Least-squares fit
  _logk = libmesh_make_unique<EquilibriumConstantFit>(_temperature_points, _logk_points);
  _logk->generate();
}

Real
EquilibriumConstantAux::computeValue()
{
  return _logk->sample(_temperature[_qp]);
}
