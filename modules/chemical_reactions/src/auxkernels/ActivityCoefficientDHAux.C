/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ActivityCoefficientDHAux.h"

template <>
InputParameters
validParams<ActivityCoefficientDHAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addCoupledVar(
      "temperature", 298.15, "The temperature of the aqueous phase (K). Default is 298.15K");
  params.addRequiredParam<Real>("z", "The charge of the chemical species");
  params.addRequiredCoupledVar("ionic_strength", "The ionic strength of the solution");
  params.addRequiredParam<Real>("a", "The effective radius of the ion (in Angstrom)");
  const std::vector<Real> a{0.51};
  params.addParam<std::vector<Real>>("A", a, "Debye-Huckel A parameter");
  params.addClassDescription("Activity coefficient using extended Debye-Huckel model");
  return params;
}

ActivityCoefficientDHAux::ActivityCoefficientDHAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _temperature(coupledValue("temperature")),
    _z(getParam<Real>("z")),
    _ionic_strength(coupledValue("ionic_strength")),
    _a(getParam<Real>("a")),
    _A(getParam<std::vector<Real>>("A")),
    _B(0.329),
    _bdot(0.041)
{
}

Real
ActivityCoefficientDHAux::computeValue()
{
  Real sqrtI = std::sqrt(_ionic_strength[_qp]);
  Real log_gamma = -(_z * _z * _A * sqrtI) / (1.0 + _a * _B * sqrtI) + _bdot * _ionic_strength[_qp];

  return std::exp(log_gamma);
}
