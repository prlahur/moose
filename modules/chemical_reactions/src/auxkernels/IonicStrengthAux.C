/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "IonicStrengthAux.h"

template <>
InputParameters
validParams<IonicStrengthAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<std::vector<Real>>("z", "The charge of each chemical species");
  params.addCoupledVar("conc", "The aqueous species in solution (primary and secondary)");
  MooseEnum units("molarity molality", "molarity");
  params.addParam<MooseEnum>("units", units, "Units of ionic strength (molarity or molality)");
  params.addCoupledVar("density", 1000.0, "Density of water");
  params.addClassDescription("Ionic strength of solution");
  return params;
}

IonicStrengthAux::IonicStrengthAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _z(getParam<std::vector<Real>>("z")),
    _units(getParam<MooseEnum>("units").getEnum<IonicStrengthUnitsEnum>()),
    _density(coupledValue("density"))
{
  const unsigned int n = coupledComponents("conc");

  // Check that n charges have been input
  if (_z.size() != n)
    mooseError("The number of charges and species concentrations must be equal in ", _name);

  _conc.resize(n);
  for (unsigned int i = 0; i < n; ++i)
    _conc[i] = &coupledValue("conc", i);
}

Real
IonicStrengthAux::computeValue()
{
  Real sum = 0.0;

  for (unsigned int i = 0; i < _conc.size(); ++i)
    sum += _z[i] * _z[i] * (*_conc[i])[_qp];

  if (_units == IonicStrengthUnitsEnum::MOLALITY)
    sum /= _density[_qp];

  return 0.5 * sum;
}
