//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PrimaryTimeDerivative.h"

template <>
InputParameters
validParams<PrimaryTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addParam<MaterialPropertyName>(
      "porosity_name", "porosity", "The name of the material property defining porosity");
  params.addClassDescription("Derivative of primary species concentration wrt time");
  return params;
}

PrimaryTimeDerivative::PrimaryTimeDerivative(const InputParameters & parameters)
  : DerivativeMaterialInterface<TimeDerivative>(parameters),
    _porosity_name(getParam<MaterialPropertyName>("porosity_name")),
    _porosity(getMaterialProperty<Real>(_porosity_name)),
    _porosity_old(getMaterialPropertyOld<Real>(_porosity_name)),
    _u_old(valueOld())
{
}

Real
PrimaryTimeDerivative::computeQpResidual()
{
  return _test[_i][_qp] * (_porosity[_qp] * _u[_qp] - _porosity_old[_qp] * _u_old[_qp]) / _dt;
}

Real
PrimaryTimeDerivative::computeQpJacobian()
{
  return _test[_i][_qp] * _phi[_j][_qp] * _porosity[_qp] / _dt;
}
