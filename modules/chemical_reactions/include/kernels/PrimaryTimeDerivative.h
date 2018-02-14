//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PRIMARYTIMEDERIVATIVE
#define PRIMARYTIMEDERIVATIVE

#include "TimeDerivative.h"
#include "DerivativeMaterialInterface.h"

class PrimaryTimeDerivative;

template <>
InputParameters validParams<PrimaryTimeDerivative>();

/**
 * Derivative of primary species concentration with respect to time
 */
class PrimaryTimeDerivative : public DerivativeMaterialInterface<TimeDerivative>
{
public:
  PrimaryTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Porosity
  const MaterialProperty<Real> & _porosity;
  /// Old value of porosity
  const MaterialProperty<Real> & _porosity_old;
  /// Old value of the primary species concentration
  const VariableValue & _u_old;
};

#endif // PRIMARYTIMEDERIVATIVE
