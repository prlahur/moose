//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COUPLEDBEKINETIC_H
#define COUPLEDBEKINETIC_H

#include "Kernel.h"

// Forward Declarations
class CoupledBEKinetic;

template <>
InputParameters validParams<CoupledBEKinetic>();

/**
 * Rate of primary species consumed in mineral reaction
 */
class CoupledBEKinetic : public Kernel
{
public:
  CoupledBEKinetic(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:
  /// Porosity
  const MaterialProperty<Real> & _porosity;
  /// Old value of porosity
  const MaterialProperty<Real> & _porosity_old;
  /// Weight of the kinetic mineral concentration in the total primary species concentration
  const std::vector<Real> _weight;
  /// Coupled kinetic mineral concentrations
  std::vector<const VariableValue *> _vals;
  /// Coupled old values of kinetic mineral concentrations
  std::vector<const VariableValue *> _vals_old;
};

#endif // COUPLEDBEKINETIC_H
