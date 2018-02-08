//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef EQUILIBRIUMCONSTANTAUX_H
#define EQUILIBRIUMCONSTANTAUX_H

#include "AuxKernel.h"
#include "EquilibriumConstantFit.h"

class EquilibriumConstantAux;

template <>
InputParameters validParams<EquilibriumConstantAux>();

/**
 * Equilibrium constant (in the form log10(K)) calculated using a least-squares
 * fit to the data provided (typically taken from a geochemical database).
 *
 * Fitted function is a Maier-Kelley type function for the equilibrium constant
 *
 * log(K)= a_0 ln(T) + a_1 + a_2 T + a_3 / T + a_4 / T^2
 *
 * where T is the temperature in Kelvin
 */
class EquilibriumConstantAux : public AuxKernel
{
public:
  EquilibriumConstantAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Temperature (in K)
  const VariableValue & _temperature;
  /// Temperature points (in K)
  const std::vector<Real> & _temperature_points;
  /// log(K) values at each temperature point
  const std::vector<Real> & _logk_points;
  /// Least-squares fit to data
  std::unique_ptr<EquilibriumConstantFit> _logk;
};

#endif // EQUILIBRIUMCONSTANTAUX_H
