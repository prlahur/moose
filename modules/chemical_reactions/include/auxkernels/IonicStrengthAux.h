/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef IONICSTRENGTHAUX_H
#define IONICSTRENGTHAUX_H

#include "AuxKernel.h"

class IonicStrengthAux;

template <>
InputParameters validParams<IonicStrengthAux>();

/**
 * Calculates the ionic strength of a solution
 *
 * I = 1/2 \sum_{i = 1}^{N} z_i^2 C_i
 *
 * where N is the total number of aqueous species (primary and secondary),
 * z_i and C_i are the charge and molar concentration of the ith species, respectively
 */
class IonicStrengthAux : public AuxKernel
{
public:
  IonicStrengthAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Charge of chemical species
  const std::vector<Real> _z;
  /// Concentration of chemical species
  std::vector<const VariableValue *> _conc;
};

#endif // IONICSTRENGTHAUX_H
