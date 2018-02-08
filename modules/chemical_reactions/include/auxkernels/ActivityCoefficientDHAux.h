/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ACTIVITYCOEFFICIENTDHAUX_H
#define ACTIVITYCOEFFICIENTDHAUX_H

#include "AuxKernel.h"

class ActivityCoefficientDHAux;

template <>
InputParameters validParams<ActivityCoefficientDHAux>();

/**
 * Calculates the activity coefficient of a chemical species, gamma,
 * using the extended Debye-Huckel model
 *
 * ln(gamma) = - (z^2) A sqrt(I) / (1 + a B sqrt(I)) + b I
 *
 * where z is the charge of the species, I is the ionic strength of the solution,
 * A and B are the Debye-Huckel parameters, a is the radius of the species, and b
 * is a constant that includes the effect of ion pair interactions.
 */
class ActivityCoefficientDHAux : public AuxKernel
{
public:
  ActivityCoefficientDHAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Temperature (in K)
  const VariableValue & _temperature;
  /// Charge of chemical species (valence)
  const Real _z;
  /// Ionic strength of solution
  const VariableValue & _ionic_strength;
  /// Effective radius of ion (Angstrom)
  const Real _a;
  /// Debye-Huckel parameter A
  const Real _A;
  /// Debye-Huckel parameter A
  const Real _B;
  /// Coefficient bdot
  const Real _bdot;
};

#endif // ACTIVITYCOEFFICIENTDHAUX_H
