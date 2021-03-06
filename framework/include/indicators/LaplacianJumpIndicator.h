/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef LAPLACIANJUMPINDICATOR_H
#define LAPLACIANJUMPINDICATOR_H

#include "JumpIndicator.h"

class LaplacianJumpIndicator;

template<>
InputParameters validParams<LaplacianJumpIndicator>();

class LaplacianJumpIndicator :
  public JumpIndicator
{
public:
  LaplacianJumpIndicator(const InputParameters & parameters);
  virtual ~LaplacianJumpIndicator(){};

protected:

  virtual Real computeQpIntegral();

  const VariableSecond & _second_u;
  const VariableSecond & _second_u_neighbor;
};

#endif /* LAPLACIANJUMPINDICATOR_H */
