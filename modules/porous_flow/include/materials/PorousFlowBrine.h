/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWBRINE_H
#define POROUSFLOWBRINE_H

#include "PorousFlowFluidPropertiesBase.h"
#include "PorousFlowBrineProperties.h"

class PorousFlowBrine;

template<>
InputParameters validParams<PorousFlowBrine>();

/**
 * Fluid properties of brine.
 * Provides density, viscosity and derivatives wrt pressure and temperature.
 */
class PorousFlowBrine : public PorousFlowFluidPropertiesBase
{
public:
  PorousFlowBrine(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();

  virtual void computeQpProperties();

  /// Fluid phase density at the nodes
  MaterialProperty<Real> & _density_nodal;

  /// Old fluid phase density at the nodes
  MaterialProperty<Real> & _density_nodal_old;

  /// Derivative of fluid density wrt phase pore pressure at the nodes
  MaterialProperty<Real> & _ddensity_nodal_dp;

  /// Fluid phase density at the qps
  MaterialProperty<Real> & _density_qp;

  /// Derivative of fluid density wrt phase pore pressure at the qps
  MaterialProperty<Real> & _ddensity_qp_dp;
  
  /// Fluid phase viscosity at the nodes
  MaterialProperty<Real> & _viscosity_nodal;

  /// NaCl mass fraction
  Real _xnacl;
};

#endif //POROUSFLOWBRINE_H
