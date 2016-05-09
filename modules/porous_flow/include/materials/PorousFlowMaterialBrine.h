/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWMATERIALBRINE_H
#define POROUSFLOWMATERIALBRINE_H

#include "PorousFlowMaterialFluidPropertiesBase.h"
#include "PorousFlowBrineProperties.h"

class PorousFlowMaterialBrine;

template<>
InputParameters validParams<PorousFlowMaterialBrine>();

/**
 * Fluid properties of brine.
 * Provides density, viscosity and derivatives wrt pressure and temperature.
 */
class PorousFlowMaterialBrine : public PorousFlowMaterialFluidPropertiesBase
{
public:
  PorousFlowMaterialBrine(const InputParameters & parameters);

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

#endif //POROUSFLOWMATERIALBRINE_H
