/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWMATERIALWATER_H
#define POROUSFLOWMATERIALWATER_H

#include "PorousFlowMaterialFluidPropertiesBase.h"
#include "PorousFlowWaterProperties.h"

class PorousFlowMaterialWater;

template<>
InputParameters validParams<PorousFlowMaterialWater>();

/**
 * Fluid properties of Water (H20).
 * Provides density, viscosity and derivatives wrt pressure and temperature.
 */
class PorousFlowMaterialWater : public PorousFlowMaterialFluidPropertiesBase
{
public:
  PorousFlowMaterialWater(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  /// Fluid phase density at the nodes
  MaterialProperty<Real> & _density_nodal;
  /// Old fluid phase density at the nodes
  MaterialProperty<Real> & _density_nodal_old;
  /// Derivative of fluid density wrt phase pore pressure at the nodes
  MaterialProperty<Real> & _ddensity_nodal_dp;
  /// Derivative of fluid density wrt temperature at the nodes
  MaterialProperty<Real> & _ddensity_nodal_dt;
  /// Fluid phase density at the qps
  MaterialProperty<Real> & _density_qp;
  /// Derivative of fluid density wrt phase pore pressure at the qps
  MaterialProperty<Real> & _ddensity_qp_dp;
  /// Derivative of fluid density wrt temperature at the qps
  MaterialProperty<Real> & _ddensity_qp_dt;
  /// Fluid phase viscosity at the nodes
  MaterialProperty<Real> & _viscosity_nodal;
};

#endif //POROUSFLOWMATERIALWATER_H
