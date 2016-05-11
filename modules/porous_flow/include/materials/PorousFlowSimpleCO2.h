/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWSIMPLECO2_H
#define POROUSFLOWSIMPLECO2_H

#include "PorousFlowFluidPropertiesBase.h"
#include "PorousFlowSimpleCO2Properties.h"

class PorousFlowSimpleCO2;

template<>
InputParameters validParams<PorousFlowSimpleCO2>();

/**
 * Simplified fluid properties of carbon dioxide (CO2).
 * Provides density, viscosity, derivatives wrt pressure and temperature,
 * and Henry Law constants.
 *
 * Note: These aren't as accurate as the full CO2 implementation, but are
 * fast and have simple analytical derivatives
 */
class PorousFlowSimpleCO2 : public PorousFlowFluidPropertiesBase
{
public:
  PorousFlowSimpleCO2(const InputParameters & parameters);

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

#endif //POROUSFLOWSIMPLECO2_H
