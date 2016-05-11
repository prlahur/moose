/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWMETHANE_H
#define POROUSFLOWMETHANE_H

#include "PorousFlowFluidPropertiesBase.h"
#include "PorousFlowMethaneProperties.h"

class PorousFlowMethane;

template<>
InputParameters validParams<PorousFlowMethane>();

/**
 * Fluid properties of Methane (CH4).
 * Provides density, viscosity, derivatives wrt pressure and temperature,
 * and Henry Law constants.
 */
class PorousFlowMethane : public PorousFlowFluidPropertiesBase
{
public:
  PorousFlowMethane(const InputParameters & parameters);

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

  /// Derivative of fluid phase viscosity wrt temperature at the nodes
  MaterialProperty<Real> & _dviscosity_nodal_dt;
  
  /// Methane molar mass (kg/mol)
  const Real _Mch4;
};

#endif //POROUSFLOWMETHANE_H
