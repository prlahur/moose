/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWIDEALGAS_H
#define POROUSFLOWIDEALGAS_H

#include "PorousFlowFluidPropertiesBase.h"
#include "PorousFlowIdealGasProperties.h"

class PorousFlowIdealGas;

template<>
InputParameters validParams<PorousFlowIdealGas>();

/**
 * This material computes fluid properties for an ideal gas
 */
class PorousFlowIdealGas : public PorousFlowFluidPropertiesBase
{
public:
  PorousFlowIdealGas(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();

  virtual void computeQpProperties();

  /// Molar mass (kg/mol)
  const Real _molar_mass;

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
};

#endif //POROUSFLOWIDEALGAS_H
