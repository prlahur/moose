/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowWater.h"
#include "Conversion.h"

template<>
InputParameters validParams<PorousFlowWater>()
{
  InputParameters params = validParams<PorousFlowFluidPropertiesBase>();
  params.addClassDescription("This Material calculates fluid properties for water (H2O)");
  return params;
}

PorousFlowWater::PorousFlowWater(const InputParameters & parameters) :
    PorousFlowFluidPropertiesBase(parameters),

  _density_nodal(declareProperty<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num))),
  _density_nodal_old(declarePropertyOld<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num))),
  _ddensity_nodal_dp(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num), _pressure_variable_name)),
  _density_qp(declareProperty<Real>("PorousFlow_fluid_phase_density_qp" + Moose::stringify(_phase_num))),
  _ddensity_qp_dp(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density_qp" + Moose::stringify(_phase_num), _pressure_variable_name)),
  _viscosity_nodal(declareProperty<Real>("PorousFlow_viscosity" + Moose::stringify(_phase_num)))
{
}

void
PorousFlowWater::initQpStatefulProperties()
{
  _density_nodal[_qp] = PorousFlowWaterProperties::density(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num]);
}

void
PorousFlowWater::computeQpProperties()
{
  /// Density and derivatives wrt pressure and temperature at the nodes
  _density_nodal[_qp] = PorousFlowWaterProperties::density(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num]);
  _ddensity_nodal_dp[_qp] = PorousFlowWaterProperties::dDensity_dP(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num]);

  /// Density and derivatives wrt pressure and temperature at the qps
  _density_qp[_qp] = PorousFlowWaterProperties::density(_porepressure_qp[_qp][_phase_num], _temperature_qp[_qp][_phase_num]);
  _ddensity_qp_dp[_qp] = PorousFlowWaterProperties::dDensity_dP(_porepressure_qp[_qp][_phase_num], _temperature_qp[_qp][_phase_num]);

  /// Viscosity and derivative wrt temperature at the nodes
  _viscosity_nodal[_qp] = PorousFlowWaterProperties::viscosity(_temperature_nodal[_qp][_phase_num], _density_nodal[_qp]);
}
