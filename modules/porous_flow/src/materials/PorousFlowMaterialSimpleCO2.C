/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowMaterialSimpleCO2.h"
#include "Conversion.h"

template<>
InputParameters validParams<PorousFlowMaterialSimpleCO2>()
{
  InputParameters params = validParams<PorousFlowMaterialFluidPropertiesBase>();
  params.addClassDescription("This Material calculates fluid properties for CO2");
  return params;
}

PorousFlowMaterialSimpleCO2::PorousFlowMaterialSimpleCO2(const InputParameters & parameters) :
    PorousFlowMaterialFluidPropertiesBase(parameters),

  _density_nodal(declareProperty<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num))),
  _density_nodal_old(declarePropertyOld<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num))),
  _ddensity_nodal_dp(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num), _pressure_variable_name)),
  _ddensity_nodal_dt(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num), _temperature_variable_name)),
  _density_qp(declareProperty<Real>("PorousFlow_fluid_phase_density_qp" + Moose::stringify(_phase_num))),
  _ddensity_qp_dp(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density_qp" + Moose::stringify(_phase_num), _pressure_variable_name)),
  _ddensity_qp_dt(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density_qp" + Moose::stringify(_phase_num), _temperature_variable_name)),
  _viscosity_nodal(declareProperty<Real>("PorousFlow_viscosity" + Moose::stringify(_phase_num)))
{
}

void
PorousFlowMaterialSimpleCO2::initQpStatefulProperties()
{
  _density_nodal[_qp] = PorousFlowSimpleCO2Properties::density(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num]);
}

void
PorousFlowMaterialSimpleCO2::computeQpProperties()
{
  /// Density and derivatives wrt pressure and temperature at the nodes
  _density_nodal[_qp] = PorousFlowSimpleCO2Properties::density(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num]);
  _ddensity_nodal_dp[_qp] = PorousFlowSimpleCO2Properties::dDensity_dP(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num]);
  _ddensity_nodal_dt[_qp] = PorousFlowSimpleCO2Properties::dDensity_dT(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num]);

  /// Density and derivatives wrt pressure and temperature at the qps
  _density_qp[_qp] = PorousFlowSimpleCO2Properties::density(_porepressure_qp[_qp][_phase_num], _temperature_qp[_qp][_phase_num]);
  _ddensity_qp_dp[_qp] = PorousFlowSimpleCO2Properties::dDensity_dP(_porepressure_qp[_qp][_phase_num], _temperature_qp[_qp][_phase_num]);
  _ddensity_qp_dt[_qp] = PorousFlowSimpleCO2Properties::dDensity_dT(_porepressure_qp[_qp][_phase_num], _temperature_qp[_qp][_phase_num]);

  /// Viscosity and derivative wrt temperature at the nodes
  _viscosity_nodal[_qp] = PorousFlowSimpleCO2Properties::viscosity(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num], _density_nodal[_qp]);
}
