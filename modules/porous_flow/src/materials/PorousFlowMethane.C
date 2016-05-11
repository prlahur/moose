/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowMethane.h"
#include "Conversion.h"

template<>
InputParameters validParams<PorousFlowMethane>()
{
  InputParameters params = validParams<PorousFlowFluidPropertiesBase>();
  params.addClassDescription("This Material calculates fluid properties for methane");
  return params;
}

PorousFlowMethane::PorousFlowMethane(const InputParameters & parameters) :
    PorousFlowFluidPropertiesBase(parameters),

  _density_nodal(declareProperty<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num))),
  _density_nodal_old(declarePropertyOld<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num))),
  _ddensity_nodal_dp(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num), _pressure_variable_name)),
  _ddensity_nodal_dt(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num), _temperature_variable_name)),
  _density_qp(declareProperty<Real>("PorousFlow_fluid_phase_density_qp" + Moose::stringify(_phase_num))),
  _ddensity_qp_dp(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density_qp" + Moose::stringify(_phase_num), _pressure_variable_name)),
  _ddensity_qp_dt(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density_qp" + Moose::stringify(_phase_num), _temperature_variable_name)),
  _viscosity_nodal(declareProperty<Real>("PorousFlow_viscosity" + Moose::stringify(_phase_num))),
  _dviscosity_nodal_dt(declarePropertyDerivative<Real>("PorousFlow_viscosity" + Moose::stringify(_phase_num), _temperature_variable_name)),
  _Mch4(PorousFlowMethaneProperties::molarMass())
{
}

void
PorousFlowMethane::initQpStatefulProperties()
{
  _density_nodal[_qp] = PorousFlowMethaneProperties::density(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num]);
}

void
PorousFlowMethane::computeQpProperties()
{
  /// Density and derivatives wrt pressure and temperature at the nodes
  _density_nodal[_qp] = PorousFlowMethaneProperties::density(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num]);
  _ddensity_nodal_dp[_qp] = PorousFlowMethaneProperties::dDensity_dP(_temperature_nodal[_qp][_phase_num]);
  _ddensity_nodal_dt[_qp] = PorousFlowMethaneProperties::dDensity_dT(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num]);

  /// Density and derivatives wrt pressure and temperature at the qps
  _density_qp[_qp] = PorousFlowMethaneProperties::density(_porepressure_qp[_qp][_phase_num], _temperature_qp[_qp][_phase_num]);
  _ddensity_qp_dp[_qp] = PorousFlowMethaneProperties::dDensity_dP(_temperature_qp[_qp][_phase_num]);
  _ddensity_qp_dt[_qp] = PorousFlowMethaneProperties::dDensity_dT(_porepressure_qp[_qp][_phase_num], _temperature_qp[_qp][_phase_num]);

  /// Viscosity and derivative wrt temperature at the nodes
  _viscosity_nodal[_qp] = PorousFlowMethaneProperties::viscosity(_temperature_nodal[_qp][_phase_num]);
  _dviscosity_nodal_dt[_qp] = PorousFlowMethaneProperties::dViscosity_dT(_temperature_nodal[_qp][_phase_num]);
}
