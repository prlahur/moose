/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowBrine.h"
#include "Conversion.h"

template<>
InputParameters validParams<PorousFlowBrine>()
{
  InputParameters params = validParams<PorousFlowFluidPropertiesBase>();
  params.addRequiredParam<Real>("xnacl", "The salt mass fraction in the brine (kg/kg)");
  params.addClassDescription("This Material calculates fluid properties for brine (H20 and NaCl)");
  return params;
}

PorousFlowBrine::PorousFlowBrine(const InputParameters & parameters) :
    PorousFlowFluidPropertiesBase(parameters),

  _density_nodal(declareProperty<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num))),
  _density_nodal_old(declarePropertyOld<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num))),
  _ddensity_nodal_dp(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density" + Moose::stringify(_phase_num), _pressure_variable_name)),
  _density_qp(declareProperty<Real>("PorousFlow_fluid_phase_density_qp" + Moose::stringify(_phase_num))),
  _ddensity_qp_dp(declarePropertyDerivative<Real>("PorousFlow_fluid_phase_density_qp" + Moose::stringify(_phase_num), _pressure_variable_name)),
  _viscosity_nodal(declareProperty<Real>("PorousFlow_viscosity" + Moose::stringify(_phase_num))),
  _xnacl(getParam<Real>("xnacl"))
{
}

void
PorousFlowBrine::initQpStatefulProperties()
{
  _density_nodal[_qp] = PorousFlowBrineProperties::density(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num], _xnacl);
}

void
PorousFlowBrine::computeQpProperties()
{
  /// Density and derivatives wrt pressure and temperature at the nodes
  _density_nodal[_qp] = PorousFlowBrineProperties::density(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num], _xnacl);
  _ddensity_nodal_dp[_qp] = PorousFlowBrineProperties::dDensity_dP(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num], _xnacl);

  /// Density and derivatives wrt pressure and temperature at the qps
  _density_qp[_qp] = PorousFlowBrineProperties::density(_porepressure_qp[_qp][_phase_num], _temperature_qp[_qp][_phase_num], _xnacl);
  _ddensity_qp_dp[_qp] = PorousFlowBrineProperties::dDensity_dP(_porepressure_qp[_qp][_phase_num], _temperature_qp[_qp][_phase_num], _xnacl);

  /// Viscosity and derivative wrt temperature at the nodes
  _viscosity_nodal[_qp] = PorousFlowBrineProperties::viscosity(_porepressure_nodal[_qp][_phase_num], _temperature_nodal[_qp][_phase_num], _xnacl);
}
