/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowFluidStateBase.h"

template<>
InputParameters validParams<PorousFlowFluidStateBase>()
{
  InputParameters params = validParams<PorousFlowMaterialVectorBase>();
  params.addClassDescription("Base class for fluid state calculations");
  return params;
}

PorousFlowFluidStateBase::PorousFlowFluidStateBase(const InputParameters & parameters) :
    PorousFlowMaterialVectorBase(parameters),

    _dporepressure_nodal_dvar(getMaterialProperty<std::vector<std::vector<Real> > >("dPorousFlow_porepressure_nodal_dvar")),
    _dporepressure_qp_dvar(getMaterialProperty<std::vector<std::vector<Real> > >("dPorousFlow_porepressure_qp_dvar")),
    _dsaturation_nodal_dvar(getMaterialProperty<std::vector<std::vector<Real> > >("dPorousFlow_saturation_nodal_dvar")),
    _dsaturation_qp_dvar(getMaterialProperty<std::vector<std::vector<Real> > >("dPorousFlow_saturation_qp_dvar")),
    _dtemperature_nodal_dvar(getMaterialProperty<std::vector<Real> >("dPorousFlow_temperature_nodal_dvar")),
    _dtemperature_qp_dvar(getMaterialProperty<std::vector<Real> >("dPorousFlow_temperature_qp_dvar")),

    _fluid_density_nodal(declareProperty<std::vector<Real> >("PorousFlow_fluid_phase_density")),
    _fluid_density_nodal_old(declarePropertyOld<std::vector<Real> >("PorousFlow_fluid_phase_density")),
    _dfluid_density_nodal_dvar(declareProperty<std::vector<std::vector<Real> > >("dPorousFlow_fluid_phase_density_dvar")),
    _fluid_density_qp(declareProperty<std::vector<Real> >("PorousFlow_fluid_phase_density_qp")),
    _dfluid_density_qp_dvar(declareProperty<std::vector<std::vector<Real> > >("dPorousFlow_fluid_phase_density_qp_dvar")),
    _fluid_viscosity_nodal(declareProperty<std::vector<Real> >("PorousFlow_viscosity")),
    _dfluid_viscosity_nodal_dvar(declareProperty<std::vector<std::vector<Real> > >("dPorousFlow_viscosity_dvar")),

    _pressure_variable_name(_dictator.pressureVariableNameDummy()),
    _saturation_variable_name(_dictator.saturationVariableNameDummy()),
    _temperature_variable_name(_dictator.temperatureVariableNameDummy()),
    _mass_fraction_variable_name(_dictator.massFractionVariableNameDummy())
{
}

void
PorousFlowFluidStateBase::initQpStatefulProperties()
{
  mooseError("initQpStatefulProperties() must be overwritten by any object inheriting from PorousFlowFluidStateBase");
}

void
PorousFlowFluidStateBase::computeQpProperties()
{
  mooseError("computeQpProperties() must be overwritten by any object inheriting from PorousFlowFluidStateBase");
}
