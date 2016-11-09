/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowChemistryTimeDerivative.h"

template<>
InputParameters validParams<PorousFlowChemistryTimeDerivative>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addParam<unsigned int>("primary_species", 0, "The index corresponding to the primary species for this kernel");
  params.addRequiredParam<UserObjectName>("PorousFlowDictator", "The UserObject that holds the list of Porous-Flow variable names.");
  params.addClassDescription("Concentration derivative wrt time for primary species given by primary_species");
  return params;
}

PorousFlowChemistryTimeDerivative::PorousFlowChemistryTimeDerivative(const InputParameters & parameters) :
    TimeKernel(parameters),
    _primary_species(getParam<unsigned int>("primary_species")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
    _num_phases(_dictator.numPhases()),
    _porosity(getMaterialProperty<Real>("PorousFlow_porosity_nodal")),
    _porosity_old(getMaterialPropertyOld<Real>("PorousFlow_porosity_nodal")),
    _dporosity_dvar(getMaterialProperty<std::vector<Real> >("dPorousFlow_porosity_nodal_dvar")),
    _dporosity_dgradvar(getMaterialProperty<std::vector<RealGradient> >("dPorousFlow_porosity_nodal_dgradvar")),
    _fluid_saturation_nodal(getMaterialProperty<std::vector<Real> >("PorousFlow_saturation_nodal")),
    _fluid_saturation_nodal_old(getMaterialPropertyOld<std::vector<Real> >("PorousFlow_saturation_nodal")),
    _dfluid_saturation_nodal_dvar(getMaterialProperty<std::vector<std::vector<Real> > >("dPorousFlow_saturation_nodal_dvar")),
    _gen_concentration(getMaterialProperty<std::vector<std::vector<Real> > >("PorousFlow_gen_conc")),
    _gen_concentration_old(getMaterialPropertyOld<std::vector<std::vector<Real> > >("PorousFlow_gen_conc")),
    _dgen_concentration_dvar(getMaterialProperty<std::vector<std::vector<std::vector<Real> > > >("dPorousFlow_gen_conc_dvar"))
{
  if (_primary_species >= _dictator.numPrimarySpecies())
    mooseError("The Dictator proclaims that the number of primary species in this simulation is " << _dictator.numPrimarySpecies() << " whereas you have used the Kernel PorousFlowComponetChemistryTimeDerivative with primary_species = " << _primary_species << ".  The Dictator does not take such mistakes lightly");
}

Real
PorousFlowChemistryTimeDerivative::computeQpResidual()
{
  Real mass = 0.0;
  Real mass_old = 0.0;
  for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
    mass += _fluid_saturation_nodal[_i][ph] * _gen_concentration[_i][ph][_primary_species];
    mass_old += _fluid_saturation_nodal_old[_i][ph] * _gen_concentration_old[_i][ph][_primary_species];
   }

  return _test[_i][_qp] * (_porosity[_i] * mass - _porosity_old[_i] * mass_old) / _dt;
}

Real
PorousFlowChemistryTimeDerivative::computeQpJacobian()
{
  /// If the variable is not a PorousFlow variable (very unusual), the diag Jacobian terms are 0
  if (!_var_is_porflow_var)
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(_var.number()));
}

Real
PorousFlowChemistryTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  /// If the variable is not a PorousFlow variable, the OffDiag Jacobian terms are 0
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(jvar));
}

Real
PorousFlowChemistryTimeDerivative::computeQpJac(unsigned int pvar)
{
  // porosity is dependent on variables that are lumped to the nodes,
  // but it can depend on the gradient
  // of variables, which are NOT lumped to the nodes, hence:
  Real dmass = 0.0;
  for (unsigned ph = 0; ph < _num_phases; ++ph)
    dmass += _fluid_saturation_nodal[_i][ph] * _gen_concentration[_i][ph][_primary_species] * _dporosity_dgradvar[_i][pvar] * _grad_phi[_j][_i];

  if (_i != _j)
    return _test[_i][_qp] * dmass/_dt;

  /// As the fluid mass is lumped to the nodes, only non-zero terms are for _i==_j
  for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
    dmass += _dfluid_saturation_nodal_dvar[_i][ph][pvar] * _gen_concentration[_i][ph][_primary_species] * _porosity[_i];
    dmass += _fluid_saturation_nodal[_i][ph] * _dgen_concentration_dvar[_i][ph][_primary_species][pvar] * _porosity[_i];
    dmass += _fluid_saturation_nodal[_i][ph] * _gen_concentration[_i][ph][_primary_species] * _dporosity_dvar[_i][pvar];
  }
  return _test[_i][_qp] * dmass / _dt;
}
