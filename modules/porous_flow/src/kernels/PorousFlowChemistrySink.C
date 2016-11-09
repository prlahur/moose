/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowChemistrySink.h"

template<>
InputParameters validParams<PorousFlowChemistrySink>()
{
  InputParameters params = validParams<Kernel>();  
  params.addClassDescription("Sink of primary concentrations in chemical reactions.  THIS IS CURRENTLY VERY SPECIALISED WITH NO ERROR CHECKING!!  Jonathons eqn229");
  params.addRequiredParam<UserObjectName>("PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addParam<unsigned int>("primary_species", 0, "The index corresponding to the primary species for this kernel");   
  params.addParam<unsigned int>("phase", 0, "The index corresponding to the primary species for this kernel");   
  params.addRequiredParam<Real>("coefficient", "The coefficient in front of C - Ceq");
  params.addParam<Real>("Ceq", 0, "Ceq");
  return params;
}

PorousFlowChemistrySink::PorousFlowChemistrySink(const InputParameters & parameters) :
    Kernel(parameters),
    _gen_concentration(getMaterialProperty<std::vector<std::vector<Real> > > ("PorousFlow_gen_conc")),
    _dgen_concentration_dvar(getMaterialProperty<std::vector<std::vector<std::vector<Real> > > > ("dPorousFlow_gen_conc_dvar")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _phase(getParam<unsigned int>("phase")),
    _primary_species(getParam<unsigned int>("primary_species")),
    _coeff(getParam<Real>("coefficient")),
    _ceq(getParam<Real>("Ceq"))
{
}

Real
PorousFlowChemistrySink::computeQpResidual()
{
  // lumped to the nodes
  return _test[_i][_qp] * _coeff * (_gen_concentration[_i][_phase][_primary_species] - _ceq);
}

Real
PorousFlowChemistrySink::computeQpJacobian()
{
  return computeQpOffDiagJacobian(_var.number());
}

Real
PorousFlowChemistrySink::computeQpOffDiagJacobian(unsigned jvar)
{
  if (_i != _j)
    return 0.0; // because of lumping
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  const unsigned pvar = _dictator.porousFlowVariableNum(jvar);
  return _test[_i][_qp] * _coeff * _dgen_concentration_dvar[_qp][_phase][_primary_species][pvar]; // no _phi[_j][_qp] because of lumping;
}
