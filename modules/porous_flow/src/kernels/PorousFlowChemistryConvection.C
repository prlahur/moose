/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowChemistryConvection.h"

template<>
InputParameters validParams<PorousFlowChemistryConvection>()
{
  InputParameters params = validParams<PorousFlowDarcyBase>();  
  params.addClassDescription("Convection of generalised concentrations in chemical reactions, with upwinding");
  params.addParam<unsigned int>("primary_species", 0, "The index corresponding to the primary species for this kernel");   
  return params;
}

PorousFlowChemistryConvection::PorousFlowChemistryConvection(const InputParameters & parameters) :
    PorousFlowDarcyBase(parameters),
    _gen_concentration(getMaterialProperty<std::vector<std::vector<Real> > > ("PorousFlow_gen_conc")),
    _dgen_concentration_dvar(getMaterialProperty<std::vector<std::vector<std::vector<Real> > > > ("dPorousFlow_gen_conc_dvar")),    
    _relative_permeability(getMaterialProperty<std::vector<Real> > ("PorousFlow_relative_permeability")),
    _drelative_permeability_dvar(getMaterialProperty<std::vector<std::vector<Real> > > ("dPorousFlow_relative_permeability_dvar")),    
    _primary_species(getParam<unsigned int> ("primary_species"))

{

}

Real PorousFlowChemistryConvection::mobility(unsigned nodenum, unsigned phase)
{

  return _gen_concentration[nodenum][phase][_primary_species] * _relative_permeability[nodenum][phase] / _fluid_viscosity[nodenum][phase];
}

Real PorousFlowChemistryConvection::dmobility(unsigned nodenum, unsigned phase, unsigned pvar)
{
  Real dm = _dgen_concentration_dvar[nodenum][phase][_primary_species][pvar] * _relative_permeability[nodenum][phase] / _fluid_viscosity[nodenum][phase];
  dm += _gen_concentration[nodenum][phase][_primary_species] * _drelative_permeability_dvar[nodenum][phase][pvar] / _fluid_viscosity[nodenum][phase];
  dm -= _gen_concentration[nodenum][phase][_primary_species] * _relative_permeability[nodenum][phase] * _dfluid_viscosity_dvar[nodenum][phase][pvar] / std::pow(_fluid_viscosity[nodenum][phase], 2);
  return dm;
}
