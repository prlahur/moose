/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ChemicalReactionsPorosity.h"

template <>
InputParameters
validParams<ChemicalReactionsPorosity>()
{
  InputParameters params = validParams<Material>();
  params.addCoupledVar("coupled_species", "The coupled mineral species");
  params.addRequiredParam<std::vector<Real>>("molar_volume",
                                             "Molar volume of coupled mineral species");
  params.addRequiredParam<Real>("intial_porosity", "Initial porosity");
  params.addClassDescription("Porosity as a function of coupled concentrations");
  return params;
}

ChemicalReactionsPorosity::ChemicalReactionsPorosity(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _num_coupled_species(coupledComponents("coupled_species")),
    _molar_volume(getParam<std::vector<Real>>("molar_volume")),
    _initial_porosity(getParam<Real>("intial_porosity")),
    _porosity(declareProperty<Real>("porosity"))
{
  _coupled_species.resize(_num_coupled_species);
  _coupled_species_names.resize(_num_coupled_species);

  for (unsigned int i = 0; i < _num_coupled_species; ++i)
  {
    _coupled_species[i] = &coupledValue("coupled_species", i);
    _coupled_species_names[i] = getVar("coupled_species", i)->name();
  }
}

void
ChemicalReactionsPorosity::computeQpProperties()
{
  Real porosity = _initial_porosity;

  for (unsigned int i = 0; i < _num_coupled_species; ++i)
    porosity += (*_coupled_species[i])[_qp];

  _porosity[_qp] = porosity;
}
