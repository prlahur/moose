/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef CHEMICALREACTIONSPOROSITY_H
#define CHEMICALREACTIONSPOROSITY_H

#include "Material.h"
#include "DerivativeMaterialInterface.h"

class ChemicalReactionsPorosity;

template <>
InputParameters validParams<ChemicalReactionsPorosity>();

/**
 * Test Material to calculate porosity as a function of coupled
 * solid kinetic species
 */
class ChemicalReactionsPorosity : public DerivativeMaterialInterface<Material>
{
public:
  ChemicalReactionsPorosity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Vector of coupled mineral species
  std::vector<const VariableValue *> _coupled_species;
  /// Number of coupled mineral species
  const unsigned int _num_coupled_species;
  /// Coupled mineral species names
  std::vector<VariableName> _coupled_species_names;
  /// Molar volume of coupled mineral species
  const std::vector<Real> _molar_volume;
  /// Initial porosity
  const Real _initial_porosity;
  /// Porosity
  MaterialProperty<Real> & _porosity;
  /// Derivative of porosity wrt coupled species concentrations
  std::vector<MaterialProperty<Real> *> _dporosity_dspecies;
};

#endif // CHEMICALREACTIONSPOROSITY_H
