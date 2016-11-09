/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWCHEMISTRYCONVECTION_H
#define POROUSFLOWCHEMISTRYCONVECTION_H

#include "PorousFlowDictator.h"
#include "PorousFlowDarcyBase.h"

// Forward Declarations
class PorousFlowChemistryConvection;

template<>
InputParameters validParams<PorousFlowChemistryConvection>();

/**
 * Kernel implements the convection of the generalised concentrations used
in the geochemistry, summed over phases. 
 */
class PorousFlowChemistryConvection : public PorousFlowDarcyBase
{
public:
  PorousFlowChemistryConvection(const InputParameters & parameters);

protected:


  virtual Real mobility(unsigned nodenum, unsigned phase) override;


  virtual Real dmobility(unsigned nodenum, unsigned phase, unsigned pvar) override;
  

  /// Generalised concentrations in each phase
  const MaterialProperty<std::vector<std::vector<Real> > > & _gen_concentration;

  /// Derivative of generalised concentrations in each phase
  const MaterialProperty<std::vector<std::vector<std::vector<Real> > > > & _dgen_concentration_dvar;  
  
  /// Relative permeability of each phase
  const MaterialProperty<std::vector<Real> > & _relative_permeability;

  /// Derivative of relative permeability of each phase
  const MaterialProperty<std::vector<std::vector<Real> > > & _drelative_permeability_dvar;
  
  
  /// Number of concentration variables supplied (should be == _num_phases)
  const unsigned _primary_species;




};

#endif //POROUSFLOWCHEMISTRYCONVECTION_H
