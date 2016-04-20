/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef POROUSFLOWVOLUMETRICSTRAIN_H
#define POROUSFLOWVOLUMETRICSTRAIN_H

#include "DerivativeMaterialInterface.h"
#include "RankTwoTensor.h"
#include "Material.h"

#include "PorousFlowDictator.h"

/**
 * PorousFlowVolumetricStrain computes volumetric strains, and derivatives thereof
 */
class PorousFlowVolumetricStrain : public DerivativeMaterialInterface<Material>
{
public:
  PorousFlowVolumetricStrain(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  /// The dictator UserObject for the Porous-Flow simulation
  const PorousFlowDictator & _dictator_UO;

  /// number of porous-flow variables
  const unsigned int _num_var;

  // number of displacements supplied
  const unsigned int _ndisp;

  // displacement variables and their derivatives
  std::vector<const VariableValue *> _disp;
  /// variable number of the displacements variables
  std::vector<unsigned int> _disp_var_num;
  std::vector<const VariableGradient *> _grad_disp;
  std::vector<const VariableGradient *> _grad_disp_old;

  MaterialProperty<Real> & _vol_strain_rate_qp;
  MaterialProperty<std::vector<RealGradient> > & _dvol_strain_rate_qp_dvar;
  MaterialProperty<Real> & _vol_total_strain_qp;
  MaterialProperty<std::vector<RealGradient> > & _dvol_total_strain_qp_dvar;
};

#endif //POROUSFLOWVOLUMETRICSTRAIN_H
