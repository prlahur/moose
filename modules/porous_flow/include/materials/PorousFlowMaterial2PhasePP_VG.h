/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#ifndef PORFLOWMATERIAL2PHASEPP_VG_H
#define PORFLOWMATERIAL2PHASEPP_VG_H

#include "PorousFlowVariableBase.h"
#include "PorousFlowCapillaryVG.h"

//Forward Declarations
class PorousFlowMaterial2PhasePP_VG;

template<>
InputParameters validParams<PorousFlowMaterial2PhasePP_VG>();

/**
 * Material designed to calculate fluid-phase porepressures and saturations at nodes and quadpoints
 */
class PorousFlowMaterial2PhasePP_VG : public PorousFlowVariableBase
{
public:
  PorousFlowMaterial2PhasePP_VG(const InputParameters & parameters);

protected:

  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  /// vanGenuchten alpha
  const Real _al;
  /// vanGenuchten m
  const Real _m;
  /// Nodal value of porepressure of the zero phase (eg, the water phase)
  const VariableValue & _phase0_porepressure_nodal;
  /// Quadpoint value of porepressure of the zero phase (eg, the water phase)
  const VariableValue & _phase0_porepressure_qp;
  /// Gradient(phase0_porepressure) at the qps
  const VariableGradient & _phase0_gradp_qp;
  /// Moose variable number of the phase0 porepressure
  const unsigned int _phase0_porepressure_varnum;
  /// Nodal value of porepressure of the one phase (eg, the gas phase)
  const VariableValue & _phase1_porepressure_nodal;
  /// Quadpoint value of porepressure of the one phase (eg, the gas phase)
  const VariableValue & _phase1_porepressure_qp;
  /// Gradient(phase1_porepressure) at the qps
  const VariableGradient & _phase1_gradp_qp;
  /// Moose variable number of the phase1 porepressure
  const unsigned int _phase1_porepressure_varnum;

 private:
  void buildQpPPSS();
};

#endif //PORFLOWMATERIAL2PHASEPP_VG_H
