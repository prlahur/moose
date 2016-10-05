/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWFLUIDSTATEWATERNCG_H
#define POROUSFLOWFLUIDSTATEWATERNCG_H

#include "PorousFlowFluidStateBase.h"
#include "Water97FluidProperties.h"

class PorousFlowFluidStateWaterNCG;

template<>
InputParameters validParams<PorousFlowFluidStateWaterNCG>();

/**
 * Fluid state class for water and a non-condensable gas. Calculates the solubility
 * of the gas phase in the water using Henry's law, and provides density, viscosity
 * and mass fractions for use in Kernels.
 */
class PorousFlowFluidStateWaterNCG : public PorousFlowFluidStateBase
{
public:
  PorousFlowFluidStateWaterNCG(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeQpProperties() override;

  virtual void thermophysicalProperties(Ecalc node_or_qp, Eprops props) const override;

  /**
   * Dissolved mass fraction of NCG in water calculated using Henry's law
   *
   * @param pressure NCG partial pressure (Pa)
   * @param temperature fluid temperature (K)
   * @param[out] xncgl mass fraction of NCG in water (-)
   * @param[out] dxncgl_dp derivative of mass fraction of NCG in water wrt pressure
   * @param[out] dxncgl_dT derivative of mass fraction of NCG in water wrt temperature
   */
  void dissolved(Real pressure, Real temperature, Real & xncgl, Real & dxncgl_dp, Real & dxncgl_dT) const;

  /**
   * Enthalpy of dissolution of NCG in water calculated using Henry's constant
   * From Himmelblau, Partial molal heats and entropies of solution for gases dissolved
   * in water from the freezing to the near critical point, J. Phys. Chem. 63 (1959)
   *
   * @param temperature fluid temperature (K)
   * @param Kh Henry's constant (Pa)
   * @param dKh_dT derivative of Henry's constant wrt temperature
   * @return enthalpy of dissolution (kJ/kg)
   */
  Real enthalpyDissolution(Real temperature, Real Kh, Real dKh_dT) const;

  /// Phase number of the aqueous phase
  const unsigned int _aqueous_phase_number;
  /// Phase number of the NCG phase
  const unsigned int _ncg_phase_number;
  /// Fluid component number of the aqueous component
  const unsigned int _aqueous_fluid_component;
  /// Fluid component number of the NCG phase
  const unsigned int _ncg_fluid_component;
  /// Fluid properties UserObject for water
  const Water97FluidProperties & _water_fp;
  /// Fluid properties UserObject for the NCG
  const SinglePhaseFluidPropertiesPT & _ncg_fp;
  /// Molar mass of water (kg/mol)
  const Real _Mh2o;
  /// Molar mass of non-condensable gas (kg/mol)
  const Real _Mncg;
};

#endif //POROUSFLOWFLUIDSTATEWATERNCG_H
