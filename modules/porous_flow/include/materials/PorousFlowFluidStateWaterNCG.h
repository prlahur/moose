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
   * @param temperature fluid temperature (C)
   * @return xncgl mass fraction of NCG in water (-)
   */
  Real dissolved(Real pressure, Real temperature) const;


  /// Pore pressure at the nodes
  const MaterialProperty<std::vector<Real> > & _porepressure_nodal;
  /// Pore pressure at the qps
  const MaterialProperty<std::vector<Real> > & _porepressure_qp;
  /// Fluid temperature at the nodes
  const MaterialProperty<Real> & _temperature_nodal;
  /// Fluid temperature at the qps
  const MaterialProperty<Real> & _temperature_qp;
  /// Mass fraction matrix
  MaterialProperty<std::vector<std::vector<Real> > > & _mass_frac;
  /// Old value of mass fraction matrix
  MaterialProperty<std::vector<std::vector<Real> > > & _mass_frac_old;
  /// Gradient of the mass fraction matrix at the quad points
  MaterialProperty<std::vector<std::vector<RealGradient> > > & _grad_mass_frac;
  /// Derivative of the mass fraction matrix with respect to the Porous Flow variables
  MaterialProperty<std::vector<std::vector<std::vector<Real> > > > & _dmass_frac_dvar;
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
  /// Conversion from degrees Celsius to degrees Kelvin
  const Real _T_c2k;
};

#endif //POROUSFLOWFLUIDSTATEWATERNCG_H
