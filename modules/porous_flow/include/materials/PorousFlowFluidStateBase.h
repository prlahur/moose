/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWFLUIDSTATEBASE_H
#define POROUSFLOWFLUIDSTATEBASE_H

#include "PorousFlowMaterialVectorBase.h"

class PorousFlowFluidStateBase;

template<>
InputParameters validParams<PorousFlowFluidStateBase>();

/**
 * Base class for fluid states
 */
class PorousFlowFluidStateBase : public PorousFlowMaterialVectorBase
{
public:
  PorousFlowFluidStateBase(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Enum to control calculation of properties at either the nodes or qps
  enum Ecalc { AT_NODE = 0, AT_QP = 1 };
  /// Enum to control calculation of only stateful properties or all properties
  enum Eprops { ALL_PROPERTIES = 0, STATEFUL_ONLY = 1};

  /**
   * Calculates the thermophysical properties and derivatives for each phase
   * and fluid component at the nodes and qps as required. If stateful is STATEFUL,
   * only those properties where old values are required are calculated
   */
  virtual void thermophysicalProperties(Ecalc node_or_qp, Eprops props) const = 0;

  /// Derivatives of porepressure variable wrt PorousFlow variables at the nodes
  const MaterialProperty<std::vector<std::vector<Real> > > & _dporepressure_nodal_dvar;

  /// Derivatives of porepressure variable wrt PorousFlow variables at the qps
  const MaterialProperty<std::vector<std::vector<Real> > > & _dporepressure_qp_dvar;

  /// Derivatives of saturation variable wrt PorousFlow variables at the nodes
  const MaterialProperty<std::vector<std::vector<Real> > > & _dsaturation_nodal_dvar;

  /// Derivatives of saturation variable wrt PorousFlow variables at the qps
  const MaterialProperty<std::vector<std::vector<Real> > > & _dsaturation_qp_dvar;

  /// Derivatives of temperature variable wrt PorousFlow variables at the nodes
  const MaterialProperty<std::vector<Real> > & _dtemperature_nodal_dvar;

  /// Derivatives of temperature variable wrt PorousFlow variables at the qps
  const MaterialProperty<std::vector<Real> > & _dtemperature_qp_dvar;

  /// Fluid density for each phase at the nodes
  MaterialProperty<std::vector<Real> > & _fluid_density_nodal;

  /// Old value of fluid density for each phase at the nodes
  MaterialProperty<std::vector<Real> > & _fluid_density_nodal_old;

  /// Derivative of the fluid density for each phase wrt PorousFlow variables at the nodes
  MaterialProperty<std::vector<std::vector<Real> > > & _dfluid_density_nodal_dvar;

  /// Fluid density for each phase at the qps
  MaterialProperty<std::vector<Real> > & _fluid_density_qp;

  /// Derivative of the fluid density for each phase wrt PorousFlow variables at the qps
  MaterialProperty<std::vector<std::vector<Real> > > & _dfluid_density_qp_dvar;

  /// Viscosity of each component in each phase at the nodes
  MaterialProperty<std::vector<Real> > & _fluid_viscosity_nodal;

  /// Derivative of the fluid viscosity for each phase wrt PorousFlow variables at the nodes
  MaterialProperty<std::vector<std::vector<Real> > > & _dfluid_viscosity_nodal_dvar;

  /// Name of (dummy) pressure variable
  const VariableName _pressure_variable_name;

  /// Name of (dummy) saturation variable
  const VariableName _saturation_variable_name;

  /// Name of (dummy) temperature variable
  const VariableName _temperature_variable_name;

  /// Name of (dummy) mass fraction variable
  const VariableName _mass_fraction_variable_name;
};

#endif //POROUSFLOWFLUIDSTATEBASE_H
