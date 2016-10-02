/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowFluidStateWaterNCG.h"

template<>
InputParameters validParams<PorousFlowFluidStateWaterNCG>()
{
  InputParameters params = validParams<PorousFlowFluidStateBase>();
  params.addParam<unsigned int>("water_phase_number", 0, "The phase number of the aqueous phase");
  params.addParam<unsigned int>("water_fluid_component", 0, "The fluid component number of the aqueous phase");
  params.addRequiredParam<UserObjectName>("water_fp", "The name of the user object for water");
  params.addRequiredParam<UserObjectName>("gas_fp", "The name of the user object for the non-condensable gas");
  params.addClassDescription("Fluid state class for water and non-condensable gas");
  return params;
}

PorousFlowFluidStateWaterNCG::PorousFlowFluidStateWaterNCG(const InputParameters & parameters) :
    PorousFlowFluidStateBase(parameters),

    _porepressure_nodal(getMaterialProperty<std::vector<Real> >("PorousFlow_porepressure_nodal")),
    _porepressure_qp(getMaterialProperty<std::vector<Real> >("PorousFlow_porepressure_qp")),
    _temperature_nodal(getMaterialProperty<Real>("PorousFlow_temperature_nodal")),
    _temperature_qp(getMaterialProperty<Real>("PorousFlow_temperature_qp")),
    _mass_frac(declareProperty<std::vector<std::vector<Real> > >("PorousFlow_mass_frac")),
    _mass_frac_old(declarePropertyOld<std::vector<std::vector<Real> > >("PorousFlow_mass_frac")),
    _grad_mass_frac(declareProperty<std::vector<std::vector<RealGradient> > >("PorousFlow_grad_mass_frac")),
    _dmass_frac_dvar(declareProperty<std::vector<std::vector<std::vector<Real> > > >("dPorousFlow_mass_frac_dvar")),
    _aqueous_phase_number(getParam<unsigned int>("water_phase_number")),
    _ncg_phase_number(1 - _aqueous_phase_number),
    _aqueous_fluid_component(getParam<unsigned int>("water_fluid_component")),
    _ncg_fluid_component(1 - _aqueous_fluid_component),
    _water_fp(getUserObject<Water97FluidProperties>("water_fp")),
    _ncg_fp(getUserObject<SinglePhaseFluidPropertiesPT>("gas_fp")),
    _Mh2o(_water_fp.molarMass()),
//  _Mncg(_ncg_fp.molarMass())
    _Mncg(44.0098e-3),//FIXME: delete when molarmass() is available
    _T_c2k(273.15)
{
}

void
PorousFlowFluidStateWaterNCG::initQpStatefulProperties()
{
  // Set the size of all the vectors in the FluidStateProperties object
  _fluid_density_nodal[_qp].resize(_num_phases);
  _fluid_density_qp[_qp].resize(_num_phases);
  _fluid_viscosity_nodal[_qp].resize(_num_phases);

  _mass_frac[_qp].resize(_num_phases);
  for (unsigned int ph = 0; ph < _num_phases; ++ph)
    _mass_frac[_qp][ph].resize(_num_components);

  // Set the size of the derivatives wrt PorousFlow variables
  _dfluid_density_nodal_dvar[_qp].resize(_num_phases);
  _dfluid_density_qp_dvar[_qp].resize(_num_phases);
  _dfluid_viscosity_nodal_dvar[_qp].resize(_num_phases);

  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    _dfluid_density_nodal_dvar[_qp][ph].resize(_num_var);
    _dfluid_density_qp_dvar[_qp][ph].resize(_num_var);
    _dfluid_viscosity_nodal_dvar[_qp][ph].resize(_num_var);
  }

  thermophysicalProperties(AT_NODE, STATEFUL_ONLY);
}

void
PorousFlowFluidStateWaterNCG::computeQpProperties()
{
  // Calculate all properties at the nodes and qps as required
  thermophysicalProperties(AT_NODE, ALL_PROPERTIES);
  thermophysicalProperties(AT_QP, ALL_PROPERTIES);
}

void
PorousFlowFluidStateWaterNCG::thermophysicalProperties(Ecalc node_or_qp, Eprops props) const
{
  // Select either nodal or qp input variables
  std::vector<Real> pressure = (node_or_qp ? _porepressure_qp[_qp] : _porepressure_nodal[_qp]);
  Real temperature = (node_or_qp ? _temperature_qp[_qp] : _temperature_nodal[_qp]);

  // The FluidProperty objects use temperature in K
  Real Tk = temperature + _T_c2k;

  // Vapor pressure
  Real pv = _water_fp.pSat(Tk);

  // Partial pressure of NCG
  Real pncg = pressure[_ncg_phase_number] - pv;

  // The density of the liquid water phase
  Real water_density = _water_fp.rho(pressure[_aqueous_phase_number], Tk);

  // The density of the water vapor phase
  Real vapor_density = _water_fp.rho(pv, Tk);

  // The density of the NCG in the gas phase
  Real ncg_density = _ncg_fp.rho(pncg, Tk);

  // The gas phase density is the sum of densities
  Real gas_density = ncg_density + vapor_density;

  // The mass fraction of NCG in liquid and gas phases
  Real xncgl = dissolved(pncg, Tk);
  Real xncgg = ncg_density / gas_density;

  // Assume that the liquid density is not modified by the solute
  // TODO: modify this to account for changes due to solute
  Real liquid_density = water_density;

  // The properties in each phase can now be stored in the appropriate std::vector
  if (node_or_qp == AT_QP)
  {
    _fluid_density_qp[_qp][_aqueous_phase_number] = liquid_density;
    _fluid_density_qp[_qp][_ncg_phase_number] = gas_density;
  }
  else
  {
    _fluid_density_nodal[_qp][_aqueous_phase_number] = liquid_density;
    _fluid_density_nodal[_qp][_ncg_phase_number] = gas_density;

    _mass_frac[_qp][_aqueous_phase_number][_aqueous_fluid_component] = 1.0 - xncgl;
    _mass_frac[_qp][_aqueous_phase_number][_ncg_fluid_component] = xncgl;
    _mass_frac[_qp][_ncg_phase_number][_aqueous_fluid_component] = 1.0 - xncgg;
    _mass_frac[_qp][_ncg_phase_number][_ncg_fluid_component] = xncgg;

    // Also calculate properties at the node that aren't stateful
    if (props == ALL_PROPERTIES)
    {
      // The viscosity of the liquid water phase
      Real water_viscosity = _water_fp.mu(water_density, Tk);
      // The viscosity of the water vapor phase
      Real vapor_viscosity = _water_fp.mu(vapor_density, Tk);
      // The viscosity of the NCG in the gas phase
      Real ncg_viscosity = _ncg_fp.mu(ncg_density, Tk);
      // Assume that the viscosities are not modified by the solute
      Real gas_viscosity = ncg_viscosity;
      Real liquid_viscosity = water_viscosity;

      _fluid_viscosity_nodal[_qp][_aqueous_phase_number] = liquid_viscosity;
      _fluid_viscosity_nodal[_qp][_ncg_phase_number] = gas_viscosity;
    }
  }
}

Real
PorousFlowFluidStateWaterNCG::dissolved(Real pressure, Real temperature) const
{
  // The dissolved mole fraction of NCG in water is given by Henry's law
  Real Xncgl = pressure / _ncg_fp.henryConstant(temperature);
  // The mass fraction is then
  Real xcncgl = Xncgl * _Mncg / (Xncgl * _Mncg + (1.0 - Xncgl) * _Mh2o);

  return xcncgl;
}
