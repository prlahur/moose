# Example of accessing properties using the PorousFlowPropertyAux AuxKernel for
# each phase and fluid component (as required).

[Mesh]
  type = GeneratedMesh
  dim = 2
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [./pwater]
    initial_condition = 1e6
  [../]
  [./sgas]
    initial_condition = 0.3
  [../]
[]

[AuxVariables]
  [./pressure_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./saturation_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./density_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./density_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./viscosity_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./viscosity_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x0_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x0_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x1_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x1_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./relperm_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./relperm_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  #[./enthalpy_water]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./enthalpy_gas]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./energy_water]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./energy_gas]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
[]

[AuxKernels]
  [./pressure_gas]
    type = PorousFlowPropertyAux
    variable = pressure_gas
    property = pressure
    phase = 1
    execute_on = timestep_end
  [../]
  [./saturation_water]
    type = PorousFlowPropertyAux
    variable = saturation_water
    property = saturation
    phase = 0
    execute_on = timestep_end
  [../]
  [./density_water]
    type = PorousFlowPropertyAux
    variable = density_water
    property = density
    phase = 0
    execute_on = timestep_end
  [../]
  [./density_gas]
    type = PorousFlowPropertyAux
    variable = density_gas
    property = density
    phase = 1
    execute_on = timestep_end
  [../]
  [./viscosity_water]
    type = PorousFlowPropertyAux
    variable = viscosity_water
    property = viscosity
    phase = 0
    execute_on = timestep_end
  [../]
  [./viscosity_gas]
    type = PorousFlowPropertyAux
    variable = viscosity_gas
    property = viscosity
    phase = 1
    execute_on = timestep_end
  [../]
  [./relperm_water]
    type = PorousFlowPropertyAux
    variable = relperm_water
    property = relperm
    phase = 0
    execute_on = timestep_end
  [../]
  [./relperm_gas]
    type = PorousFlowPropertyAux
    variable = relperm_gas
    property = relperm
    phase = 1
    execute_on = timestep_end
  [../]
  [./x1_water]
    type = PorousFlowPropertyAux
    variable = x1_water
    property = mass_fraction
    phase = 0
    fluid_component = 1
    execute_on = timestep_end
  [../]
  [./x1_gas]
    type = PorousFlowPropertyAux
    variable = x1_gas
    property = mass_fraction
    phase = 1
    fluid_component = 1
    execute_on = timestep_end
  [../]
  [./x0_water]
    type = PorousFlowPropertyAux
    variable = x0_water
    property = mass_fraction
    phase = 0
    fluid_component = 0
    execute_on = timestep_end
  [../]
  [./x0_gas]
    type = PorousFlowPropertyAux
    variable = x0_gas
    property = mass_fraction
    phase = 1
    fluid_component = 0
    execute_on = timestep_end
  [../]
  #[./enthalpy_water]
  #  type = PorousFlowPropertyAux
  #  variable = enthalpy_water
  #  property = enthalpy
  #  phase = 0
  #  execute_on = timestep_end
  #[../]
  #[./enthalpy_gas]
  #  type = PorousFlowPropertyAux
  #  variable = enthalpy_gas
  #  property = enthalpy
  #  phase = 1
  #  execute_on = timestep_end
  #[../]
  #[./energy_water]
  #  type = PorousFlowPropertyAux
  #  variable = energy_water
  #  property = internal_energy
  #  phase = 0
  #  execute_on = timestep_end
  #[../]
  #[./energy_gas]
  #  type = PorousFlowPropertyAux
  #  variable = energy_gas
  #  property = internal_energy
  #  phase = 1
  #  execute_on = timestep_end
  #[../]
[]

[Kernels]
  [./mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pwater
  [../]
  [./flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pwater
  [../]
  [./mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = sgas
  [../]
  [./flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = sgas
  [../]
  #[./energy_dot]
  #  type = PorousFlowEnergyTimeDerivative
  #  variable = temperature
  #[../]
  #[./convection]
  #  type = PorousFlowConvectiveFlux
  #  variable = temperature
  #[../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater sgas'
    number_fluid_phases = 2
    number_fluid_components = 2
  [../]
[]

[Modules]
  [./FluidProperties]
    [./co2]
      type = MethaneFluidProperties
    [../]
    [./water]
      type = Water97FluidProperties
    [../]
  [../]
[]

[Materials]
  [./temperature]
    type = PorousFlowTemperature
    temperature = 50
  [../]
  [./nnn]
    type = PorousFlowNodeNumber
    on_initial_only = true
  [../]
  [./ppss]
    type = PorousFlow2PhasePS_VG
    phase0_porepressure = pwater
    phase1_saturation = sgas
    m = 0.5
    p0 = 1e5
    pc_max = 1e7
    sat_lr = 0.1
    sat_ls = 1
  [../]
  [./waterncg]
    type = PorousFlowFluidStateWaterNCG
    gas_fp = co2
    water_fp = water
  [../]
  #[./massfrac]
  #  type = PorousFlowMassFraction
  #  mass_fraction_vars = 'x0_water x0_gas'
  #[../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0 0 1e-12 0 0 0 1e-12'
  [../]
  [./relperm0]
    type = PorousFlowRelativePermeabilityCorey
    n_j = 2
    phase = 0
  [../]
  [./relperm1]
    type = PorousFlowRelativePermeabilityCorey
    n_j = 3
    phase = 1
  [../]
  [./relperm_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_relative_permeability
  [../]
  [./porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  [../]
  [./rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 1.0
    density = 125
  [../]
  #[./fluid_energy0]
  #  type = PorousFlowInternalEnergyIdeal
  #  specific_heat_capacity = 2
  #  phase = 0
  #[../]
  #[./fluid_energy1]
  #  type = PorousFlowInternalEnergyIdeal
  #  specific_heat_capacity = 1
  #  phase = 1
  #[../]
  #[./energy_all]
  #  type = PorousFlowJoiner
  #  include_old = true
  #  material_property = PorousFlow_fluid_phase_internal_energy_nodal
  #[../]
  #[./enthalpy0]
  #  type = PorousFlowEnthalpy
  #  phase = 0
  #[../]
  #[./enthalpy1]
  #  type = PorousFlowEnthalpy
  #  phase = 1
  #[../]
  #[./enthalpy_all]
  #  type = PorousFlowJoiner
  #  material_property = PorousFlow_fluid_phase_enthalpy_nodal
  #[../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 1
  end_time = 1
  nl_abs_tol = 1e-12
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Outputs]
  exodus = true
[]

[Problem]
  solve = false
[]
