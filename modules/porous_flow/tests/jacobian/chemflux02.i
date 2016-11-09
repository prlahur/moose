# 3 primary species, chemistry time derivative
# vanGenuchten, unsaturated
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0.1 -2 3'
[]

[Variables]
  [./pp]
  [../]
  [./gc0]
  [../]
  [./gc1]
  [../]
  [./gc2]
  [../]
[]

[ICs]
  [./pp]
    type = RandomIC
    variable = pp
    min = -1
    max = 0
  [../]
  [./gc0]
    type = RandomIC
    variable = gc0
    min = 0
    max = 1
  [../]
  [./gc1]
    type = RandomIC
    variable = gc1
    min = 0
    max = 1
  [../]
  [./gc2]
    type = RandomIC
    variable = gc2
    min = 0
    max = 1
  [../]
[]

[Kernels]
  [./dummy_pp]
    type = TimeDerivative
    variable = pp
  [../]
  [./gc0]
    type = PorousFlowChemistryConvection
    variable = gc0
    primary_species = 0
  [../]
  [./gc1]
    type = PorousFlowChemistryConvection
    variable = gc1
    primary_species = 1
  [../]
  [./gc2]
    type = PorousFlowChemistryConvection
    variable = gc2
    primary_species = 2
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'gc0 pp gc1 gc2'
    number_fluid_phases = 1
    number_fluid_components = 1
    number_primary_species = 3
  [../]
[]

[Materials]
  [./dens0]
    type = PorousFlowDensityConstBulk
    density_P0 = 1
    bulk_modulus = 1.5
    phase = 0
  [../]
  [./dens_all]
    type = PorousFlowJoiner
    include_old = true
    material_property = PorousFlow_fluid_phase_density
  [../]
  [./dens_qp_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_fluid_phase_density_qp
    at_qps = true
  [../]
  [./visc0]
    type = PorousFlowViscosityConst
    viscosity = 1
    phase = 0
  [../]
  [./visc_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_viscosity
  [../]
  [./temperature]
    type = PorousFlowTemperature
  [../]
  [./nnn]
    type = PorousFlowNodeNumber
    on_initial_only = true
  [../]
  [./ppss]
    type = PorousFlow1PhaseP_VG
    porepressure = pp
    al = 1
    m = 0.5
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
  [../]
  [./concs]
    type = PorousFlowGeneralisedConcentrations
    generalised_conc_vars = 'gc0 gc1 gc2'
  [../]
  [./relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
  [../]
  [./relperm_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_relative_permeability
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1 0 0 0 2 0 0 0 3'
  [../]
[]

[Preconditioning]
  active = check
  [./andy]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'bcgs bjacobi 1E-15 1E-10 10000'
  [../]
  [./check]
    type = SMP
    full = true
    petsc_options = ''
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -snes_type'
    petsc_options_value = 'bcgs bjacobi 1E-15 1E-10 10000 test'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 1
  end_time = 2
[]

[Outputs]
  exodus = false
[]
