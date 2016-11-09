# 3 primary species, chemistry sink
# vanGenuchten, unsaturated
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
[]

[GlobalParams]
  PorousFlowDictator = dictator
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
    type = PorousFlowChemistrySink
    variable = gc0
    primary_species = 0
    phase = 0
    coefficient = 10
    Ceq = 1
  [../]
  [./gc1]
    type = PorousFlowChemistrySink
    variable = gc1
    primary_species = 1
    coefficient = -20
    Ceq = -2
  [../]
  [./gc2]
    type = PorousFlowChemistrySink
    variable = gc2
    primary_species = 2
    phase = 0
    coefficient = -20
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
  [./porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
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
