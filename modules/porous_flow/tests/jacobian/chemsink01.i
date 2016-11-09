# 1 primary species, chemistry sink
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
    coefficient = 1
    Ceq = 1
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'gc0 pp'
    number_fluid_phases = 1
    number_fluid_components = 1
    number_primary_species = 1
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
    generalised_conc_vars = 'gc0'
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
    #petsc_options = '-snes_test_display'
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
