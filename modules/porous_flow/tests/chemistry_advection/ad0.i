# Jonathon's eqns 229 and 230
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
  primary_species = 0
  phase = 0
[]

[Variables]
  [./gc]
  [../]
  [./phis]
    initial_condition = 0.1
  [../]
[]


[BCs]
  [./gc]
    type = PresetBC
    boundary = left
    variable = gc
    value = 1
  [../]
[]

[Kernels]
  [./gc_time]
    type = PorousFlowChemistryTimeDerivative
    variable = gc
  [../]
  [./gc_conv]
    type = PorousFlowChemistryConvection
    variable = gc
  [../]
  [./gc_diff]
    type = PorousFlowChemistryDispersiveFlux
    disp_long = 0.2
    disp_trans = 0.2
    variable = gc
  [../]
  [./gc_source]
    type = PorousFlowChemistrySink
    variable = gc
    coefficient = 1E-1
  [../]
  [./phist]
    type = TimeDerivative
    variable = phis
  [../]
  [./phis_source]
    type = PorousFlowChemistrySink
    variable = phis
    coefficient = -1E-2
  [../]
[]

[AuxVariables]
  [./pp]
  [../]
[]

[AuxKernels]
  [./pp]
    type = FunctionAux
    variable = pp
    function = '1-x'
    execute_on = initial
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'gc'
    number_fluid_phases = 1
    number_fluid_components = 1
    number_primary_species = 1
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
     viscosity = 1E3
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
    generalised_conc_vars = gc
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
  [./porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  [../]
  [./diffusion]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = 0.2
    tortuosity = 0.2
  [../]
[]

[Preconditioning]
  [./andy]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'bcgs bjacobi 1E-15 1E-10 10000'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 0.5
  end_time = 10
[]

[Outputs]
  exodus = true
[]
