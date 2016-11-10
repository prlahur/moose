# Test the Jacobian of the PorousFlowChemistryDisperiveFlux kernel

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  xmin = -1
  xmax = 1
  ny = 1
  ymin = 0
  ymax = 1
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [./gc0]
  [../]
[]

[AuxVariables]
  [./pp]
  [../]
  [./mf]
  [../]
[]

[AuxKernels]
  [./pp]
    type = FunctionAux
    variable = pp
    function = x
  [../]
[]

[ICs]
  [./gc0]
    type = RandomIC
    variable = gc0
    min = 0
    max = 1
  [../]
[]

[Kernels]
  [./diff1]
    type = PorousFlowChemistryDispersiveFlux
    variable = gc0
    gravity = '1 0 0'
    disp_long = 0.2
    disp_trans = 0.1
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'gc0'
    number_fluid_phases = 1
    number_fluid_components = 2
    number_primary_species = 1
  [../]
[]

[Materials]
  [./nn]
    type = PorousFlowNodeNumber
  [../]
  [./temp]
    type = PorousFlowTemperature
  [../]
  [./ppss]
    type = PorousFlow1PhaseP
    porepressure = pp
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = mf
  [../]
  [./concs]
    type = PorousFlowGeneralisedConcentrations
    generalised_conc_vars = gc0
  [../]
  [./dens0]
    type = PorousFlowDensityConstBulk
    density_P0 = 10
    bulk_modulus = 1e7
    phase = 0
  [../]
  [./dens_qp_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_fluid_phase_density_qp
    at_qps = true
  [../]
  [./poro]
    type = PorousFlowPorosityConst
    porosity = 0.1
  [../]
  [./diff]
    type = PorousFlowDiffusivityConst
     diffusion_coeff = '1e-2 1e-1'
     tortuosity = '0.1'
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
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1 0 0 0 2 0 0 0 3'
  [../]
  [./relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  [../]
  [./relperm_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_relative_permeability
  [../]
[]

[Preconditioning]
  active = smp
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_test_display'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -snes_type'
    petsc_options_value = 'bcgs bjacobi 1E-15 1E-10 10000 test'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 1
  end_time = 1
[]

[Outputs]
  exodus = false
[]
