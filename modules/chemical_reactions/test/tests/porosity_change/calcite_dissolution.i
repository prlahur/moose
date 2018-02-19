# Porosity increase due to calcite dissolution.
# Calcium (Ca2+) and bicarbonate (HCO3-) react to form calcite (CaCO3) via
# the kinetic reaction
#
# Ca2+ + HCO3- = H+ + CaCO3(s)
#
# with reactive surface are A = 0.461 m^2/L, kinetic rate constant
# k = 6.456542e-2 mol/m^2 s, equilibrium constant Keq = 10^(1.8487)
#
# Note: Calcite molar volume is 36.934e-6 m^3/mol, but a greatly exaggerated
# value of 36.92e1 is used to exaggerate the size of porosity decrease in this test

[Mesh]
  type = GeneratedMesh
  dim = 2
[]

[Variables]
  [./ca2+]
    initial_condition = 1e-5
  [../]
  [./h+]
    initial_condition = 1.0e-6
  [../]
  [./hco3-]
    initial_condition = 1.0e-5
  [../]
[]

[AuxVariables]
  [./porosity]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./caco3_s]
    initial_condition = 8.1e-6
  [../]
  [./mineral_frac]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./porosity]
    type = MaterialRealAux
    variable = porosity
    property = porosity
    execute_on = 'initial linear'
  [../]
  [./mineral_frac]
    type = MaterialRealAux
    variable = mineral_frac
    property = mineral_volume_frac
    execute_on = 'initial linear'
  [../]
[]

[ReactionNetwork]
  primary_species = 'ca2+ hco3- h+'
  [./SolidKineticReactions]
    primary_species = 'ca2+ hco3- h+'
    kin_reactions = 'ca2+ + hco3- - h+ = caco3_s'
    secondary_species = 'caco3_s'
    log10_keq = '1.8487'
    reference_temperature = '298.15'
    system_temperature = '298.15'
    specific_reactive_surface_area = '0.1'
    kinetic_rate_constant = '6.456542e-7'
    activation_energy = '1.5e4'
  [../]
[]

[Kernels]
  [./ca2+_ie]
    type = PrimaryTimeDerivative
    variable = ca2+
  [../]
  [./h+_ie]
    type = PrimaryTimeDerivative
    variable = h+
  [../]
  [./hco3-_ie]
    type = PrimaryTimeDerivative
    variable = hco3-
  [../]
[]

[Materials]
  [./mineral_fraction]
    type = MineralVolumeFraction
    mineral_species = 'caco3_s'
    molar_volume = '36.934e2'
  [../]
  [./porous]
    type = SolidKineticPorosity
    base_porosity = 0.2
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  nl_abs_tol = 1e-10
  end_time = 100
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./porosity]
    type = ElementAverageValue
    variable = 'porosity'
    execute_on = 'initial timestep_end'
  [../]
  [./mineral_frac]
    type = ElementAverageValue
    variable = 'mineral_frac'
    execute_on = 'initial timestep_end'
  [../]
  [./caco3_s]
    type = ElementIntegralVariablePostprocessor
    variable = 'caco3_s'
    execute_on = 'initial timestep_end'
  [../]
[]

[Outputs]
  print_perf_log = true
  csv = true
[]
