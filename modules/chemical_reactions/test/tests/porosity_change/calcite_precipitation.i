# Porosity increase due to calcite precipitation.
# Calcium (Ca2+) and bicarbonate (HCO3-) react to form calcite (CaCO3) via
# the kinetic reaction
#
# Ca2+ + HCO3- = H+ + CaCO3(s)
#
# with reactive surface are A = 0.461 m^2/L, kinetic rate constant
# k = 6.456542e-2 mol/m^2 s, equilibrium constant Keq = 10^(1.8487)
#
# Calcite molar volume is calculated as molar mass / density, where
# the molar mass is 100.0869 g/mol, and density is 2.711 g/cm^3,
# giving a molar volume of 36.92 cm^3/mol

[Mesh]
  type = GeneratedMesh
  dim = 2
  xmax = 2
[]

[Variables]
  [./ca2+]
    initial_condition = 0.1
  [../]
  [./h+]
    initial_condition = 1.0e-6
  [../]
  [./hco3-]
    initial_condition = 0.1
  [../]
[]

[AuxVariables]
  [./porosity]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./caco3_s]
    initial_condition = 0.1
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
    secondary_species = caco3_s
    log10_keq = 1.8487
    reference_temperature = 298.15
    system_temperature = 298.15
    gas_constant = 8.314
    specific_reactive_surface_area = 4.61e-4
    kinetic_rate_constant = 6.456542e-2
    activation_energy = 1.5e4
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
    mineral_species = caco3_s
    molar_volume = 36.92e-6
  [../]
  [./porous]
    type = SolidKineticPorosity
    base_porosity = 0.2
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  nl_abs_tol = 1e-12
  dt = 1e2
  end_time = 1e3
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
     variable = porosity
     execute_on = 'initial timestep_end'
  [../]
  [./mineral_frac]
    type = ElementAverageValue
    variable = mineral_frac
    execute_on = 'initial timestep_end'
  [../]
[]

[Outputs]
  print_perf_log = true
  exodus = true
[]
