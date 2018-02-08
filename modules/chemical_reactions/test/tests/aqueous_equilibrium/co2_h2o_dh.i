# Batch CO2 - H2O equilibrium reaction at 25C using extended Debye-Huckel
# activity coefficients.
#
# Aqueous equilibrium reactions:
# a)  H+ + HCO3- = CO2(aq),         Keq = 10^(6.3447)
# b)  HCO3- = H+ + CO3--,           Keq = 10^(-10.3288)
# c)  - H+ = OH-,                   Keq = 10^(-13.9951)
#
# The primary chemical species are h+ and hco3-, and the secondary equilibrium
# species are co2(aq), co3-- and oh-.
#
# The charge z and ion radii a (in Angstrom) are:
# H+: z = 1, a = 9
# HCO3-: z = -1, a = 4
# OH-: z = -1, a = 3.5
# CO3-: z = -1, a = 4.5
# CO2(aq): z = 0, a = 3

[Mesh]
  type = GeneratedMesh
  dim = 2
[]

[AuxVariables]
  [./oh-]
  [../]
  [./co3--]
  [../]
  [./co2_aq]
  [../]
  [./ph]
  [../]
  [./ionic_strength]
  [../]
  [./gamma_h+]
  [../]
  [./gamma_hco3-]
  [../]
  [./gamma_oh-]
  [../]
  [./gamma_co3--]
  [../]
  [./gamma_co2_aq]
  [../]
[]

[AuxKernels]
  [./ph]
    type = PHAux
    variable = ph
    h_conc = h+
    activity_coeff = gamma_h+
    execute_on = 'initial linear'
  [../]
  [./oh-]
    type = AqueousEquilibriumRxnAux
    variable = oh-
    v = h+
    log_k = -13.9951
    sto_v = -1
    execute_on = 'initial linear'
  [../]
  [./co3--]
    type = AqueousEquilibriumRxnAux
    variable = co3--
    v = 'hco3- h+'
    log_k = -10.3288
    sto_v = '1 -1'
    execute_on = 'initial linear'
  [../]
  [./co2_aq]
    type = AqueousEquilibriumRxnAux
    variable = co2_aq
    v = 'hco3- h+'
    log_k = 6.3447
    sto_v = '1 1'
    execute_on = 'initial linear'
  [../]
  [./ionic_strength]
    type = IonicStrengthAux
    conc = 'h+ hco3- oh- co3-- co2_aq'
    z = '1 -1 -1 -1 0'
    variable = ionic_strength
    execute_on = 'initial linear'
  [../]
  [./gamma_h+]
    type = ActivityCoefficientDHAux
    variable = gamma_h+
    a = 9
    z = 1
    ionic_strength = ionic_strength
    execute_on = 'initial linear'
  [../]
  [./gamma_hco3-]
    type = ActivityCoefficientDHAux
    variable = gamma_hco3-
    a = 4
    z = -1
    ionic_strength = ionic_strength
    execute_on = 'initial linear'
  [../]
  [./gamma_oh-]
    type = ActivityCoefficientDHAux
    variable = gamma_oh-
    a = 3.5
    z = -1
    ionic_strength = ionic_strength
    execute_on = 'initial linear'
  [../]
  [./gamma_co3--]
    type = ActivityCoefficientDHAux
    variable = gamma_co3--
    a = 4.5
    z = -1
    ionic_strength = ionic_strength
    execute_on = 'initial linear'
  [../]
  [./gamma_co2_aq]
    type = ActivityCoefficientDHAux
    variable = gamma_co2_aq
    a = 3
    z = 0
    ionic_strength = ionic_strength
    execute_on = 'initial linear'
  [../]
[]

[Variables]
  [./h+]
    initial_condition = 1e-5
  [../]
  [./hco3-]
    initial_condition = 1e-5
  [../]
[]

[Kernels]
  [./h+_td]
    type = PrimaryTimeDerivative
    variable = h+
  [../]
  [./hco3-_td]
    type = PrimaryTimeDerivative
    variable = hco3-
  [../]
  [./h+_oh-_td]
    type = CoupledBEEquilibriumSub
    variable = h+
    sto_u = -1
    log_k = -13.9951
    sto_v = ' '
    gamma_eq = gamma_oh-
    gamma_u = gamma_h+
  [../]
  [./h+_oh-_diff]
    type = CoupledDiffusionReactionSub
    variable = h+
    sto_u = -1
    log_k = -13.9951
    sto_v = ' '
    gamma_eq = gamma_oh-
    gamma_u = gamma_h+
  [../]
  [./h+_co3--_td]
    type = CoupledBEEquilibriumSub
    variable = h+
    v = hco3-
    sto_v = 1
    sto_u = -1
    log_k = -10.3288
    gamma_eq = gamma_co3--
    gamma_u = gamma_h+
    gamma_v = gamma_hco3-
  [../]
  [./h+_co3--_diff]
    type = CoupledDiffusionReactionSub
    variable = h+
    v = hco3-
    sto_v = 1
    sto_u = -1
    log_k = -10.3288
    gamma_eq = gamma_co3--
    gamma_u = gamma_h+
    gamma_v = gamma_hco3-
  [../]
  [./h+_co2_aq_td]
    type = CoupledBEEquilibriumSub
    variable = h+
    v = hco3-
    sto_v = 1
    sto_u = 1
    log_k = 6.3447
    gamma_eq = gamma_co2_aq
    gamma_u = gamma_h+
    gamma_v = gamma_hco3-
  [../]
  [./h+_co2_aq_diff]
    type = CoupledDiffusionReactionSub
    variable = h+
    v = hco3-
    sto_v = 1
    sto_u = 1
    log_k = 6.3447
    gamma_eq = gamma_co2_aq
    gamma_u = gamma_h+
    gamma_v = gamma_hco3-
  [../]
  [./hco3-_co3--_td]
    type = CoupledBEEquilibriumSub
    variable = hco3-
    v = h+
    sto_v = -1
    sto_u = 1
    log_k = -10.3288
    gamma_eq = gamma_co3--
    gamma_u = gamma_hco3-
    gamma_v = gamma_h+
  [../]
  [./hco3-_co3--_diff]
    type = CoupledDiffusionReactionSub
    variable = hco3-
    v = h+
    sto_v = -1
    sto_u = 1
    log_k = -10.3288
    gamma_eq = gamma_co3--
    gamma_u = gamma_hco3-
    gamma_v = gamma_h+
  [../]
  [./hco3-_co2_aq_td]
    type = CoupledBEEquilibriumSub
    variable = hco3-
    v = h+
    sto_v = 1
    sto_u = 1
    log_k = 6.3447
    gamma_eq = gamma_co2_aq
    gamma_u = gamma_hco3-
    gamma_v = gamma_h+
  [../]
  [./hco3-_co2_aq_diff]
    type = CoupledDiffusionReactionSub
    variable = hco3-
    v = h+
    sto_v = 1
    sto_u = 1
    log_k = 6.3447
    gamma_eq = gamma_co2_aq
    gamma_u = gamma_hco3-
    gamma_v = gamma_h+
  [../]
[]

[Materials]
  [./porous]
    type = GenericConstantMaterial
    prop_names = 'diffusivity porosity'
    prop_values = '1e-7 0.25'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  nl_abs_tol = 1e-12
  end_time = 1
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./h+]
    type = ElementIntegralVariablePostprocessor
    variable = h+
    execute_on = 'initial timestep_end'
  [../]
  [./hco3-]
    type = ElementIntegralVariablePostprocessor
    variable = hco3-
    execute_on = 'initial timestep_end'
  [../]
  [./co2_aq]
    type = ElementIntegralVariablePostprocessor
    variable = co2_aq
    execute_on = 'initial timestep_end'
  [../]
  [./oh-]
    type = ElementIntegralVariablePostprocessor
    variable = oh-
    execute_on = 'initial timestep_end'
  [../]
  [./co3--]
    type = ElementIntegralVariablePostprocessor
    variable = co3--
    execute_on = 'initial timestep_end'
  [../]
  [./ph]
    type = ElementIntegralVariablePostprocessor
    variable = ph
    execute_on = 'initial timestep_end'
  [../]
  [./I]
    type = ElementIntegralVariablePostprocessor
    variable = ionic_strength
    execute_on = 'initial timestep_end'
  [../]
  [./gamma_h+]
    type = ElementIntegralVariablePostprocessor
    variable = gamma_h+
    execute_on = 'initial timestep_end'
  [../]
  [./gamma_hco3-]
    type = ElementIntegralVariablePostprocessor
    variable = gamma_hco3-
    execute_on = 'initial timestep_end'
  [../]
  [./gamma_oh-]
    type = ElementIntegralVariablePostprocessor
    variable = gamma_oh-
    execute_on = 'initial timestep_end'
  [../]
[]

[Outputs]
  print_perf_log = true
  csv = true
[]
