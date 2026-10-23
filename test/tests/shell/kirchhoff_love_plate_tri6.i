# Verification of TRI6 shell elements using a uniformly loaded plate.
# The 4 x 4 m model represents half of an 8 x 4 m simply supported plate.
# Analytical reference: Navier solution (Kirchhoff-Love plate theory).
#
# Quantity       Analytical       TRI6 (MOOSE)      Error
# w_max          8.3896 mm        8.5900 mm         2.39%
# M11_max        390.46 kN*m/m     387.20 kN*m/m      0.83%
# M22_max        178.85 kN*m/m     178.17 kN*m/m      0.38%
#
# Moment components follow the local shell-axis convention.
# Good agreement with analytical results (maximum error < 2.5%).

[Mesh]
    [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    xmax = 4
    ymax = 4
    nx = 6
    ny = 6
    elem_type = TRI6
  []
[]

[Variables]
  [disp_x]
    order = SECOND
    family = LAGRANGE
  []
  [disp_y]
    order = SECOND
    family = LAGRANGE
  []
  [disp_z]
    order = SECOND
    family = LAGRANGE
  []
  [rot_x]
    order = SECOND
    family = LAGRANGE
  []
  [rot_y]
    order = SECOND
    family = LAGRANGE
  []
[]

[BCs]
  [symm_left_rot]
    type = DirichletBC
    variable = rot_y
    boundary = 'right' #symmetric BC
    value = 0.0
  []

  [simply_support_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'right top bottom left'
    value = 0.0
  []
  [simply_support_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'top bottom left'
    value = 0.0
  []
  [simply_support_z]
    type = DirichletBC
    variable = disp_z
    boundary = ' bottom top left '
    value = 0.0
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options = '-snes_ksp_ew'

  # best overall
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = ' lu       mumps'
  line_search = 'none'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  dt = 1
  end_time = 1
[]

[Kernels]
  [solid_disp_x]
    type = ADStressDivergenceTri6Shell

    component = 0
    variable = disp_x
    through_thickness_order = SECOND
  []
  [solid_disp_y]
    type = ADStressDivergenceTri6Shell

    component = 1
    variable = disp_y
    through_thickness_order = SECOND
  []
  [solid_disp_z]
    type = ADStressDivergenceTri6Shell

    component = 2
    variable = disp_z
    through_thickness_order = SECOND
  []
  [solid_rot_x]
    type = ADStressDivergenceTri6Shell

    component = 3
    variable = rot_x
    through_thickness_order = SECOND
  []
  [solid_rot_y]
    type = ADStressDivergenceTri6Shell

    component = 4
    variable = rot_y
    through_thickness_order = SECOND
  []

  [self_weight]
    type = ShellDistributedLoad
    value = 240000
    variable = disp_z
    function = 1
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensorShell
    youngs_modulus = 30e9
    poissons_ratio = 0.3
    through_thickness_order = SECOND
  []
  [strain]
    type = ADComputeIncrementalTri6ShellStrain
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y'
    thickness = 0.3
    reference_first_local_direction = '1 0 0'
    through_thickness_order = SECOND
  []
  [stress]
    type = ADComputeShellStress
    through_thickness_order = SECOND
  []
[]

[AuxVariables]
  [stress_00]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_11]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_22]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_01]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_02]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_12]
    order = CONSTANT
    family = MONOMIAL
  []
  [force_1]
    order = CONSTANT
    family = MONOMIAL
  []
  [force_2]
    order = CONSTANT
    family = MONOMIAL
  []
  [moment_11]
    order = CONSTANT
    family = MONOMIAL
  []
  [moment_22]
    order = CONSTANT
    family = MONOMIAL
  []
  [moment_12]
    order = CONSTANT
    family = MONOMIAL
  []
  [shear_12]
    order = CONSTANT
    family = MONOMIAL
  []
  [shear_13]
    order = CONSTANT
    family = MONOMIAL
  []
  [shear_23]
    order = CONSTANT
    family = MONOMIAL
  []
  [first_axis_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [first_axis_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [first_axis_z]
    order = CONSTANT
    family = MONOMIAL
  []
  [second_axis_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [second_axis_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [second_axis_z]
    order = CONSTANT
    family = MONOMIAL
  []
  [normal_axis_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [normal_axis_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [normal_axis_z]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [stress_00]
    type = RankTwoAux
    variable = stress_00
    rank_two_tensor = local_stress_t_points_0
    index_i = 0
    index_j = 0
    execute_on = TIMESTEP_END
  []
  [stress_11]
    type = RankTwoAux
    variable = stress_11
    rank_two_tensor = local_stress_t_points_0
    index_i = 1
    index_j = 1
    execute_on = TIMESTEP_END
  []
  [stress_22]
    type = RankTwoAux
    variable = stress_22
    rank_two_tensor = local_stress_t_points_0
    index_i = 2
    index_j = 2
    execute_on = TIMESTEP_END
  []
  [stress_01]
    type = RankTwoAux
    variable = stress_01
    rank_two_tensor = local_stress_t_points_0
    index_i = 0
    index_j = 1
    execute_on = TIMESTEP_END
  []
  [stress_02]
    type = RankTwoAux
    variable = stress_02
    rank_two_tensor = local_stress_t_points_0
    index_i = 0
    index_j = 2
    execute_on = TIMESTEP_END
  []
  [stress_12]
    type = RankTwoAux
    variable = stress_12
    rank_two_tensor = local_stress_t_points_0
    index_i = 1
    index_j = 2
    execute_on = TIMESTEP_END
  []

  [force_1]
    type = ShellResultantsAux
    variable = force_1
    stress_resultant = axial_force_0
    thickness = 0.3
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [force_2]
    type = ShellResultantsAux
    variable = force_2
    stress_resultant = axial_force_1
    thickness = 0.3
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [moment_11]
    type = ShellResultantsAux
    variable = moment_11
    stress_resultant = bending_moment_0
    thickness = 0.3
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [moment_22]
    type = ShellResultantsAux
    variable = moment_22
    stress_resultant = bending_moment_1
    thickness = 0.3
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [moment_12]
    type = ShellResultantsAux
    variable = moment_12
    stress_resultant = bending_moment_01
    thickness = 0.3
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [shear_12]
    type = ShellResultantsAux
    variable = shear_12
    stress_resultant = shear_force_01
    thickness = 0.3
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [shear_13]
    type = ShellResultantsAux
    variable = shear_13
    stress_resultant = shear_force_02
    thickness = 0.3
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [shear_23]
    type = ShellResultantsAux
    variable = shear_23
    stress_resultant = shear_force_12
    thickness = 0.3
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []

  [first_axis_x]
    type = ShellLocalCoordinatesAux
    variable = first_axis_x
    property = first_local_vector
    component = 0
  []
  [first_axis_y]
    type = ShellLocalCoordinatesAux
    variable = first_axis_y
    property = first_local_vector
    component = 1
  []
  [first_axis_z]
    type = ShellLocalCoordinatesAux
    variable = first_axis_z
    property = first_local_vector
    component = 2
  []

  [second_axis_x]
    type = ShellLocalCoordinatesAux
    variable = second_axis_x
    property = second_local_vector
    component = 0
  []
  [second_axis_y]
    type = ShellLocalCoordinatesAux
    variable = second_axis_y
    property = second_local_vector
    component = 1
  []
  [second_axis_z]
    type = ShellLocalCoordinatesAux
    variable = second_axis_z
    property = second_local_vector
    component = 2
  []

  [normal_axis_x]
    type = ShellLocalCoordinatesAux
    variable = normal_axis_x
    property = normal_local_vector
    component = 0
  []
  [normal_axis_y]
    type = ShellLocalCoordinatesAux
    variable = normal_axis_y
    property = normal_local_vector
    component = 1
  []
  [normal_axis_z]
    type = ShellLocalCoordinatesAux
    variable = normal_axis_z
    property = normal_local_vector
    component = 2
  []
[]

[Outputs]
  exodus = true
[]
