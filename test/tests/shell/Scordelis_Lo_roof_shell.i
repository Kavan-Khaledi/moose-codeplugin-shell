# This model is a widely used benchmark model denoted the Scordelis-Lo roof.
# The maximum z-deformation is compared with the value given in "Proposed Standard Set of Problems to Test Finite Element Accuracy, Finite Elements in Analysis and Design, 1985".
# Based on the existing analytical Solutions, maximum deflection of the roof should be 0.3086
# The model results in a maximum vertical deflection of 0.3068 (0.5% error) using 6-noded shell elements.

[Mesh]
  [file]
    type = FileMeshGenerator
    file = Scordelis_Lo_roof_shell_tri6.msh
  []

  second_order = true
[]

[Variables]
  [disp_x]
    order = SECOND
    family = LAGRANGE
    block = 'shell'
  []
  [disp_y]
    order = SECOND
    family = LAGRANGE
    block = 'shell'
  []
  [disp_z]
    order = SECOND
    family = LAGRANGE
    block = 'shell'
  []
  [rot_x]
    order = SECOND
    family = LAGRANGE
    block = 'shell'
  []
  [rot_y]
    order = SECOND
    family = LAGRANGE
    block = 'shell'
  []
[]

[BCs]

  [simply_support_y]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'back'
    value = 0.0
  []
  [simply_support_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  []

  [simply_support_x]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'front'
    value = 0.0
  []

  [simply_rot_x]
    type = ADDirichletBC
    variable = rot_y
    boundary = 'front'
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
  nl_abs_tol = 1e-7
  dt = 1
  end_time = 1
[]

[Kernels]

  [solid_disp_x]
    type = ADStressDivergenceTri6Shell
    component = 0
    variable = disp_x
    through_thickness_order = SECOND
    block = 'shell'
  []
  [solid_disp_y]
    type = ADStressDivergenceTri6Shell
    component = 1
    variable = disp_y
    through_thickness_order = SECOND
    block = 'shell'
  []
  [solid_disp_z]
    type = ADStressDivergenceTri6Shell
    component = 2
    variable = disp_z
    through_thickness_order = SECOND
    block = 'shell'
  []
  [solid_rot_x]
    type = ADStressDivergenceTri6Shell
    component = 3
    variable = rot_x
    through_thickness_order = SECOND
    block = 'shell'
  []
  [solid_rot_y]
    type = ADStressDivergenceTri6Shell
    component = 4
    variable = rot_y
    through_thickness_order = SECOND
    block = 'shell'
  []

  [self_weight]
    type = ADDistributedLoadShell
    variable = disp_z
    function = 90
    block = 'shell'
    displacements= 'disp_x disp_y disp_z'
  []
[]

[Materials]

  [elasticity_tshell]
    type = ADComputeIsotropicElasticityTensorShell
    youngs_modulus = 4.32e8
    poissons_ratio = 0.0
    through_thickness_order = SECOND
    block = 'shell'
  []
  [strain_shell]
    type =ADComputeIncrementalTri6ShellStrain
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y'
    thickness = 0.25
    through_thickness_order = SECOND
    reference_first_local_direction = '0 1 0'
    block = 'shell'
  []
  [stress_shell]
    type = ADComputeShellStress
    through_thickness_order = SECOND
    block = 'shell'
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
    thickness = 0.01
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [force_2]
    type = ShellResultantsAux
    variable = force_2
    stress_resultant = axial_force_1
    thickness = 0.01
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [moment_11]
    type = ShellResultantsAux
    variable = moment_11
    stress_resultant = bending_moment_0
    thickness = 0.01
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [moment_22]
    type = ShellResultantsAux
    variable = moment_22
    stress_resultant = bending_moment_1
    thickness = 0.01
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []

  [moment_12]
    type = ShellResultantsAux
    variable = moment_12
    stress_resultant = bending_moment_01
    thickness = 0.01
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [shear_12]
    type = ShellResultantsAux
    variable = shear_12
    stress_resultant = shear_force_01
    thickness = 0.01
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [shear_13]
    type = ShellResultantsAux
    variable = shear_13
    stress_resultant = shear_force_02
    thickness = 0.01
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []
  [shear_23]
    type = ShellResultantsAux
    variable = shear_23
    stress_resultant = shear_force_12
    thickness = 0.01
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

[Postprocessors]
  [disp_z2]
    type = PointValue
    point = '-16.7 0  19.2'
    variable = disp_z
  []
[]

[Outputs]
  exodus = true
[]
