# Verification of TRI6 shell elements in x-z plane under two concentrated loads.
# Simply supported plate: L = 9 m, width = 1 m, thickness = 0.133887 m.
# Two loads of 10 kN applied at x = 3 m and x = 6 m.
# Analytical reference: Euler-Bernoulli beam theory.
#
# Quantity              Analytical       TRI6 (MOOSE)     
# Center displacement   6.4684e-3 m      6.44e-3 m        
# Max. bending moment   29998.5 N*m/m    30000 N*m/m      
# Moment at center      29998.5 N*m/m    30000 N*m/m     
# Max. shear force      9999.5 N/m       10000 N/m       
# Shear at center       0 N/m            0 N/m            
# Max. membrane force   0 N/m            0 N/m            
#
# Each TRI6 edge has 3 nodes. For ny elements, the loading line
# contains (2*ny + 1) unique nodes, each receiving 10000/(2*ny+1) N.


ny = 3
force_per_node = ${fparse -10000 / (2 * ${ny} + 1)}
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    xmax = 9
    ymax= 1
    nx = 12
    ny = ${ny}
    elem_type = TRI6
  []


  [p1]
    type = BoundingBoxNodeSetGenerator
    input = gmg
    bottom_left = '2.99  -0.1 0.0'
    top_right = '3.1  1.1 0.0'
    new_boundary = 100
  []

  [p2]
    type = BoundingBoxNodeSetGenerator
    input = p1
    bottom_left = '5.99  -0.1 0.0'
    top_right = '6.1  1.1 0.0'
    new_boundary = 101
  []
  [rotate]
    type = TransformGenerator
    transform = ROTATE
    vector_value = '0 90 0'
    input =p2
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
  [simply_support_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'right left'
    value = 0.0
  []
  [simply_support_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'right left'
    value = 0.0
  []
[]

[NodalKernels]
  [force_p1]
    type = ConstantRate
    variable = disp_y
    boundary = 100
    rate = ${force_per_node}
  []

  [force_p2]
    type = ConstantRate
    variable = disp_y
    boundary = 101
    rate = ${force_per_node}
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

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu mumps'

  line_search = 'none'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-6
  dt = 1
  dtmin = 1
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
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensorShell
    youngs_modulus = 200e9
    poissons_ratio = 0.0
    through_thickness_order = SECOND
  []

  [strain]
    type = ADComputeIncrementalTri6ShellStrain
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y'
    thickness = 0.133887
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
    thickness = 0.133887
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []

  [force_2]
    type = ShellResultantsAux
    variable = force_2
    stress_resultant = axial_force_1
    thickness = 0.133887
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []

  [moment_11]
    type = ShellResultantsAux
    variable = moment_11
    stress_resultant = bending_moment_0
    thickness = 0.133887
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []

  [moment_22]
    type = ShellResultantsAux
    variable = moment_22
    stress_resultant = bending_moment_1
    thickness = 0.133887
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []

  [moment_12]
    type = ShellResultantsAux
    variable = moment_12
    stress_resultant = bending_moment_01
    thickness = 0.133887
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []

  [shear_12]
    type = ShellResultantsAux
    variable = shear_12
    stress_resultant = shear_force_01
    thickness = 0.133887
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []

  [shear_13]
    type = ShellResultantsAux
    variable = shear_13
    stress_resultant = shear_force_02
    thickness = 0.133887
    through_thickness_order = SECOND
    execute_on = TIMESTEP_END
  []

  [shear_23]
    type = ShellResultantsAux
    variable = shear_23
    stress_resultant = shear_force_12
    thickness = 0.133887
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
