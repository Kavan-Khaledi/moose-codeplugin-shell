# Constant bending benchmark for a TRI6 shell cantilever.
# Prescribed tip rotation: 0.05 rad.
# Analytical: theta = M L / (E I), w = M L^2 / (2 E I).
# E = 2e11 Pa, I = b h^3 / 12 ~= 2e-4 m^4, L = 10 m.
# Expected bending moment = 2e5 N*m; tip displacement = 0.25 m.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 1
  xmin = 0.0
  xmax = 10
  zmin = 0.0
  zmax = 1
  elem_type = TRI6
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
  [fix_left_rotation]
    type = DirichletBC
    variable = rot_y
    boundary = left
    value = 0
  []
  [fix_left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_left_z]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0
  []
  [tip_rotation]
    type = DirichletBC
    variable = rot_y
    boundary = right
    value = 0.05
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
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
  line_search = none
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-5
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
    youngs_modulus = 2e11
    poissons_ratio = 0
    through_thickness_order = SECOND
  []
  [strain]
    type = ADComputeIncrementalTri6ShellStrain
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y'
    thickness = 0.133887
    through_thickness_order = SECOND
  []
  [stress]
    type = ADComputeShellStress
    through_thickness_order = SECOND
  []
[]

[Postprocessors]
  [tip_disp_z]
    type = PointValue
    point = '10 0 0'
    variable = disp_z
  []
[]

[AuxVariables]
  [moment_22]
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
  [moment_22]
    type = ShellResultantsAux
    variable = moment_22
    stress_resultant = bending_moment_1
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
