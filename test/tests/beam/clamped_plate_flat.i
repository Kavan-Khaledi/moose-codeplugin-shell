# Test for simply supported plate under uniform pressure

# One quarter of a 50 m x 50 m x 1m plate is modeled in this test.
# Pressure loading is applied on the top surface using nodal forces
# of magnitude -10 N on all nodes. This corresponds to a pressure (q) of
# -10.816 N/m^2.

# The FEM solution at (0,0), which is at the center of the full plate
# is -2.997084e-03 m.

# The analytical solution for displacement at center of plate obtained
# using a thin plate assumption for a square plate is
# w = 16 q a^4/(D*pi^6) \sum_{m = 1,3,5, ..}^\inf \sum_{n = 1,3,5, ..}^\inf  (-1)^{(m+n-2)/2}/(mn*(m^2+n^2)^2)

# The above solution is the Naviers series solution from the "Theory of plates
# and shells" by Timoshenko and Woinowsky-Krieger (1959).

# where a = 50 m, q = -10.816 N/m^2 and D = E/(12(1-v^2))
# The analytical series solution converges to 2.998535904e-03 m
# when the first 16 terms of the series are considered (i.e., until
# m & n = 7).

# The resulting relative error between FEM and analytical solution is
# 0.048%.

[Mesh]
  [gmg]
    type = FileMeshGenerator
    file = clamped_plate_flat.msh
  []
  second_order = true
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
    boundary = 'left'
    value = 0.0
  []
  [symm_bottom_rot]
    type = DirichletBC
    variable = rot_x
    boundary = 'bottom'
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
    boundary = 'right top bottom left'
    value = 0.0
  []
  [simply_support_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'right top'
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
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-5
  dt = 1
  dtmin = 1
  end_time = 1.
[]

[Kernels]
  [solid_disp_x]
    type = ADStressDivergenceTriangularShell

    component = 0
    variable = disp_x
    through_thickness_order = SECOND
  []
  [solid_disp_y]
    type = ADStressDivergenceTriangularShell

    component = 1
    variable = disp_y
    through_thickness_order = SECOND
  []
  [solid_disp_z]
    type = ADStressDivergenceTriangularShell

    component = 2
    variable = disp_z
    through_thickness_order = SECOND
  []
  [solid_rot_x]
    type = ADStressDivergenceTriangularShell

    component = 3
    variable = rot_x
    through_thickness_order = SECOND
  []
  [solid_rot_y]
    type = ADStressDivergenceTriangularShell

    component = 4
    variable = rot_y
    through_thickness_order = SECOND
  []

  [self_weight]
    type = ShellDistributedLoad
    value = 10.816
    variable = disp_z
    function = 1
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensorTriangularShell
    youngs_modulus = 1e9
    poissons_ratio = 0.3

    through_thickness_order = SECOND
  []
  [strain]
    type = ADComputeIncrementalTriangularShellStrain
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y'
    thickness = 1
    #user_defined_first_local_vector=' 1 0 0'
    through_thickness_order = SECOND
  []
  [stress]
    type = ADComputeTriangularShellStress
    thickness = 1
    through_thickness_order = SECOND
  []
[]

[Postprocessors]
  [disp_z2]
    type = PointValue
    point = '0.0 0.0 0.0'
    variable = disp_z
  []
[]

[AuxVariables]
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
  [force_1]
    type = ADMaterialRealAux
    variable = force_1
    property = shell_force_1
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [force_2]
    type = ADMaterialRealAux
    variable = force_2
    property = shell_force_2
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [moment_11]
    type = ADMaterialRealAux
    variable = moment_11
    property = shell_moment_11
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [moment_22]
    type = ADMaterialRealAux
    variable = moment_22
    property = shell_moment_22
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [shear_12]
    type = ADMaterialRealAux
    variable = shear_12
    property = shell_shear_12
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [shear_13]
    type = ADMaterialRealAux
    variable = shear_13
    property = shell_shear_13
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [shear_23]
    type = ADMaterialRealAux
    variable = shear_23
    property = shell_shear_23
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [first_axis_x]
    type = MaterialRealAux
    variable = first_axis_x
    property = first_local_axis_x
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [first_axis_y]
    type = MaterialRealAux
    variable = first_axis_y
    property = first_local_axis_y
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [first_axis_z]
    type = MaterialRealAux
    variable = first_axis_z
    property = first_local_axis_z
    execute_on = TIMESTEP_END
    block = 'shell'
  []

  [second_axis_x]
    type = MaterialRealAux
    variable = second_axis_x
    property = second_local_axis_x
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [second_axis_y]
    type = MaterialRealAux
    variable = second_axis_y
    property = second_local_axis_y
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [second_axis_z]
    type = MaterialRealAux
    variable = second_axis_z
    property = second_local_axis_z
    execute_on = TIMESTEP_END
    block = 'shell'
  []

  [normal_axis_x]
    type = MaterialRealAux
    variable = normal_axis_x
    property = normal_local_axis_x
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [normal_axis_y]
    type = MaterialRealAux
    variable = normal_axis_y
    property = normal_local_axis_y
    execute_on = TIMESTEP_END
    block = 'shell'
  []
  [normal_axis_z]
    type = MaterialRealAux
    variable = normal_axis_z
    property = normal_local_axis_z
    execute_on = TIMESTEP_END
    block = 'shell'
  []
[]
[Outputs]
  exodus = true
[]
