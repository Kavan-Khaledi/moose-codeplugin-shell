[Mesh]
  [gmg]
    type = FileMeshGenerator
    file = Flac3d_Isotropic_Plate.msh
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
  #dtmin = 1
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
    value = 240000
    variable = disp_z
    function = 1
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensorTriangularShell
    youngs_modulus = 30e9
    poissons_ratio = 0.3
    through_thickness_order = SECOND
  []
  [strain]
    type = ADComputeIncrementalTriangularShellStrain
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y'
    thickness = 0.3
    #user_defined_first_local_vector=' 1 0 0'
    through_thickness_order = SECOND
  []
  [stress]
    type = ADComputeTriangularShellStress
    thickness = 0.3
    through_thickness_order = SECOND
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
  []
  [force_2]
    type = ADMaterialRealAux
    variable = force_2
    property = shell_force_2
    execute_on = TIMESTEP_END
  []
  [moment_11]
    type = ADMaterialRealAux
    variable = moment_11
    property = shell_moment_11
    execute_on = TIMESTEP_END
  []
  [moment_22]
    type = ADMaterialRealAux
    variable = moment_22
    property = shell_moment_22
    execute_on = TIMESTEP_END
  []
  [shear_12]
    type = ADMaterialRealAux
    variable = shear_12
    property = shell_shear_12
    execute_on = TIMESTEP_END
  []
  [shear_13]
    type = ADMaterialRealAux
    variable = shear_13
    property = shell_shear_13
    execute_on = TIMESTEP_END
  []
  [shear_23]
    type = ADMaterialRealAux
    variable = shear_23
    property = shell_shear_23
    execute_on = TIMESTEP_END
  []
  [first_axis_x]
    type = MaterialRealAux
    variable = first_axis_x
    property = first_local_axis_x
    execute_on = TIMESTEP_END
  []
  [first_axis_y]
    type = MaterialRealAux
    variable = first_axis_y
    property = first_local_axis_y
    execute_on = TIMESTEP_END
  []
  [first_axis_z]
    type = MaterialRealAux
    variable = first_axis_z
    property = first_local_axis_z
    execute_on = TIMESTEP_END
  []

  [second_axis_x]
    type = MaterialRealAux
    variable = second_axis_x
    property = second_local_axis_x
    execute_on = TIMESTEP_END
  []
  [second_axis_y]
    type = MaterialRealAux
    variable = second_axis_y
    property = second_local_axis_y
    execute_on = TIMESTEP_END
  []
  [second_axis_z]
    type = MaterialRealAux
    variable = second_axis_z
    property = second_local_axis_z
    execute_on = TIMESTEP_END
  []

  [normal_axis_x]
    type = MaterialRealAux
    variable = normal_axis_x
    property = normal_local_axis_x
    execute_on = TIMESTEP_END
  []
  [normal_axis_y]
    type = MaterialRealAux
    variable = normal_axis_y
    property = normal_local_axis_y
    execute_on = TIMESTEP_END
  []
  [normal_axis_z]
    type = MaterialRealAux
    variable = normal_axis_z
    property = normal_local_axis_z
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  exodus = true
[]
