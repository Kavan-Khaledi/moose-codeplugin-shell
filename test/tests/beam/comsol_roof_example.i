# based on the existing analytical Solutions, maximum deflection of the roof is 0.3086

[Mesh]
  [file]
    type = FileMeshGenerator
    file = comsol_roof_example.msh
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
    variable = rot_x
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
  nl_abs_tol = 1e-10
  dt = 1
  end_time = 1
[]

[Kernels]

  [solid_disp_x]
    type = ADStressDivergenceTriangularShell
    component = 0
    variable = disp_x
    through_thickness_order = SECOND
    block = 'shell'
  []
  [solid_disp_y]
    type = ADStressDivergenceTriangularShell
    component = 1
    variable = disp_y
    through_thickness_order = SECOND
    block = 'shell'
  []
  [solid_disp_z]
    type = ADStressDivergenceTriangularShell
    component = 2
    variable = disp_z
    through_thickness_order = SECOND
    block = 'shell'
  []
  [solid_rot_x]
    type = ADStressDivergenceTriangularShell
    component = 3
    variable = rot_x
    through_thickness_order = SECOND
    block = 'shell'
  []
  [solid_rot_y]
    type = ADStressDivergenceTriangularShell
    component = 4
    variable = rot_y
    through_thickness_order = SECOND
    block = 'shell'
  []

  [self_weight]
    type = ShellDistributedLoad
    value = 90
    variable = disp_z
    function = 1
    block = 'shell'
  []
[]

[Materials]

  [elasticity_tshell]
    type = ADComputeIsotropicElasticityTensorTriangularShell
    youngs_modulus = 4.32e8
    poissons_ratio = 0.0
    through_thickness_order = SECOND
    block = 'shell'
  []
  [strain_shell]
    type = ADComputeIncrementalTriangularShellStrain
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y'
    thickness = 0.25
    through_thickness_order = SECOND
    #user_defined_first_local_vector=' 1 0 0'
    block = 'shell'
  []
  [stress_shell]
    type = ADComputeTriangularShellStress
    through_thickness_order = SECOND
    thickness = 0.25
    block = 'shell'
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
  [moment_12]
    type = ADMaterialRealAux
    variable = moment_12
    property = shell_moment_12
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
