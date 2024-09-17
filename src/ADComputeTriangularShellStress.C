//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeTriangularShellStress.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "ADComputeIsotropicElasticityTensorShell.h"
#include "ADComputeOrthotropicElasticityTensorShell.h"
#include "ADComputeIsotropicElasticityTensorTriangularShell.h"

#include "libmesh/quadrature.h"
#include "libmesh/utility.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_type.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature_gauss.h"

registerMooseObject(MOOSEAPPNAME, ADComputeTriangularShellStress);

InputParameters
ADComputeTriangularShellStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute in-plane stress using elasticity for shell");
  params.addRequiredParam<std::string>("through_thickness_order",
                                       "Quadrature order in out of plane direction");
  params.addRequiredCoupledVar(
      "thickness",
      "Thickness of the shell. Can be supplied as either a number or a variable name.");
  return params;
}

ADComputeTriangularShellStress::ADComputeTriangularShellStress(const InputParameters & parameters)
  : Material(parameters),

    _thickness(coupledValue("thickness"))
{
  // get number of quadrature points along thickness based on order
  std::unique_ptr<QGauss> t_qrule = std::make_unique<QGauss>(
      1, Utility::string_to_enum<Order>(getParam<std::string>("through_thickness_order")));
  _t_points = t_qrule->get_points();
  _elasticity_tensor.resize(_t_points.size());
  _stress.resize(_t_points.size());
  _stress_old.resize(_t_points.size());
  _strain_increment.resize(_t_points.size());
  _element_transformation_matrix.resize(_t_points.size());
  _covariant_transformation_matrix.resize(_t_points.size());
  _global_stress.resize(_t_points.size());
  _local_shell_stress.resize(_t_points.size());

  for (unsigned int t = 0; t < _t_points.size(); ++t)
  {
    _elasticity_tensor[t] =
        &getADMaterialProperty<RankFourTensor>("elasticity_tensor_t_points_" + std::to_string(t));
    _stress[t] = &declareADProperty<RankTwoTensor>("stress_t_points_" + std::to_string(t));
    _stress_old[t] =
        &getMaterialPropertyOldByName<RankTwoTensor>("stress_t_points_" + std::to_string(t));
    _strain_increment[t] =
        &getADMaterialProperty<RankTwoTensor>("strain_increment_t_points_" + std::to_string(t));
    // rotation matrix and stress for output purposes only
    _covariant_transformation_matrix[t] = &getMaterialProperty<RankTwoTensor>(
        "covariant_transformation_t_points_" + std::to_string(t));
    _element_transformation_matrix[t] =
        &getMaterialProperty<RankTwoTensor>("element_transformation_t_points_" + std::to_string(t));
    _global_stress[t] =
        &declareProperty<RankTwoTensor>("global_stress_t_points_" + std::to_string(t));
    _local_shell_stress[t] =
        &declareProperty<RankTwoTensor>("local_stress_t_points_" + std::to_string(t));
  }
  _local_moment_x = &declareADProperty<Real>("local_moment_x");
  _shell_force_1 = &declareADProperty<Real>("shell_force_1");
  _shell_force_2 = &declareADProperty<Real>("shell_force_2");
  _shell_shear_12 = &declareADProperty<Real>("shell_shear_12");
  _shell_shear_13 = &declareADProperty<Real>("shell_shear_13");
  _shell_shear_23 = &declareADProperty<Real>("shell_shear_23");
  _shell_moment_11 = &declareADProperty<Real>("shell_moment_11");
  _shell_moment_22 = &declareADProperty<Real>("shell_moment_22");
  _shell_moment_12 = &declareADProperty<Real>("shell_moment_12");
}

void
ADComputeTriangularShellStress::initQpStatefulProperties()
{
  // initialize stress tensor to zero
  for (unsigned int i = 0; i < _t_points.size(); ++i)
    (*_stress[i])[_qp].zero();

  (*_local_moment_x)[_qp] = 0.0;
  (*_shell_force_1)[_qp] = 0.0;
  (*_shell_force_2)[_qp] = 0.0;
  (*_shell_shear_12)[_qp] = 0.0;
  (*_shell_shear_13)[_qp] = 0.0;
  (*_shell_shear_23)[_qp] = 0.0;
  (*_shell_moment_11)[_qp] = 0.0;
  (*_shell_moment_22)[_qp] = 0.0;
  (*_shell_moment_12)[_qp] = 0.0;
}

void
ADComputeTriangularShellStress::computeQpProperties()
{
  for (unsigned int i = 0; i < _t_points.size(); ++i)
  {
    (*_stress[i])[_qp] =
        (*_stress_old[i])[_qp] + (*_elasticity_tensor[i])[_qp] * (*_strain_increment[i])[_qp];

    for (unsigned int ii = 0; ii < 3; ++ii)
      for (unsigned int jj = 0; jj < 3; ++jj)
        _unrotated_stress(ii, jj) = MetaPhysicL::raw_value((*_stress[i])[_qp](ii, jj));
    (*_global_stress[i])[_qp] = (*_covariant_transformation_matrix[i])[_qp].transpose() *
                                _unrotated_stress * (*_covariant_transformation_matrix[i])[_qp];
    (*_local_shell_stress[i])[_qp] = (*_element_transformation_matrix[i])[_qp] *
                                     (*_global_stress[i])[_qp] *
                                     (*_element_transformation_matrix[i])[_qp].transpose();
  }

  Real p11 = MetaPhysicL::raw_value((*_local_shell_stress[0])[_qp](0, 0));
  Real s11 = MetaPhysicL::raw_value((*_local_shell_stress[1])[_qp](0, 0));
  Real p22 = MetaPhysicL::raw_value((*_local_shell_stress[0])[_qp](1, 1));
  Real s22 = MetaPhysicL::raw_value((*_local_shell_stress[1])[_qp](1, 1));
  Real p12 = MetaPhysicL::raw_value((*_local_shell_stress[0])[_qp](0, 1));
  Real s12 = MetaPhysicL::raw_value((*_local_shell_stress[1])[_qp](0, 1));
  Real p13 = MetaPhysicL::raw_value((*_local_shell_stress[0])[_qp](0, 2));
  Real s13 = MetaPhysicL::raw_value((*_local_shell_stress[1])[_qp](0, 2));
  Real p23 = MetaPhysicL::raw_value((*_local_shell_stress[0])[_qp](1, 2));
  Real s23 = MetaPhysicL::raw_value((*_local_shell_stress[1])[_qp](1, 2));
  //(*_local_moment_x)[_qp]= std::pow(std::pow((p1+p2)*0.21/2,2),0.5);

  //(*_local_moment_x)[_qp]= (p1*-0.57735+p2*0.57735)*0.21*0.21/4;

  (*_shell_force_1)[_qp] = (p11 + s11) * (_thickness[_qp] / 2);
  (*_shell_force_2)[_qp] = (p22 + s22) * (_thickness[_qp] / 2);
  (*_shell_shear_12)[_qp] = (p12 + s12) * (_thickness[_qp] / 2);
  (*_shell_shear_13)[_qp] = (p13 + s13) * (_thickness[_qp] / 2);
  (*_shell_shear_23)[_qp] = (p23 + s23) * (_thickness[_qp] / 2);
  (*_shell_moment_11)[_qp] =
      -(p11 * -0.57735 + s11 * 0.57735) * (_thickness[_qp] / 2) * (_thickness[_qp] / 2);
  (*_shell_moment_22)[_qp] =
      -(p22 * -0.57735 + s22 * 0.57735) * (_thickness[_qp] / 2) * (_thickness[_qp] / 2);
  (*_shell_moment_12)[_qp] =
      -(p12 * -0.57735 + s12 * 0.57735) * (_thickness[_qp] / 2) * (_thickness[_qp] / 2);
}
