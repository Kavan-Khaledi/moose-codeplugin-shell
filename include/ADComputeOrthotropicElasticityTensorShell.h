//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

#define usingComputeIsotropicElasticityTensorShellMembers usingMaterialMembers

namespace libMesh
{
class QGauss;
}

class ADComputeOrthotropicElasticityTensorShell : public Material
{
public:
  static InputParameters validParams();

  ADComputeOrthotropicElasticityTensorShell(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Lame' paramters
  Real _poissons_ratio_12;
  Real _shear_modulus_12;
  Real _shear_modulus_13;
  Real _shear_modulus_23;
  Real _youngs_modulus_1;
  Real _youngs_modulus_2;

  /// Individual elasticity tensor
  RankFourTensor _Cijkl;

  /// Quadrature points along thickness
  std::vector<Point> _t_points;

  /// Material property elasticity tensor
  std::vector<ADMaterialProperty<RankFourTensor> *> _elasticity_tensor;

  /// Prefactor function used to modify (i.e., multiply) the material stiffness
  const Function * const _prefactor_function;

  /// Material property for ge matrix
  std::vector<const ADMaterialProperty<RankTwoTensor> *> _ge;
};
