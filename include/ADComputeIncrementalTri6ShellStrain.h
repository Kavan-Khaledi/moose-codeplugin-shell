//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DenseMatrix.h"
#include "Material.h"
#include "ADRankTwoTensorForward.h"

namespace libMesh
{
class QGauss;
}

class ADComputeIncrementalTri6ShellStrain : public Material
{
public:
  static InputParameters validParams();

  ADComputeIncrementalTri6ShellStrain(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeProperties() override;

  /// Computes the 30x1 solution vector and its derivatives for each six-node shell element
  virtual void computeSolnVector();

  /// Computes the B matrix that connects strains to nodal displacements and rotations
  virtual void computeBMatrix();

  /// Computes the node normal at each node
  virtual void computeNodeNormal();

  /// Updates the vectors required for shear-locking computation for finite rotations
  virtual void updateGVectors() {}

  /// Updates covariant vectors at each quadrature point for finite rotations
  virtual void updatedxyz() {}

  /// Computes the transformation matrices and covariant/contravariant base vectors
  virtual void computeGMatrix();

  /// Computes the element local coordinate system using the interpolated thickness direction
  virtual void computeLocalCoordinates(const ADRealVectorValue & g3);

  /// Number of coupled rotational variables
  unsigned int _nrot;

  /// Number of coupled displacement variables
  unsigned int _ndisp;

  /// Variable numbers corresponding to the rotational variables
  std::vector<unsigned int> _rot_num;

  /// Variable numbers corresponding to the displacement variables
  std::vector<unsigned int> _disp_num;

  /// Coupled variable for the shell thickness
  const VariableValue & _thickness;

  /// Flag to compute large strains
  const bool _large_strain;

  /// Strain increment in the covariant coordinate system
  std::vector<ADMaterialProperty<RankTwoTensor> *> _strain_increment;

  /// Total strain in the covariant coordinate system
  std::vector<ADMaterialProperty<RankTwoTensor> *> _total_strain;

  /// Old total strain in the covariant coordinate system
  std::vector<const MaterialProperty<RankTwoTensor> *> _total_strain_old;

  /// Reference to the nonlinear system object
  NonlinearSystemBase & _nonlinear_sys;

  /// Displacement DOF indices for three directions at six nodes
  std::vector<std::vector<unsigned int>> _soln_disp_index;

  /// Rotation DOF indices for two directions at six nodes
  std::vector<std::vector<unsigned int>> _soln_rot_index;

  /// Incremental solution vector for the 30 element DOFs
  ADDenseVector _soln_vector;

  /// Strain vector containing two axial and three shear components
  ADDenseVector _strain_vector;

  /// Pointers to the six shell-element nodes
  std::vector<const Node *> _nodes;

  /// Normal vectors at the six shell nodes
  ADMaterialProperty<RealVectorValue> & _node_normal;

  /// Old normal vectors at the six shell nodes
  const MaterialProperty<RealVectorValue> & _node_normal_old;

  /// Quadrature rule in the through-thickness direction
  std::unique_ptr<libMesh::QGauss> _t_qrule;

  /// Through-thickness quadrature points in isoparametric coordinates
  std::vector<Point> _t_points;

  /// In-plane quadrature points in isoparametric coordinates
  std::vector<Point> _2d_points;

  /// Shape-function derivatives with respect to xi
  std::vector<std::vector<Real>> _dphidxi_map;

  /// Shape-function derivatives with respect to eta
  std::vector<std::vector<Real>> _dphideta_map;

  /// Shape-function values
  std::vector<std::vector<Real>> _phi_map;

  /// Derivative of global coordinates with respect to xi
  std::vector<ADMaterialProperty<RealVectorValue> *> _dxyz_dxi;

  /// Derivative of global coordinates with respect to eta
  std::vector<ADMaterialProperty<RealVectorValue> *> _dxyz_deta;

  /// Derivative of global coordinates with respect to zeta
  std::vector<ADMaterialProperty<RealVectorValue> *> _dxyz_dzeta;

  /// Old derivative of global coordinates with respect to xi
  std::vector<const MaterialProperty<RealVectorValue> *> _dxyz_dxi_old;

  /// Old derivative of global coordinates with respect to eta
  std::vector<const MaterialProperty<RealVectorValue> *> _dxyz_deta_old;

  /// Old derivative of global coordinates with respect to zeta
  std::vector<const MaterialProperty<RealVectorValue> *> _dxyz_dzeta_old;

  /// First tangential vector at each node
  std::vector<ADRealVectorValue> _v1;

  /// Second tangential vector at each node
  std::vector<ADRealVectorValue> _v2;

  /// Small-strain B matrix
  std::vector<ADMaterialProperty<DenseMatrix<Real>> *> _B;

  /// Old small-strain B matrix
  std::vector<const MaterialProperty<DenseMatrix<Real>> *> _B_old;

  /// Transformation matrix used for elasticity-tensor conversion
  std::vector<ADMaterialProperty<RankTwoTensor> *> _ge;

  /// Old transformation matrix used for elasticity-tensor conversion
  std::vector<const MaterialProperty<RankTwoTensor> *> _ge_old;

  /// Jacobian of the isoparametric mapping
  std::vector<ADMaterialProperty<Real> *> _J_map;

  /// Old Jacobian of the isoparametric mapping
  std::vector<const MaterialProperty<Real> *> _J_map_old;

  /// Transformation from global coordinates to the element local coordinate system
  std::vector<MaterialProperty<RankTwoTensor> *> _local_transformation_matrix;
  std::vector<const MaterialProperty<RankTwoTensor> *> _local_transformation_matrix_old;

  /// Matrix of covariant base vectors
  std::vector<MaterialProperty<RankTwoTensor> *> _covariant_transformation_matrix;
  std::vector<const MaterialProperty<RankTwoTensor> *> _covariant_transformation_matrix_old;

  /// Matrix of contravariant base vectors
  std::vector<MaterialProperty<RankTwoTensor> *> _contravariant_transformation_matrix;
  std::vector<const MaterialProperty<RankTwoTensor> *> _contravariant_transformation_matrix_old;

  /// Total strain expressed in global coordinates
  std::vector<MaterialProperty<RankTwoTensor> *> _total_global_strain;

  /// Element transformation matrix material property
  ADMaterialProperty<RankTwoTensor> * _transformation_matrix;

  /// Current and old nonlinear-system solution vectors
  const NumericVector<Number> * const & _sol;
  const NumericVector<Number> & _sol_old;

  /// Reference direction used to define the first local in-plane axis
  const RealVectorValue _ref_first_local_dir;

  /// Temporary tensor used when writing the total strain in global coordinates
  RankTwoTensor _unrotated_total_strain;

  /// Element local coordinate axes
  ADRealVectorValue _e1;
  ADRealVectorValue _e2;
  ADRealVectorValue _e3;
};
