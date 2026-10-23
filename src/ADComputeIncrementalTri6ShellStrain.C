//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeIncrementalTri6ShellStrain.h"
#include "MooseMesh.h"
#include "Assembly.h"
#include "NonlinearSystem.h"
#include "MooseVariable.h"
#include "ArbitraryQuadrature.h"
#include "DenseMatrix.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/fe.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include <vector>

#include "libmesh/quadrature.h"
#include "libmesh/utility.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_type.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature_gauss.h"

registerMooseObject(MOOSEAPPNAME, ADComputeIncrementalTri6ShellStrain);

InputParameters
ADComputeIncrementalTri6ShellStrain::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute a small strain increment for a six-node triangular shell.");
  params.addRequiredCoupledVar(
      "rotations", "The rotations appropriate for the simulation geometry and coordinate system");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addRequiredCoupledVar(
      "thickness",
      "Thickness of the shell. Can be supplied as either a number or a variable name.");
  params.addRequiredParam<std::string>("through_thickness_order",
                                       "Quadrature order in out-of-plane direction");
  params.addParam<bool>(
      "large_strain", false, "Set to true to turn on finite-strain calculations.");
  params.addParam<RealVectorValue>("reference_first_local_direction",
                                   RealVectorValue(1, 0, 0),
                                   "Vector that is projected onto an element to define the "
                                   "direction of the first local coordinate axis e1.");
  return params;
}

ADComputeIncrementalTri6ShellStrain::ADComputeIncrementalTri6ShellStrain(
    const InputParameters & parameters)
  : Material(parameters),
    _nrot(coupledComponents("rotations")),
    _ndisp(coupledComponents("displacements")),
    _rot_num(_nrot),
    _disp_num(_ndisp),
    _thickness(coupledValue("thickness")),
    _large_strain(getParam<bool>("large_strain")),
    _strain_increment(),
    _total_strain(),
    _total_strain_old(),
    _nonlinear_sys(_fe_problem.getNonlinearSystemBase(/*nl_sys_num=*/0)),
    _soln_disp_index(6),
    _soln_rot_index(6),
    _soln_vector(30, 1),
    _strain_vector(5, 1),
    _nodes(6),
    _node_normal(declareADProperty<RealVectorValue>("node_normal")),
    _node_normal_old(getMaterialPropertyOldByName<RealVectorValue>("node_normal")),
    _dxyz_dxi(),
    _dxyz_deta(),
    _dxyz_dzeta(),
    _dxyz_dxi_old(),
    _dxyz_deta_old(),
    _dxyz_dzeta_old(),
    _v1(6),
    _v2(6),
    _B(),
    _B_old(),
    _ge(),
    _ge_old(),
    _J_map(),
    _J_map_old(),
    _local_transformation_matrix(),
    _local_transformation_matrix_old(),
    _covariant_transformation_matrix(),
    _covariant_transformation_matrix_old(),
    _contravariant_transformation_matrix(),
    _contravariant_transformation_matrix_old(),
    _total_global_strain(),
    _sol(_nonlinear_sys.currentSolution()),
    _sol_old(_nonlinear_sys.solutionOld()),
    _ref_first_local_dir(getParam<RealVectorValue>("reference_first_local_direction"))
{
  if (_ndisp != 3 || _nrot != 2)
    mooseError("ADComputeIncrementalTri6ShellStrain: The number of variables supplied in "
               "'displacements' must be 3 and that in 'rotations' must be 2.");

  if (!_mesh.hasSecondOrderElements())
    mooseError("ADComputeIncrementalTri6ShellStrain: The six-node triangular shell "
               "requires second-order elements.");

  if (MooseUtils::absoluteFuzzyEqual(_ref_first_local_dir.norm(), 0.0, 1e-3))
    mooseError("The norm of the defined first local axis is less than the acceptable tolerance "
               "(1e-3).");

  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    MooseVariable * disp_variable = getVar("displacements", i);
    _disp_num[i] = disp_variable->number();

    if (i < _nrot)
    {
      MooseVariable * rot_variable = getVar("rotations", i);
      _rot_num[i] = rot_variable->number();
    }
  }

  _t_qrule = std::make_unique<libMesh::QGauss>(
      1, Utility::string_to_enum<Order>(getParam<std::string>("through_thickness_order")));
  _t_points = _t_qrule->get_points();

  _strain_increment.resize(_t_points.size());
  _total_strain.resize(_t_points.size());
  _total_strain_old.resize(_t_points.size());
  _B.resize(_t_points.size());
  _B_old.resize(_t_points.size());
  _ge.resize(_t_points.size());
  _ge_old.resize(_t_points.size());
  _J_map.resize(_t_points.size());
  _J_map_old.resize(_t_points.size());
  _dxyz_dxi.resize(_t_points.size());
  _dxyz_deta.resize(_t_points.size());
  _dxyz_dzeta.resize(_t_points.size());
  _dxyz_dxi_old.resize(_t_points.size());
  _dxyz_deta_old.resize(_t_points.size());
  _dxyz_dzeta_old.resize(_t_points.size());
  _local_transformation_matrix.resize(_t_points.size());
  _local_transformation_matrix_old.resize(_t_points.size());
  _covariant_transformation_matrix.resize(_t_points.size());
  _covariant_transformation_matrix_old.resize(_t_points.size());
  _contravariant_transformation_matrix.resize(_t_points.size());
  _contravariant_transformation_matrix_old.resize(_t_points.size());
  _total_global_strain.resize(_t_points.size());

  _transformation_matrix = &declareADProperty<RankTwoTensor>("transformation_matrix_element");

  for (unsigned int i = 0; i < _t_points.size(); ++i)
  {
    _strain_increment[i] =
        &declareADProperty<RankTwoTensor>("strain_increment_t_points_" + std::to_string(i));
    _total_strain[i] =
        &declareADProperty<RankTwoTensor>("total_strain_t_points_" + std::to_string(i));
    _total_strain_old[i] =
        &getMaterialPropertyOldByName<RankTwoTensor>("total_strain_t_points_" + std::to_string(i));
    _B[i] = &declareADProperty<DenseMatrix<Real>>("B_t_points_" + std::to_string(i));
    _B_old[i] = &getMaterialPropertyOldByName<DenseMatrix<Real>>("B_t_points_" + std::to_string(i));
    _ge[i] = &declareADProperty<RankTwoTensor>("ge_t_points_" + std::to_string(i));
    _ge_old[i] = &getMaterialPropertyOldByName<RankTwoTensor>("ge_t_points_" + std::to_string(i));
    _J_map[i] = &declareADProperty<Real>("J_mapping_t_points_" + std::to_string(i));
    _J_map_old[i] = &getMaterialPropertyOldByName<Real>("J_mapping_t_points_" + std::to_string(i));
    _dxyz_dxi[i] = &declareADProperty<RealVectorValue>("dxyz_dxi_t_points_" + std::to_string(i));
    _dxyz_dxi_old[i] =
        &getMaterialPropertyOldByName<RealVectorValue>("dxyz_dxi_t_points_" + std::to_string(i));
    _dxyz_deta[i] = &declareADProperty<RealVectorValue>("dxyz_deta_t_points_" + std::to_string(i));
    _dxyz_deta_old[i] =
        &getMaterialPropertyOldByName<RealVectorValue>("dxyz_deta_t_points_" + std::to_string(i));
    _dxyz_dzeta[i] =
        &declareADProperty<RealVectorValue>("dxyz_dzeta_t_points_" + std::to_string(i));
    _dxyz_dzeta_old[i] =
        &getMaterialPropertyOldByName<RealVectorValue>("dxyz_dzeta_t_points_" + std::to_string(i));

    _local_transformation_matrix[i] =
        &declareProperty<RankTwoTensor>("local_transformation_t_points_" + std::to_string(i));
    _local_transformation_matrix_old[i] = &getMaterialPropertyOldByName<RankTwoTensor>(
        "local_transformation_t_points_" + std::to_string(i));
    _covariant_transformation_matrix[i] =
        &declareProperty<RankTwoTensor>("covariant_transformation_t_points_" + std::to_string(i));
    _covariant_transformation_matrix_old[i] = &getMaterialPropertyOldByName<RankTwoTensor>(
        "covariant_transformation_t_points_" + std::to_string(i));
    _contravariant_transformation_matrix[i] = &declareProperty<RankTwoTensor>(
        "contravariant_transformation_t_points_" + std::to_string(i));
    _contravariant_transformation_matrix_old[i] = &getMaterialPropertyOldByName<RankTwoTensor>(
        "contravariant_transformation_t_points_" + std::to_string(i));
    _total_global_strain[i] =
        &declareProperty<RankTwoTensor>("total_global_strain_t_points_" + std::to_string(i));
  }

  _e1(0) = 1;
  _e2(1) = 1;
  _e3(2) = 1;
}

void
ADComputeIncrementalTri6ShellStrain::initQpStatefulProperties()
{
  const unsigned int dim = _current_elem->dim();

  if (dim != 2)
    mooseError("ADComputeIncrementalTri6ShellStrain: Shell element is implemented only for "
               "2D elements.");

  if (_current_elem->n_nodes() != 6)
    mooseError("ADComputeIncrementalTri6ShellStrain: Shell element needs to have exactly "
               "six nodes.");

  computeGMatrix();
  computeBMatrix();
}

void
ADComputeIncrementalTri6ShellStrain::computeProperties()
{
  _2d_points = _qrule->get_points();

  const unsigned int dim = _current_elem->dim();
  FEType fe_type(Utility::string_to_enum<Order>("SECOND"),
                 Utility::string_to_enum<FEFamily>("LAGRANGE"));
  auto & fe = _fe_problem.assembly(_tid, _nonlinear_sys.number()).getFE(fe_type, dim);
  _dphidxi_map = fe->get_fe_map().get_dphidxi_map();
  _dphideta_map = fe->get_fe_map().get_dphideta_map();
  _phi_map = fe->get_fe_map().get_phi_map();

  for (unsigned int i = 0; i < _nodes.size(); ++i)
    _nodes[i] = _current_elem->node_ptr(i);

  for (unsigned int i = 0; i < _2d_points.size(); ++i)
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      (*_ge[j])[i] = (*_ge_old[j])[i];
      (*_J_map[j])[i] = (*_J_map_old[j])[i];
      (*_dxyz_dxi[j])[i] = (*_dxyz_dxi_old[j])[i];
      (*_dxyz_deta[j])[i] = (*_dxyz_deta_old[j])[i];
      (*_dxyz_dzeta[j])[i] = (*_dxyz_dzeta_old[j])[i];
      (*_B[j])[i] = (*_B_old[j])[i];
      (*_local_transformation_matrix[j])[i] = (*_local_transformation_matrix_old[j])[i];
      (*_covariant_transformation_matrix[j])[i] = (*_covariant_transformation_matrix_old[j])[i];
      (*_contravariant_transformation_matrix[j])[i] =
          (*_contravariant_transformation_matrix_old[j])[i];
    }

  computeSolnVector();
  computeNodeNormal();

  for (unsigned int i = 0; i < _2d_points.size(); ++i)
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      for (unsigned int row = 0; row < 5; ++row)
      {
        _strain_vector(row) = 0.0;
        for (unsigned int column = 0; column < 30; ++column)
          _strain_vector(row) += (*_B[j])[i](row, column) * _soln_vector(column);
      }

      (*_strain_increment[j])[i](0, 0) = _strain_vector(0);
      (*_strain_increment[j])[i](1, 1) = _strain_vector(1);
      (*_strain_increment[j])[i](0, 1) = _strain_vector(2);
      (*_strain_increment[j])[i](0, 2) = _strain_vector(3);
      (*_strain_increment[j])[i](1, 2) = _strain_vector(4);
      (*_strain_increment[j])[i](1, 0) = (*_strain_increment[j])[i](0, 1);
      (*_strain_increment[j])[i](2, 0) = (*_strain_increment[j])[i](0, 2);
      (*_strain_increment[j])[i](2, 1) = (*_strain_increment[j])[i](1, 2);

      (*_total_strain[j])[i] = (*_total_strain_old[j])[i] + (*_strain_increment[j])[i];

      for (unsigned int ii = 0; ii < Moose::dim; ++ii)
        for (unsigned int jj = 0; jj < Moose::dim; ++jj)
          _unrotated_total_strain(ii, jj) = MetaPhysicL::raw_value((*_total_strain[j])[i](ii, jj));

      (*_total_global_strain[j])[i] = (*_contravariant_transformation_matrix[j])[i] *
                                      _unrotated_total_strain *
                                      (*_contravariant_transformation_matrix[j])[i].transpose();
    }
}

// Evaluate the unnormalized quadratic-surface normal at a TRI6 node.
// Use libMesh's SECOND-order LAGRANGE shape-function derivatives rather
// than hard-coded TRI6 derivative coefficients.
namespace
{
RealVectorValue
tri6NormalAtNode(const libMesh::Elem & elem, const unsigned int k)
{
  if (elem.dim() != 2 || elem.n_nodes() != 6 || k >= 6)
    mooseError("tri6NormalAtNode requires a TRI6 element and local node index 0..5");

  const libMesh::Point p = elem.master_point(k);

  RealVectorValue g1(0., 0., 0.);
  RealVectorValue g2(0., 0., 0.);

  for (unsigned int a = 0; a < 6; ++a)
  {
    const Real dN_dxi =
        libMesh::FE<2, libMesh::LAGRANGE>::shape_deriv(&elem, libMesh::SECOND, a, 0, p);
    const Real dN_deta =
        libMesh::FE<2, libMesh::LAGRANGE>::shape_deriv(&elem, libMesh::SECOND, a, 1, p);

    const auto & X = elem.point(a);
    for (unsigned int c = 0; c < 3; ++c)
    {
      g1(c) += dN_dxi * X(c);
      g2(c) += dN_deta * X(c);
    }
  }

  // Magnitude provides the local surface-Jacobian weighting.
  return g1.cross(g2);
}
} // namespace

void
ADComputeIncrementalTri6ShellStrain::computeGMatrix()
{
  using std::sqrt;

  _2d_points = _qrule->get_points();

  const unsigned int dim = _current_elem->dim();
  FEType fe_type(Utility::string_to_enum<Order>("SECOND"),
                 Utility::string_to_enum<FEFamily>("LAGRANGE"));
  auto & fe = _fe_problem.assembly(_tid, _nonlinear_sys.number()).getFE(fe_type, dim);
  _dphidxi_map = fe->get_fe_map().get_dphidxi_map();
  _dphideta_map = fe->get_fe_map().get_dphideta_map();
  _phi_map = fe->get_fe_map().get_phi_map();

  for (unsigned int i = 0; i < _nodes.size(); ++i)
    _nodes[i] = _current_elem->node_ptr(i);

  // Deliberately build adjacency and recompute nodal normals for every element.
  // This avoids a material-level cache, at the expense of repeated work.
  // The adjacency must contain ALL elements touching each node (MPI ghosting).
  const auto & mesh = _mesh.getMesh();
  std::vector<std::vector<dof_id_type>> nodes_to_elem_map;
  libMesh::MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

  for (unsigned int k = 0; k < _nodes.size(); ++k)
  {
    const dof_id_type node_id = _nodes[k]->id();
    if (node_id >= nodes_to_elem_map.size())
      mooseError("No node-to-element adjacency for TRI6 node ", node_id);

    RealVectorValue normal_sum(0., 0., 0.);
    unsigned int contributions = 0;

    for (const dof_id_type elem_id : nodes_to_elem_map[node_id])
    {
      const libMesh::Elem * neighbor = mesh.query_elem_ptr(elem_id);
      if (!neighbor)
        mooseError("Missing adjacent element ", elem_id, " at TRI6 node ", node_id,
                   ". Ensure complete MPI ghosting of the nodal patch.");

      // Do not mix normals of distinct shell blocks at a shared node.
      if (!neighbor->active() ||
          neighbor->subdomain_id() != _current_elem->subdomain_id())
        continue;

      if (neighbor->dim() != 2 || neighbor->n_nodes() != 6)
        mooseError("Expected TRI6 neighbor at node ", node_id,
                   ", but found element ", elem_id);

      unsigned int local_node = 0;
      for (; local_node < 6; ++local_node)
        if (neighbor->node_id(local_node) == node_id)
          break;
      if (local_node == 6)
        mooseError("Invalid node-to-element adjacency at node ", node_id);

      const RealVectorValue contribution = tri6NormalAtNode(*neighbor, local_node);
      if (contribution.norm() < 1e-14)
        mooseError("Degenerate TRI6 normal at node ", node_id,
                   " in element ", elem_id);

      normal_sum += contribution;
      ++contributions;
    }

    if (!contributions || normal_sum.norm() < 1e-14)
      mooseError("Cannot construct averaged TRI6 normal at node ", node_id,
                 ". Check element orientation and nodal adjacency.");

    _node_normal[k] = normal_sum / normal_sum.norm();
  }

  ADRankTwoTensor zero_ad_tensor;
  ADDenseMatrix zero_b_matrix(5, 30);
  ADRealVectorValue zero_ad_vector;
  RankTwoTensor zero_tensor;

  for (unsigned int t = 0; t < _t_points.size(); ++t)
  {
    (*_strain_increment[t])[_qp] = zero_ad_tensor;
    (*_total_strain[t])[_qp] = zero_ad_tensor;
    (*_B[t])[_qp] = zero_b_matrix;
    (*_ge[t])[_qp] = zero_ad_tensor;
    (*_J_map[t])[_qp] = 0;
    (*_dxyz_dxi[t])[_qp] = zero_ad_vector;
    (*_dxyz_deta[t])[_qp] = zero_ad_vector;
    (*_dxyz_dzeta[t])[_qp] = zero_ad_vector;
    (*_local_transformation_matrix[t])[_qp] = zero_tensor;
    (*_covariant_transformation_matrix[t])[_qp] = zero_tensor;
    (*_contravariant_transformation_matrix[t])[_qp] = zero_tensor;
  }

  for (unsigned int i = 0; i < _2d_points.size(); ++i)
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      (*_dxyz_dxi[j])[i].zero();
      (*_dxyz_deta[j])[i].zero();
      (*_dxyz_dzeta[j])[i].zero();

      for (unsigned int component = 0; component < Moose::dim; ++component)
        for (unsigned int k = 0; k < _nodes.size(); ++k)
        {
          (*_dxyz_dxi[j])[i](component) += _dphidxi_map[k][i] * (*_nodes[k])(component) +
                                           _t_points[j](0) / 2.0 * _thickness[i] *
                                               _dphidxi_map[k][i] * _node_normal[k](component);
          (*_dxyz_deta[j])[i](component) += _dphideta_map[k][i] * (*_nodes[k])(component) +
                                            _t_points[j](0) / 2.0 * _thickness[i] *
                                                _dphideta_map[k][i] * _node_normal[k](component);
          (*_dxyz_dzeta[j])[i](component) +=
              _thickness[i] * _phi_map[k][i] * _node_normal[k](component) / 2.0;
        }
    }

  for (unsigned int i = 0; i < _2d_points.size(); ++i)
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      ADRankTwoTensor gmn;
      RankTwoTensor J;

      for (unsigned int component = 0; component < Moose::dim; ++component)
      {
        gmn(0, 0) += (*_dxyz_dxi[j])[i](component) * (*_dxyz_dxi[j])[i](component);
        gmn(1, 1) += (*_dxyz_deta[j])[i](component) * (*_dxyz_deta[j])[i](component);
        gmn(2, 2) += (*_dxyz_dzeta[j])[i](component) * (*_dxyz_dzeta[j])[i](component);
        gmn(0, 1) += (*_dxyz_dxi[j])[i](component) * (*_dxyz_deta[j])[i](component);
        gmn(0, 2) += (*_dxyz_dxi[j])[i](component) * (*_dxyz_dzeta[j])[i](component);
        gmn(1, 2) += (*_dxyz_deta[j])[i](component) * (*_dxyz_dzeta[j])[i](component);

        J(0, component) = MetaPhysicL::raw_value((*_dxyz_dxi[j])[i](component));
        J(1, component) = MetaPhysicL::raw_value((*_dxyz_deta[j])[i](component));
        J(2, component) = MetaPhysicL::raw_value((*_dxyz_dzeta[j])[i](component));
      }

      gmn(1, 0) = gmn(0, 1);
      gmn(2, 0) = gmn(0, 2);
      gmn(2, 1) = gmn(1, 2);

      (*_J_map[j])[i] = sqrt(gmn.det());
      (*_covariant_transformation_matrix[j])[i] = J;
      (*_contravariant_transformation_matrix[j])[i] = J.inverse();

      const Real norm_x = sqrt(J(0, 0) * J(0, 0) + J(0, 1) * J(0, 1) + J(0, 2) * J(0, 2));
      const Real norm_y = sqrt(J(1, 0) * J(1, 0) + J(1, 1) * J(1, 1) + J(1, 2) * J(1, 2));
      const Real norm_z = sqrt(J(2, 0) * J(2, 0) + J(2, 1) * J(2, 1) + J(2, 2) * J(2, 2));

      J(0, 0) /= norm_x;
      J(0, 1) /= norm_x;
      J(0, 2) /= norm_x;

      J(1, 0) /= norm_y;
      J(1, 1) /= norm_y;
      J(1, 2) /= norm_y;

      J(2, 0) /= norm_z;
      J(2, 1) /= norm_z;
      J(2, 2) /= norm_z;

      computeLocalCoordinates((*_dxyz_dzeta[j])[i]);

      const auto local_rotation_mat = ADRankTwoTensor::initializeFromRows(_e1, _e2, _e3);
      for (const auto ii : make_range(Moose::dim))
        for (const auto jj : make_range(Moose::dim))
          (*_local_transformation_matrix[j])[i](ii, jj) =
              MetaPhysicL::raw_value(local_rotation_mat(ii, jj));

      const auto denom = (*_dxyz_dxi[j])[i] * (*_dxyz_deta[j])[i].cross((*_dxyz_dzeta[j])[i]);
      const auto gi0 = (*_dxyz_deta[j])[i].cross((*_dxyz_dzeta[j])[i]) / denom;
      const auto gi1 = (*_dxyz_dzeta[j])[i].cross((*_dxyz_dxi[j])[i]) / denom;
      const auto gi2 = (*_dxyz_dxi[j])[i].cross((*_dxyz_deta[j])[i]) / denom;

      (*_ge[j])[i](0, 0) = gi0 * _e1;
      (*_ge[j])[i](0, 1) = gi0 * _e2;
      (*_ge[j])[i](0, 2) = gi0 * _e3;
      (*_ge[j])[i](1, 0) = gi1 * _e1;
      (*_ge[j])[i](1, 1) = gi1 * _e2;
      (*_ge[j])[i](1, 2) = gi1 * _e3;
      (*_ge[j])[i](2, 0) = gi2 * _e1;
      (*_ge[j])[i](2, 1) = gi2 * _e2;
      (*_ge[j])[i](2, 2) = gi2 * _e3;
    }
}

void
ADComputeIncrementalTri6ShellStrain::computeLocalCoordinates(
    const ADRealVectorValue & g3)
{
  using std::abs;

  // Unit director in the thickness direction at this integration point.
  if (MetaPhysicL::raw_value(g3.norm()) < 1e-14)
    mooseError("Zero thickness-direction vector in shell element ", _current_elem->id());
  _e3 = g3 / g3.norm();

  _e1 = _ref_first_local_dir;
  _e1 /= _e1.norm();

  if (MooseUtils::absoluteFuzzyEqual(abs(_e1 * _e3), 1.0, 0.05 * libMesh::pi / 180.0))
    mooseError("The reference first local direction must not be parallel to the shell thickness direction.");

  _e2 = _e3.cross(_e1);
  _e2 /= _e2.norm();

  _e1 = _e2.cross(_e3);
  _e1 /= _e1.norm();
}

void
ADComputeIncrementalTri6ShellStrain::computeNodeNormal()
{
  for (unsigned int k = 0; k < _nodes.size(); ++k)
    _node_normal[k] = _node_normal_old[k];
}

void
ADComputeIncrementalTri6ShellStrain::computeBMatrix()
{

  for (unsigned int k = 0; k < _nodes.size(); ++k)
  {
    _v2[k] = _node_normal[k].cross(_ref_first_local_dir);

    if (MooseUtils::absoluteFuzzyEqual(_v2[k].norm(), 0.0, 1e-6))
      mooseError("The reference first local direction is parallel to a nodal shell normal.");

    _v2[k] /= _v2[k].norm();
    _v1[k] = _v2[k].cross(_node_normal[k]);
    _v1[k] /= _v1[k].norm();
  }

  for (unsigned int i = 0; i < _2d_points.size(); ++i)
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      (*_B[j])[i].resize(5, 30);
      (*_B[j])[i].zero();

      for (unsigned int k = 0; k < _nodes.size(); ++k)
      {
        (*_B[j])[i](0, k) = _dphidxi_map[k][i] * (*_dxyz_dxi[j])[i](0);
        (*_B[j])[i](0, 6 + k) = _dphidxi_map[k][i] * (*_dxyz_dxi[j])[i](1);
        (*_B[j])[i](0, 12 + k) = _dphidxi_map[k][i] * (*_dxyz_dxi[j])[i](2);
        (*_B[j])[i](0, 18 + k) = _dphidxi_map[k][i] * _t_points[j](0) / 2.0 * _thickness[i] *
                                 (-_v2[k] * (*_dxyz_dxi[j])[i]);
        (*_B[j])[i](0, 24 + k) = _dphidxi_map[k][i] * _t_points[j](0) / 2.0 * _thickness[i] *
                                 (_v1[k] * (*_dxyz_dxi[j])[i]);

        (*_B[j])[i](1, k) = _dphideta_map[k][i] * (*_dxyz_deta[j])[i](0);
        (*_B[j])[i](1, 6 + k) = _dphideta_map[k][i] * (*_dxyz_deta[j])[i](1);
        (*_B[j])[i](1, 12 + k) = _dphideta_map[k][i] * (*_dxyz_deta[j])[i](2);
        (*_B[j])[i](1, 18 + k) = _dphideta_map[k][i] * _t_points[j](0) / 2.0 * _thickness[i] *
                                 (-_v2[k] * (*_dxyz_deta[j])[i]);
        (*_B[j])[i](1, 24 + k) = _dphideta_map[k][i] * _t_points[j](0) / 2.0 * _thickness[i] *
                                 (_v1[k] * (*_dxyz_deta[j])[i]);

        (*_B[j])[i](2, k) = 0.5 * (_dphideta_map[k][i] * (*_dxyz_dxi[j])[i](0) +
                                   _dphidxi_map[k][i] * (*_dxyz_deta[j])[i](0));
        (*_B[j])[i](2, 6 + k) = 0.5 * (_dphideta_map[k][i] * (*_dxyz_dxi[j])[i](1) +
                                       _dphidxi_map[k][i] * (*_dxyz_deta[j])[i](1));
        (*_B[j])[i](2, 12 + k) = 0.5 * (_dphideta_map[k][i] * (*_dxyz_dxi[j])[i](2) +
                                        _dphidxi_map[k][i] * (*_dxyz_deta[j])[i](2));
        (*_B[j])[i](2, 18 + k) =
            0.25 * _t_points[j](0) * _thickness[i] * -_v2[k] *
            (_dphideta_map[k][i] * (*_dxyz_dxi[j])[i] + _dphidxi_map[k][i] * (*_dxyz_deta[j])[i]);
        (*_B[j])[i](2, 24 + k) =
            0.25 * _t_points[j](0) * _thickness[i] * _v1[k] *
            ((*_dxyz_deta[j])[i] * _dphidxi_map[k][i] + (*_dxyz_dxi[j])[i] * _dphideta_map[k][i]);

        (*_B[j])[i](3, k) = 0.5 * _dphidxi_map[k][i] * (*_dxyz_dzeta[j])[i](0);
        (*_B[j])[i](3, 6 + k) = 0.5 * _dphidxi_map[k][i] * (*_dxyz_dzeta[j])[i](1);
        (*_B[j])[i](3, 12 + k) = 0.5 * _dphidxi_map[k][i] * (*_dxyz_dzeta[j])[i](2);
        (*_B[j])[i](3, 18 + k) = 0.25 * _thickness[i] * -_v2[k] *
                                 (_phi_map[k][i] * (*_dxyz_dxi[j])[i] +
                                  _t_points[j](0) * _dphidxi_map[k][i] * (*_dxyz_dzeta[j])[i]);
        (*_B[j])[i](3, 24 + k) = 0.25 * _thickness[i] * _v1[k] *
                                 (_phi_map[k][i] * (*_dxyz_dxi[j])[i] +
                                  _t_points[j](0) * _dphidxi_map[k][i] * (*_dxyz_dzeta[j])[i]);

        (*_B[j])[i](4, k) = 0.5 * _dphideta_map[k][i] * (*_dxyz_dzeta[j])[i](0);
        (*_B[j])[i](4, 6 + k) = 0.5 * _dphideta_map[k][i] * (*_dxyz_dzeta[j])[i](1);
        (*_B[j])[i](4, 12 + k) = 0.5 * _dphideta_map[k][i] * (*_dxyz_dzeta[j])[i](2);
        (*_B[j])[i](4, 18 + k) = 0.25 * _thickness[i] * -_v2[k] *
                                 (_phi_map[k][i] * (*_dxyz_deta[j])[i] +
                                  _t_points[j](0) * _dphideta_map[k][i] * (*_dxyz_dzeta[j])[i]);
        (*_B[j])[i](4, 24 + k) = 0.25 * _thickness[i] * _v1[k] *
                                 (_phi_map[k][i] * (*_dxyz_deta[j])[i] +
                                  _t_points[j](0) * _dphideta_map[k][i] * (*_dxyz_dzeta[j])[i]);
      }

      updateGVectors();
    }
}

void
ADComputeIncrementalTri6ShellStrain::computeSolnVector()
{
  _soln_vector.zero();

  for (unsigned int j = 0; j < _nodes.size(); ++j)
  {
    _soln_disp_index[j].resize(_ndisp);
    _soln_rot_index[j].resize(_nrot);

    for (unsigned int i = 0; i < _ndisp; ++i)
    {
      _soln_disp_index[j][i] = _nodes[j]->dof_number(_nonlinear_sys.number(), _disp_num[i], 0);
      _soln_vector(j + i * _nodes.size()) =
          (*_sol)(_soln_disp_index[j][i]) - _sol_old(_soln_disp_index[j][i]);

      if (ADReal::do_derivatives)
        Moose::derivInsert(
            _soln_vector(j + i * _nodes.size()).derivatives(), _soln_disp_index[j][i], 1.0);
    }

    for (unsigned int i = 0; i < _nrot; ++i)
    {
      _soln_rot_index[j][i] = _nodes[j]->dof_number(_nonlinear_sys.number(), _rot_num[i], 0);
      _soln_vector(j + 18 + i * _nodes.size()) =
          (*_sol)(_soln_rot_index[j][i]) - _sol_old(_soln_rot_index[j][i]);

      if (ADReal::do_derivatives)
        Moose::derivInsert(
            _soln_vector(j + 18 + i * _nodes.size()).derivatives(), _soln_rot_index[j][i], 1.0);
    }
  }
}
