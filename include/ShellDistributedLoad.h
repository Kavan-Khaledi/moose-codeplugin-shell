//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GenericKernel.h"

class Function;

/**
 * ShellDistributedLoad computes the body force (force/volume) given the acceleration of
 * ShellDistributedLoad (value) and the density
 */
template <bool is_ad>
class ShellDistributedLoadTempl : public GenericKernel<is_ad>
{
public:
  static InputParameters validParams();

  ShellDistributedLoadTempl(const InputParameters & parameters);

protected:
  virtual GenericReal<is_ad> computeQpResidual() override;

  const Real _value;
  const Function & _function;

  // _alpha parameter for HHT time integration scheme

  usingTransientInterfaceMembers;
  using GenericKernel<is_ad>::_i;
  using GenericKernel<is_ad>::_q_point;
  using GenericKernel<is_ad>::_qp;
  using GenericKernel<is_ad>::_test;
};

using ShellDistributedLoad = ShellDistributedLoadTempl<false>;
using ADShellDistributedLoad = ShellDistributedLoadTempl<true>;
