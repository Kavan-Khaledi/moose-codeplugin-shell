//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ShellDistributedLoad.h"
#include "Function.h"

registerMooseObject(MOOSEAPPNAME, ShellDistributedLoad);
registerMooseObject(MOOSEAPPNAME, ADShellDistributedLoad);

template <bool is_ad>
InputParameters
ShellDistributedLoadTempl<is_ad>::validParams()
{
  InputParameters params = GenericKernel<is_ad>::validParams();
  params.addClassDescription("Apply ShellDistributedLoad. Value is in units of MPa/m2");
  params.addParam<bool>("use_displaced_mesh", true, "Displaced mesh defaults to true");
  params.addRequiredParam<Real>(
      "value", "Value multiplied against the residual");
  params.addParam<FunctionName>(
      "function", "1", "A function that describes the gravitational force");
  return params;
}

template <bool is_ad>
ShellDistributedLoadTempl<is_ad>::ShellDistributedLoadTempl(const InputParameters & parameters)
  : GenericKernel<is_ad>(parameters),
    _value(this->template getParam<Real>("value")),
    _function(this->getFunction("function"))
{
}

template <bool is_ad>
GenericReal<is_ad>
ShellDistributedLoadTempl<is_ad>::computeQpResidual()
{
  Real factor = _value * _function.value(_t , _q_point[_qp]);
  return _test[_i][_qp] * factor;
}

template class ShellDistributedLoadTempl<false>;
template class ShellDistributedLoadTempl<true>;
