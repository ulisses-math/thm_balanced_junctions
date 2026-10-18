//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Adapted from ADVolumeJunction1PhaseBC of the MOOSE thermal_hydraulics module (Battelle Energy
//* Alliance, LLC, LGPL 2.1).

#include "ADStagnationJunction1PhaseBC.h"
#include "ADVolumeJunctionBaseUserObject.h"
#include "THMIndicesVACE.h"

registerMooseObject("ThmBalancedJunctionsApp", ADStagnationJunction1PhaseBC);

InputParameters
ADStagnationJunction1PhaseBC::validParams()
{
  InputParameters params = ADOneDIntegratedBC::validParams();

  params.addRequiredParam<unsigned int>("connection_index", "Index of the connected flow channel");
  params.addRequiredParam<UserObjectName>("stagnation_junction_uo",
                                          "Stagnation junction user object name");

  params.addRequiredCoupledVar("A_elem", "Cross-sectional area, elemental");
  params.addRequiredCoupledVar("A_linear", "Cross-sectional area, linear");

  params.addRequiredCoupledVar("rhoA", "Flow channel variable: rho*A");
  params.addRequiredCoupledVar("rhouA", "Flow channel variable: rho*u*A");
  params.addRequiredCoupledVar("rhoEA", "Flow channel variable: rho*E*A");

  params.addClassDescription(
      "Adds boundary fluxes for flow channels connected to a 1-phase stagnation junction");

  return params;
}

ADStagnationJunction1PhaseBC::ADStagnationJunction1PhaseBC(const InputParameters & params)
  : ADOneDIntegratedBC(params),

    _connection_index(getParam<unsigned int>("connection_index")),
    _junction_uo(getUserObject<ADVolumeJunctionBaseUserObject>("stagnation_junction_uo")),

    _A_elem(adCoupledValue("A_elem")),
    _A_linear(adCoupledValue("A_linear")),

    _rhoA_jvar(coupled("rhoA")),
    _rhouA_jvar(coupled("rhouA")),
    _rhoEA_jvar(coupled("rhoEA")),

    _flow_channel_jvar_map(getFlowChannelIndexMapping()),
    _equation_index(_flow_channel_jvar_map.at(_var.number()))
{
}

ADReal
ADStagnationJunction1PhaseBC::computeQpResidual()
{
  const auto & flux = _junction_uo.getFlux(_connection_index);

  // the flux function uses A_elem, the boundary uses A_linear
  return flux[_equation_index] * _A_linear[_qp] / _A_elem[_qp] * _normal * _test[_i][_qp];
}

std::map<unsigned int, unsigned int>
ADStagnationJunction1PhaseBC::getFlowChannelIndexMapping() const
{
  std::map<unsigned int, unsigned int> jvar_map;
  jvar_map.insert(std::pair<unsigned int, unsigned int>(_rhoA_jvar, THMVACE1D::MASS));
  jvar_map.insert(std::pair<unsigned int, unsigned int>(_rhouA_jvar, THMVACE1D::MOMENTUM));
  jvar_map.insert(std::pair<unsigned int, unsigned int>(_rhoEA_jvar, THMVACE1D::ENERGY));

  return jvar_map;
}
