//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Parts of this file are adapted from ADVolumeJunction1PhaseUserObject of the MOOSE
//* thermal_hydraulics module (Battelle Energy Alliance, LLC, LGPL 2.1).

#include "ADStagnationJunction1PhaseUserObject.h"
#include "SinglePhaseFluidProperties.h"
#include "THMIndicesVACE.h"
#include "HydrostaticFace.h"
#include "StagnationJunction1Phase.h"
#include "ADNumericalFlux3EqnBase.h"
#include "Numerics.h"
#include "THMUtils.h"
#include "metaphysicl/parallel_numberarray.h"
#include "metaphysicl/parallel_dualnumber.h"
#include "metaphysicl/parallel_semidynamicsparsenumberarray.h"
#include "libmesh/parallel_algebra.h"

registerMooseObject("ThmBalancedJunctionsApp", ADStagnationJunction1PhaseUserObject);

InputParameters
ADStagnationJunction1PhaseUserObject::validParams()
{
  InputParameters params = ADVolumeJunctionBaseUserObject::validParams();

  params.addRequiredCoupledVar("A", "Cross-sectional area of connected flow channels");
  params.addRequiredCoupledVar("rhoA", "rho*A of the connected flow channels");
  params.addRequiredCoupledVar("rhouA", "rhou*A of the connected flow channels");
  params.addRequiredCoupledVar("rhoEA", "rhoE*A of the connected flow channels");

  params.addRequiredCoupledVar("rhoV", "rho*V of the junction");
  params.addRequiredCoupledVar("rhoEV", "rho*E*V of the junction");

  params.addRequiredParam<UserObjectName>("fp", "Fluid properties user object name");

  params.addRequiredParam<Real>("K", "Form loss coefficient [-], not used (no momentum equation)");
  params.addRequiredParam<Real>("A_ref", "Reference area [m^2], not used");

  params.addParam<RealVectorValue>(
      "gravity_vector", RealVectorValue(0, 0, -9.81), "Gravitational acceleration vector [m/s^2]");
  params.addParam<bool>("hydrostatic_correction",
                        true,
                        "Extrapolate the channel cell state hydrostatically to the junction face");

  params.addClassDescription(
      "Computes and caches flux and residual vectors for a 1-phase stagnation junction");

  params.declareControllable("K");
  return params;
}

ADStagnationJunction1PhaseUserObject::ADStagnationJunction1PhaseUserObject(
    const InputParameters & params)
  : ADVolumeJunctionBaseUserObject(params),

    _A(adCoupledValue("A")),
    _rhoA(adCoupledValue("rhoA")),
    _rhouA(adCoupledValue("rhouA")),
    _rhoEA(adCoupledValue("rhoEA")),

    _K(getParam<Real>("K")),
    _A_ref(getParam<Real>("A_ref")),

    _gravity_vector(getParam<RealVectorValue>("gravity_vector")),
    _hydrostatic_correction(getParam<bool>("hydrostatic_correction")),

    _fp(getUserObject<SinglePhaseFluidProperties>("fp"))
{
  _flow_variable_names.resize(THMVACE1D::N_FLUX_OUTPUTS);
  _flow_variable_names[THMVACE1D::RHOA] = "rhoA";
  _flow_variable_names[THMVACE1D::RHOUA] = "rhouA";
  _flow_variable_names[THMVACE1D::RHOEA] = "rhoEA";

  _scalar_variable_names.resize(StagnationJunction1Phase::N_EQ);
  _scalar_variable_names[StagnationJunction1Phase::RHOV_INDEX] = "rhoV";
  _scalar_variable_names[StagnationJunction1Phase::RHOEV_INDEX] = "rhoEV";

  _junction_var_values.resize(StagnationJunction1Phase::N_EQ);
  _junction_var_values[StagnationJunction1Phase::RHOV_INDEX] = &coupledJunctionValue("rhoV");
  _junction_var_values[StagnationJunction1Phase::RHOEV_INDEX] = &coupledJunctionValue("rhoEV");

  _numerical_flux_uo.resize(_n_connections);
  for (std::size_t i = 0; i < _n_connections; i++)
    _numerical_flux_uo[i] = &getUserObjectByName<ADNumericalFlux3EqnBase>(_numerical_flux_names[i]);

  // K is checked only here: a K changed later by a Controls block is not reported, but the junction
  // never reads it
  if (_K != 0.0)
    mooseWarning("the form loss coefficient K = ",
                 _K,
                 " is accepted but not applied by StagnationJunction1Phase, because the junction "
                 "has no momentum equation; put the loss in a connected flow channel instead.");
}

void
ADStagnationJunction1PhaseUserObject::computeFluxesAndResiduals(const unsigned int & c)
{
  const Real din = _normal[c];
  const RealVectorValue di = _dir[0];
  const RealVectorValue ni = di * din;
  const Real nJi_dot_di = -din;

  // channel state at the junction face; without the correction this is the cell-average
  // (cell-center) state, as in the official VolumeJunction1Phase
  ADReal rhoEA_f = _rhoEA[0];
  if (_hydrostatic_correction)
  {
    // x_face is the quadrature point of the boundary side, x_c the centroid of the cell
    const RealVectorValue dx = _q_point[0] - _current_elem->vertex_average();
    HydrostaticFace::extrapolate(_fp, _A[0], _gravity_vector * dx, _rhoA[0], _rhouA[0], rhoEA_f);
  }

  std::vector<ADReal> Ui(THMVACE3D::N_FLUX_INPUTS, 0.);
  Ui[THMVACE3D::RHOA] = _rhoA[0];
  Ui[THMVACE3D::RHOUA] = _rhouA[0] * di(0);
  Ui[THMVACE3D::RHOVA] = _rhouA[0] * di(1);
  Ui[THMVACE3D::RHOWA] = _rhouA[0] * di(2);
  Ui[THMVACE3D::RHOEA] = rhoEA_f;
  Ui[THMVACE3D::AREA] = _A[0];

  const auto flux_3d = compute3DFlux(*(_numerical_flux_uo[c]), Ui, ni);

  // flux applied to the connected channel; the momentum flux carries the HLLC contact pressure,
  // so the channel sees the junction pressure at its boundary
  _flux[c].resize(THMVACE1D::N_FLUX_OUTPUTS);
  _flux[c][THMVACE1D::MASS] = flux_3d[THMVACE3D::MASS] * nJi_dot_di;
  _flux[c][THMVACE1D::MOMENTUM] = flux_3d[THMVACE3D::MOM_NORM];
  _flux[c][THMVACE1D::ENERGY] = flux_3d[THMVACE3D::ENERGY] * nJi_dot_di;

  // junction residuals: mass and energy only
  _residual[StagnationJunction1Phase::RHOV_INDEX] -= -flux_3d[THMVACE3D::RHOA];
  _residual[StagnationJunction1Phase::RHOEV_INDEX] -= -flux_3d[THMVACE3D::RHOEA];
}

std::vector<ADReal>
ADStagnationJunction1PhaseUserObject::compute3DFlux(const ADNumericalFlux3EqnBase & numerical_flux,
                                                    const std::vector<ADReal> & Ui,
                                                    const RealVectorValue & ni) const
{
  const RealVectorValue nJi = -ni;

  RealVectorValue t1, t2;
  THM::computeOrthogonalDirections(nJi, t1, t2);

  const auto & Ai = Ui[THMVACE3D::AREA];

  const auto & rhoV = _cached_junction_var_values[StagnationJunction1Phase::RHOV_INDEX];
  const auto & rhoEV = _cached_junction_var_values[StagnationJunction1Phase::RHOEV_INDEX];

  // the junction state in the numerical flux has the junction density and internal energy, but the
  // velocity of the connected cell:
  // - a junction state at rest acts as an artificial drag on a weak natural circulation flow;
  // - with no velocity jump the Riemann problem has no acoustic pressure jump rho * c * du, so
  // there is also no acoustic brake, and the hydrostatic face correction is needed to keep a fluid
  // at rest
  const ADReal rhoJ = rhoV / _volume;
  const ADReal eJ = rhoEV / rhoV;
  const ADRealVectorValue veci(Ui[THMVACE3D::RHOUA] / Ui[THMVACE3D::RHOA],
                               Ui[THMVACE3D::RHOVA] / Ui[THMVACE3D::RHOA],
                               Ui[THMVACE3D::RHOWA] / Ui[THMVACE3D::RHOA]);

  std::vector<ADReal> UJi(THMVACE3D::N_FLUX_INPUTS, 0.);
  UJi[THMVACE3D::RHOA] = rhoJ * Ai;
  UJi[THMVACE3D::RHOUA] = rhoJ * veci(0) * Ai;
  UJi[THMVACE3D::RHOVA] = rhoJ * veci(1) * Ai;
  UJi[THMVACE3D::RHOWA] = rhoJ * veci(2) * Ai;
  UJi[THMVACE3D::RHOEA] = rhoJ * (eJ + 0.5 * veci * veci) * Ai;
  UJi[THMVACE3D::AREA] = Ai;

  std::vector<ADReal> FL, FR;
  numerical_flux.calcFlux(UJi, Ui, nJi, t1, t2, FL, FR);

  return FL;
}

void
ADStagnationJunction1PhaseUserObject::finalize()
{
  ADVolumeJunctionBaseUserObject::finalize();
  for (unsigned int i = 0; i < _n_scalar_eq; i++)
    comm().sum(_residual[i]);
}

std::vector<const MooseVariableBase *>
ADStagnationJunction1PhaseUserObject::getFlowChannelVariables() const
{
  std::vector<const MooseVariableBase *> vars(THMVACE1D::N_FLUX_OUTPUTS);
  vars[THMVACE1D::RHOA] = getVar("rhoA", 0);
  vars[THMVACE1D::RHOUA] = getVar("rhouA", 0);
  vars[THMVACE1D::RHOEA] = getVar("rhoEA", 0);
  return vars;
}

std::vector<const MooseVariableBase *>
ADStagnationJunction1PhaseUserObject::getJunctionVariables() const
{
  std::vector<const MooseVariableBase *> vars(StagnationJunction1Phase::N_EQ);
  vars[StagnationJunction1Phase::RHOV_INDEX] = getJunctionVar("rhoV", 0);
  vars[StagnationJunction1Phase::RHOEV_INDEX] = getJunctionVar("rhoEV", 0);
  return vars;
}
