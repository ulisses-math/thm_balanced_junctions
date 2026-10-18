//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* computeFluxesAndResiduals() and finalize() are adapted from ADVolumeJunction1PhaseUserObject and
//* ADJunctionParallelChannels1PhaseUserObject of the MOOSE thermal_hydraulics module (Battelle
//* Energy Alliance, LLC, LGPL 2.1).

#include "ADHydrostaticJunctionParallelChannels1PhaseUserObject.h"
#include "SinglePhaseFluidProperties.h"
#include "THMIndicesVACE.h"
#include "HydrostaticFace.h"
#include "VolumeJunction1Phase.h"
#include "ADNumericalFlux3EqnBase.h"
#include "Numerics.h"
#include "THMUtils.h"
#include "metaphysicl/parallel_numberarray.h"
#include "metaphysicl/parallel_dualnumber.h"
#include "metaphysicl/parallel_semidynamicsparsenumberarray.h"
#include "libmesh/parallel_algebra.h"

registerMooseObject("ThmBalancedJunctionsApp",
                    ADHydrostaticJunctionParallelChannels1PhaseUserObject);

InputParameters
ADHydrostaticJunctionParallelChannels1PhaseUserObject::validParams()
{
  InputParameters params = ADJunctionParallelChannels1PhaseUserObject::validParams();

  params.addParam<RealVectorValue>(
      "gravity_vector", RealVectorValue(0, 0, -9.81), "Gravitational acceleration vector [m/s^2]");
  params.addParam<bool>("hydrostatic_correction",
                        true,
                        "Extrapolate the channel cell state hydrostatically to the junction face");
  params.addParam<Real>("mdot_tolerance",
                        1e-3,
                        "Band on rho*u*A of connection 0 [kg/s]: below it, the flow direction used "
                        "by the inlet/outlet classification and by the wall term is held at its "
                        "previous value");

  params.addClassDescription("Computes and caches flux and residual vectors for a 1-phase junction "
                             "that connects parallel flow channels, with the channel state "
                             "extrapolated hydrostatically to the junction face and one held flow "
                             "direction for the wall term");
  return params;
}

ADHydrostaticJunctionParallelChannels1PhaseUserObject::
    ADHydrostaticJunctionParallelChannels1PhaseUserObject(const InputParameters & params)
  : ADJunctionParallelChannels1PhaseUserObject(params),
    _gravity_vector(getParam<RealVectorValue>("gravity_vector")),
    _hydrostatic_correction(getParam<bool>("hydrostatic_correction")),
    _mdot_tol(getParam<Real>("mdot_tolerance")),
    _dirs(_n_connections),
    _stored_rhouA(_n_connections, 0.0),
    _d_flow_held(0, 0, 0),
    _d_flow_initialized(false)
{
}

void
ADHydrostaticJunctionParallelChannels1PhaseUserObject::computeFluxesAndResiduals(
    const unsigned int & c)
{
  // from ADVolumeJunction1PhaseUserObject::computeFluxesAndResiduals()
  const Real din = _normal[c];
  const RealVectorValue di = _dir[0];
  const RealVectorValue ni = di * din;
  const Real nJi_dot_di = -din;

  // change 1: face state instead of the cell state
  ADReal rhoEA_f = _rhoEA[0];
  if (_hydrostatic_correction)
  {
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

  _flux[c].resize(THMVACE1D::N_FLUX_OUTPUTS);
  _flux[c][THMVACE1D::MASS] = flux_3d[THMVACE3D::MASS] * nJi_dot_di;
  _flux[c][THMVACE1D::MOMENTUM] = flux_3d[THMVACE3D::MOM_NORM];
  _flux[c][THMVACE1D::ENERGY] = flux_3d[THMVACE3D::ENERGY] * nJi_dot_di;

  const bool is_primary_connection = (c == 0);
  const auto residual = computeResidual(flux_3d, Ui, ni, is_primary_connection);

  for (const auto i : index_range(_residual))
    _residual[i] += residual[i];

  // from ADJunctionParallelChannels1PhaseUserObject::computeFluxesAndResiduals()
  const auto & rhoV = _cached_junction_var_values[VolumeJunction1Phase::RHOV_INDEX];
  const auto & rhouV = _cached_junction_var_values[VolumeJunction1Phase::RHOUV_INDEX];
  const auto & rhovV = _cached_junction_var_values[VolumeJunction1Phase::RHOVV_INDEX];
  const auto & rhowV = _cached_junction_var_values[VolumeJunction1Phase::RHOWV_INDEX];
  const auto & rhoEV = _cached_junction_var_values[VolumeJunction1Phase::RHOEV_INDEX];

  const ADReal vJ = THM::v_from_rhoA_A(rhoV, _volume);
  const ADRealVectorValue rhouV_vec(rhouV, rhovV, rhowV);
  const ADReal rhouV2 = rhouV_vec * rhouV_vec;
  const ADReal eJ = rhoEV / rhoV - 0.5 * rhouV2 / (rhoV * rhoV);
  const ADReal pJ = _fp.p_from_v_e(vJ, eJ);

  _residual[VolumeJunction1Phase::RHOUV_INDEX] -= pJ * ni(0) * _A[0];
  _residual[VolumeJunction1Phase::RHOVV_INDEX] -= pJ * ni(1) * _A[0];
  _residual[VolumeJunction1Phase::RHOWV_INDEX] -= pJ * ni(2) * _A[0];

  if (!THM::areParallelVectors(_dir_c0, _dir[0]))
    mooseError(_component_name,
               ": Connected flow channels are not parallel, use VolumeJunction1Phase "
               "component instead.");

  _areas[c] = _A[0];

  // change 3: only the data are stored here; the flow direction and the inlet/outlet flags are
  // set in finalize()
  _dirs[c] = _dir[0];
  _stored_rhouA[c] = MetaPhysicL::raw_value(_rhouA[0]);

  // change 2: wall pressure at the junction elevation (face pressure), not the cell-center
  // material property p
  {
    const ADReal rho_f = _rhoA[0] / _A[0];
    const ADReal u_f = _rhouA[0] / _rhoA[0];
    const ADReal v_f = 1.0 / rho_f;
    const ADReal e_f = rhoEA_f / _rhoA[0] - 0.5 * u_f * u_f;
    const ADReal p_f = _fp.p_from_v_e(v_f, e_f);
    _stored_pA[c] = p_f * _A[0];
  }
}

void
ADHydrostaticJunctionParallelChannels1PhaseUserObject::threadJoin(const UserObject & uo)
{
  ADJunctionParallelChannels1PhaseUserObject::threadJoin(uo);

  const auto & other =
      static_cast<const ADHydrostaticJunctionParallelChannels1PhaseUserObject &>(uo);
  for (unsigned int i = 0; i < other._connection_indices.size(); i++)
  {
    const unsigned int c = other._connection_indices[i];
    _dirs[c] = other._dirs[c];
    _stored_rhouA[c] = other._stored_rhouA[c];
  }
}

void
ADHydrostaticJunctionParallelChannels1PhaseUserObject::finalize()
{
  // from ADJunctionParallelChannels1PhaseUserObject::finalize(): broadcasts
  for (unsigned int i = 0; i < _n_connections; i++)
  {
    processor_id_type owner_proc = _processor_ids[i];
    comm().broadcast(_areas[i], owner_proc, true);
    comm().broadcast(_stored_pA[i], owner_proc, true);
    comm().broadcast(_dirs[i], owner_proc, true);
    comm().broadcast(_stored_rhouA[i], owner_proc, true);
  }

  // change 3: one flow direction for the whole junction, from connection 0, held while
  // |rho*u*A| of connection 0 is below the tolerance. The official class uses
  // _d_flow = dir_0 * sign(rhouA_0) (or dir_0 when rhouA_0 is zero) and sets the inlet flags from
  // the sign of rhouA of each connection.
  if (!_d_flow_initialized)
  {
    _d_flow_held = _dirs[0];
    _d_flow_initialized = true;
  }
  if (std::abs(_stored_rhouA[0]) > _mdot_tol)
    _d_flow_held = _dirs[0] * (_stored_rhouA[0] > 0 ? 1.0 : -1.0);
  _d_flow = _d_flow_held;

  for (unsigned int c = 0; c < _n_connections; c++)
  {
    // fluid moving along _d_flow enters the junction through connection c when the outward
    // normal of the channel at the junction (dir_c * normal_c) points along _d_flow
    const Real along = (_dirs[c] * _d_flow) * _normal[c];
    _is_inlet[c] = (along >= 0);
  }

  // from ADJunctionParallelChannels1PhaseUserObject::finalize(): wall term
  ADReal Ain_total = 0;
  ADReal pAin_total = 0;
  ADReal pAout_total = 0;
  ADReal Aout_total = 0;
  _c_in.clear();
  _c_out.clear();
  for (unsigned int c = 0; c < _n_connections; c++)
  {
    if (_is_inlet[c] == true)
    {
      Ain_total += _areas[c];
      pAin_total += _stored_pA[c];
      _c_in.push_back(c);
    }
    else
    {
      Aout_total += _areas[c];
      pAout_total += _stored_pA[c];
      _c_out.push_back(c);
    }
  }

  ADReal p_wall = 0;
  ADReal A_wall = 0;
  RealVectorValue d_wall;

  if (Aout_total > Ain_total)
  {
    if (Aout_total > 1e-15)
    {
      p_wall = pAout_total / Aout_total;
      _c_wall = _c_out;
    }
    else
    {
      p_wall = 0;
    }
    A_wall = Aout_total - Ain_total;
    d_wall = _d_flow;
  }
  else
  {
    if (Ain_total > 1e-15)
    {
      p_wall = pAin_total / Ain_total;
      _c_wall = _c_in;
    }
    else
    {
      p_wall = 0;
    }
    A_wall = Ain_total - Aout_total;
    d_wall = -_d_flow;
  }

  for (unsigned int i = 0; i < _n_scalar_eq; i++)
    comm().sum(_residual[i]);

  _residual[VolumeJunction1Phase::RHOUV_INDEX] -= d_wall(0) * p_wall * A_wall;
  _residual[VolumeJunction1Phase::RHOVV_INDEX] -= d_wall(1) * p_wall * A_wall;
  _residual[VolumeJunction1Phase::RHOWV_INDEX] -= d_wall(2) * p_wall * A_wall;
}
