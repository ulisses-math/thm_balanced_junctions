//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* buildVolumeJunctionUserObject() is adapted from JunctionParallelChannels1Phase of the MOOSE
//* thermal_hydraulics module (Battelle Energy Alliance, LLC, LGPL 2.1).

#include "HydrostaticJunctionParallelChannels1Phase.h"
#include "FlowModelSinglePhase.h"
#include "FlowChannel1PhaseBase.h"

registerMooseObject("ThmBalancedJunctionsApp", HydrostaticJunctionParallelChannels1Phase);

InputParameters
HydrostaticJunctionParallelChannels1Phase::validParams()
{
  InputParameters params = JunctionParallelChannels1Phase::validParams();

  params.addParam<RealVectorValue>("gravity_vector",
                                   RealVectorValue(0, 0, -9.81),
                                   "Gravitational acceleration vector [m/s^2], used to extrapolate "
                                   "the channel cell state hydrostatically to the junction face");
  params.addParam<bool>("hydrostatic_correction",
                        true,
                        "Extrapolate the channel cell state hydrostatically from the cell center "
                        "to the junction face, also for the wall pressure of the area change "
                        "term; false gives the cell-center values of the official "
                        "JunctionParallelChannels1Phase (the held flow direction of the wall term "
                        "stays active)");
  params.addParam<Real>("mdot_tolerance",
                        1e-3,
                        "Band on rho*u*A of the first connection [kg/s]: below it, the flow "
                        "direction used by the wall term is held at its previous value");

  params.addClassDescription(
      "Junction between 1-phase flow channels that are parallel, hydrostatically balanced");
  return params;
}

HydrostaticJunctionParallelChannels1Phase::HydrostaticJunctionParallelChannels1Phase(
    const InputParameters & params)
  : JunctionParallelChannels1Phase(params)
{
}

void
HydrostaticJunctionParallelChannels1Phase::check() const
{
  JunctionParallelChannels1Phase::check();

  // preconditions of the hydrostatic face correction:
  // - the connected channels use slope reconstruction, otherwise the interior faces are not
  // balanced and the face correction cannot balance them;
  // - the junction and the connected channels use the same gravity vector;
  // - the channel ends are at the same elevation (the junction has one pressure for all faces)
  if (getParam<bool>("hydrostatic_correction"))
  {
    const RealVectorValue g = getParam<RealVectorValue>("gravity_vector");
    for (const auto & conn : _connections)
    {
      if (hasComponentByName<FlowChannel1PhaseBase>(conn._component_name))
      {
        const auto & ch = getComponentByName<FlowChannel1PhaseBase>(conn._component_name);
        if (ch.getParam<MooseEnum>("rdg_slope_reconstruction") == "none")
          logWarning("hydrostatic_correction = true but the connected channel '",
                     conn._component_name,
                     "' uses rdg_slope_reconstruction = none: the interior faces are not "
                     "reconstructed and the junction face correction cannot balance them.");
        if (!ch.getParam<RealVectorValue>("gravity_vector").absolute_fuzzy_equals(g, 1e-12))
          logError("gravity_vector of the junction differs from that of the connected channel '",
                   conn._component_name,
                   "'; set gravity_vector in [GlobalParams] or on both.");
      }
    }
    if (_positions.size() > 1 && g.norm() > 0)
    {
      const RealVectorValue gu = g / g.norm();
      const Real z0 = _positions[0] * gu;
      for (const auto & x : _positions)
        if (std::abs(x * gu - z0) > 1e-6)
          logError("the channel ends connected to this junction are not at the same elevation "
                   "along the gravity vector; the hydrostatic face correction assumes they are.");
    }
  }
}

void
HydrostaticJunctionParallelChannels1Phase::buildVolumeJunctionUserObject()
{
  ExecFlagEnum execute_on(MooseUtils::getDefaultExecFlagEnum());
  execute_on = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};

  {
    const std::string class_name = "ADHydrostaticJunctionParallelChannels1PhaseUserObject";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<bool>("use_scalar_variables") = false;
    params.set<subdomain_id_type>("junction_subdomain_id") = _junction_subdomain_id;
    params.set<std::vector<BoundaryName>>("boundary") = _boundary_names;
    params.set<std::vector<Real>>("normals") = _normals;
    params.set<std::vector<processor_id_type>>("processor_ids") = getConnectedProcessorIDs();
    params.set<std::vector<UserObjectName>>("numerical_flux_names") = _numerical_flux_names;
    params.set<Real>("volume") = _volume;
    params.set<std::string>("component_name") = name();
    params.set<std::vector<VariableName>>("A") = {FlowModel::AREA};
    params.set<std::vector<VariableName>>("rhoA") = {FlowModelSinglePhase::RHOA};
    params.set<std::vector<VariableName>>("rhouA") = {FlowModelSinglePhase::RHOUA};
    params.set<std::vector<VariableName>>("rhoEA") = {FlowModelSinglePhase::RHOEA};
    params.set<std::vector<VariableName>>("rhoV") = {_rhoV_var_name};
    params.set<std::vector<VariableName>>("rhouV") = {_rhouV_var_name};
    params.set<std::vector<VariableName>>("rhovV") = {_rhovV_var_name};
    params.set<std::vector<VariableName>>("rhowV") = {_rhowV_var_name};
    params.set<std::vector<VariableName>>("rhoEV") = {_rhoEV_var_name};
    params.set<RealVectorValue>("dir_c0") = _directions[0];
    params.set<Real>("K") = _K;
    params.set<Real>("A_ref") = _A_ref;
    params.set<UserObjectName>("fp") = _fp_name;
    params.set<bool>("apply_velocity_scaling") = getParam<bool>("apply_velocity_scaling");
    params.set<RealVectorValue>("gravity_vector") = getParam<RealVectorValue>("gravity_vector");
    params.set<bool>("hydrostatic_correction") = getParam<bool>("hydrostatic_correction");
    params.set<Real>("mdot_tolerance") = getParam<Real>("mdot_tolerance");
    params.set<ExecFlagEnum>("execute_on") = execute_on;
    getTHMProblem().addUserObject(class_name, _junction_uo_name, params);
    connectObject(params, _junction_uo_name, "K");
  }
}
