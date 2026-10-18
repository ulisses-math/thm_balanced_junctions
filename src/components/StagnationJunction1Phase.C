//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Parts of this file are adapted from VolumeJunction1Phase of the MOOSE thermal_hydraulics module
//* (Battelle Energy Alliance, LLC, LGPL 2.1).

#include "StagnationJunction1Phase.h"
#include "FlowModelSinglePhase.h"
#include "FlowChannel1PhaseBase.h"
#include "THMMesh.h"

registerMooseObject("ThmBalancedJunctionsApp", StagnationJunction1Phase);

const unsigned int StagnationJunction1Phase::N_EQ = 2;

InputParameters
StagnationJunction1Phase::validParams()
{
  InputParameters params = FlowJunction1Phase::validParams();

  params.addRequiredParam<Real>("volume", "Volume of the junction [m^3]");
  params.addRequiredParam<Point>("position", "Spatial position of the center of the junction [m]");

  params.addParam<FunctionName>("initial_p", "Initial pressure [Pa]");
  params.addParam<FunctionName>("initial_T", "Initial temperature [K]");
  params.addParam<FunctionName>("initial_vel_x", "Initial velocity in x-direction [m/s]");
  params.addParam<FunctionName>("initial_vel_y", "Initial velocity in y-direction [m/s]");
  params.addParam<FunctionName>("initial_vel_z", "Initial velocity in z-direction [m/s]");

  params.addParam<Real>("scaling_factor_rhoV", 1.0, "Scaling factor for rho*V [-]");
  params.addParam<Real>("scaling_factor_rhoEV", 1.0, "Scaling factor for rho*E*V [-]");

  params.addParam<Real>(
      "K", 0., "Form loss factor [-], accepted but not applied (no momentum equation)");
  params.addParam<Real>("A_ref", "Reference area [m^2], not used");

  params.addParam<RealVectorValue>("gravity_vector",
                                   RealVectorValue(0, 0, -9.81),
                                   "Gravitational acceleration vector [m/s^2], used to extrapolate "
                                   "the channel cell state hydrostatically to the junction face");
  params.addParam<bool>("hydrostatic_correction",
                        true,
                        "Extrapolate the channel cell state hydrostatically from the cell center "
                        "to the junction face; false gives the cell-center coupling of the "
                        "official junctions");

  params.declareControllable("K");
  params.addClassDescription(
      "Junction between 1-phase flow channels that conserves mass and energy "
      "only (no 0-D momentum equation), hydrostatically balanced");

  return params;
}

StagnationJunction1Phase::StagnationJunction1Phase(const InputParameters & params)
  : FlowJunction1Phase(params),

    _volume(getParam<Real>("volume")),
    _position(getParam<Point>("position")),

    _scaling_factor_rhoV(getParam<Real>("scaling_factor_rhoV")),
    _scaling_factor_rhoEV(getParam<Real>("scaling_factor_rhoEV")),

    _rhoV_var_name("rhoV"),
    _rhoEV_var_name("rhoEV"),

    _K(getParam<Real>("K")),
    _A_ref(isParamValid("A_ref") ? getParam<Real>("A_ref") : _zero)
{
}

void
StagnationJunction1Phase::setupMesh()
{
  FlowJunction1Phase::setupMesh();

  auto * node = addNode(_position);
  auto * elem = addNodeElement(node->id());
  _junction_subdomain_id = mesh().getNextSubdomainId();
  elem->subdomain_id() = _junction_subdomain_id;
  setSubdomainInfo(_junction_subdomain_id, name());

  const auto & elem_ids = getConnectedElementIDs();
  for (unsigned int i = 0; i < elem_ids.size(); i++)
    getTHMProblem().augmentSparsity(elem_ids[i], elem->id());
}

void
StagnationJunction1Phase::check() const
{
  FlowJunction1Phase::check();

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

  bool ics_set =
      getTHMProblem().hasInitialConditionsFromFile() ||
      (isParamValid("initial_p") && isParamValid("initial_T") && isParamValid("initial_vel_x") &&
       isParamValid("initial_vel_y") && isParamValid("initial_vel_z"));

  if (!ics_set && !_app.isRestarting())
  {
    const std::vector<std::string> ic_params{
        "initial_p", "initial_T", "initial_vel_x", "initial_vel_y", "initial_vel_z"};
    std::ostringstream oss;
    for (const auto & ic_param : ic_params)
      if (!isParamValid(ic_param))
        oss << " " << ic_param;

    logError("The following initial condition parameters are missing:", oss.str());
  }
}

void
StagnationJunction1Phase::addVariables()
{
  addJunctionVariable(true, _rhoV_var_name, _scaling_factor_rhoV);
  addJunctionVariable(true, _rhoEV_var_name, _scaling_factor_rhoEV);

  if (isParamValid("initial_p") && isParamValid("initial_T") && isParamValid("initial_vel_x") &&
      isParamValid("initial_vel_y") && isParamValid("initial_vel_z"))
  {
    addStagnationJunctionIC(_rhoV_var_name, "rhoV");
    addStagnationJunctionIC(_rhoEV_var_name, "rhoEV");
  }
}

void
StagnationJunction1Phase::buildStagnationJunctionUserObject()
{
  ExecFlagEnum execute_on(MooseUtils::getDefaultExecFlagEnum());
  execute_on = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};

  const std::string class_name = "ADStagnationJunction1PhaseUserObject";
  InputParameters params = _factory.getValidParams(class_name);
  params.set<bool>("use_scalar_variables") = false;
  params.set<subdomain_id_type>("junction_subdomain_id") = _junction_subdomain_id;
  params.set<std::vector<BoundaryName>>("boundary") = _boundary_names;
  params.set<std::vector<Real>>("normals") = _normals;
  params.set<std::vector<processor_id_type>>("processor_ids") = getConnectedProcessorIDs();
  params.set<std::vector<UserObjectName>>("numerical_flux_names") = _numerical_flux_names;
  params.set<Real>("volume") = _volume;
  params.set<std::vector<VariableName>>("A") = {FlowModel::AREA};
  params.set<std::vector<VariableName>>("rhoA") = {FlowModelSinglePhase::RHOA};
  params.set<std::vector<VariableName>>("rhouA") = {FlowModelSinglePhase::RHOUA};
  params.set<std::vector<VariableName>>("rhoEA") = {FlowModelSinglePhase::RHOEA};
  params.set<std::vector<VariableName>>("rhoV") = {_rhoV_var_name};
  params.set<std::vector<VariableName>>("rhoEV") = {_rhoEV_var_name};
  params.set<Real>("K") = _K;
  params.set<Real>("A_ref") = _A_ref;
  params.set<UserObjectName>("fp") = _fp_name;
  params.set<RealVectorValue>("gravity_vector") = getParam<RealVectorValue>("gravity_vector");
  params.set<bool>("hydrostatic_correction") = getParam<bool>("hydrostatic_correction");
  params.set<ExecFlagEnum>("execute_on") = execute_on;
  getTHMProblem().addUserObject(class_name, _junction_uo_name, params);
  connectObject(params, _junction_uo_name, "K");
}

void
StagnationJunction1Phase::addMooseObjects()
{
  buildStagnationJunctionUserObject();

  // boundary condition on each connected flow channel
  for (std::size_t i = 0; i < _boundary_names.size(); i++)
  {
    const std::vector<NonlinearVariableName> var_names = {
        FlowModelSinglePhase::RHOA, FlowModelSinglePhase::RHOUA, FlowModelSinglePhase::RHOEA};
    for (std::size_t j = 0; j < var_names.size(); j++)
    {
      const std::string class_name = "ADStagnationJunction1PhaseBC";
      InputParameters params = _factory.getValidParams(class_name);
      params.set<std::vector<BoundaryName>>("boundary") = {_boundary_names[i]};
      params.set<Real>("normal") = _normals[i];
      params.set<NonlinearVariableName>("variable") = var_names[j];
      params.set<UserObjectName>("stagnation_junction_uo") = _junction_uo_name;
      params.set<unsigned int>("connection_index") = i;
      params.set<std::vector<VariableName>>("A_elem") = {FlowModel::AREA};
      params.set<std::vector<VariableName>>("A_linear") = {FlowModel::AREA_LINEAR};
      params.set<std::vector<VariableName>>("rhoA") = {FlowModelSinglePhase::RHOA};
      params.set<std::vector<VariableName>>("rhouA") = {FlowModelSinglePhase::RHOUA};
      params.set<std::vector<VariableName>>("rhoEA") = {FlowModelSinglePhase::RHOEA};
      params.set<bool>("implicit") = getTHMProblem().getImplicitTimeIntegrationFlag();
      getTHMProblem().addBoundaryCondition(
          class_name, genName(name(), i, var_names[j] + ":" + class_name), params);
    }
  }

  // kernels of the junction equations (mass and energy)
  std::vector<NonlinearVariableName> var_names(N_EQ);
  var_names[RHOV_INDEX] = _rhoV_var_name;
  var_names[RHOEV_INDEX] = _rhoEV_var_name;
  for (std::size_t i = 0; i < N_EQ; i++)
  {
    {
      const std::string class_name = "ADTimeDerivative";
      InputParameters params = _factory.getValidParams(class_name);
      params.set<NonlinearVariableName>("variable") = var_names[i];
      params.set<std::vector<SubdomainName>>("block") = getSubdomainNames();
      getTHMProblem().addKernel(class_name, genName(name(), var_names[i], "td"), params);
    }
    {
      const std::string class_name = "ADVolumeJunctionAdvectionKernel";
      InputParameters params = _factory.getValidParams(class_name);
      params.set<NonlinearVariableName>("variable") = var_names[i];
      params.set<UserObjectName>("volume_junction_uo") = _junction_uo_name;
      params.set<unsigned int>("equation_index") = i;
      params.set<std::vector<SubdomainName>>("block") = getSubdomainNames();
      getTHMProblem().addKernel(class_name, genName(name(), var_names[i], "vja_sk"), params);
    }
  }

  // empty material, because THM needs a material on every block
  {
    const std::string class_name = "GenericConstantMaterial";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<std::vector<SubdomainName>>("block") = getSubdomainNames();
    params.set<std::vector<std::string>>("prop_names") = {};
    params.set<std::vector<Real>>("prop_values") = {};
    getTHMProblem().addMaterial(class_name, genName(name(), "dummy_mat"), params);
  }
}

void
StagnationJunction1Phase::addJunctionVariable(bool is_nonlinear,
                                              const VariableName & var,
                                              Real scaling_factor)
{
  auto & problem = getTHMProblem();

  const libMesh::FEType fe_type(CONSTANT, MONOMIAL);
  const auto & subdomains = getSubdomainNames();

  if (is_nonlinear)
    problem.addSimVariable(is_nonlinear, var, fe_type, subdomains, scaling_factor);
  else
    problem.addSimVariable(is_nonlinear, var, fe_type, subdomains);
}

void
StagnationJunction1Phase::addStagnationJunctionIC(const VariableName & var,
                                                  const std::string & quantity)
{
  const std::string class_name = "VolumeJunction1PhaseIC";
  InputParameters params = _factory.getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = getSubdomainNames();
  params.set<VariableName>("variable") = var;
  params.set<MooseEnum>("quantity") = quantity;
  params.applySpecificParameters(parameters(),
                                 {"initial_p",
                                  "initial_T",
                                  "initial_vel_x",
                                  "initial_vel_y",
                                  "initial_vel_z",
                                  "volume",
                                  "position"});
  params.set<UserObjectName>("fluid_properties") = _fp_name;
  getTHMProblem().addSimInitialCondition(class_name, genName(name(), var, "ic"), params);
}
