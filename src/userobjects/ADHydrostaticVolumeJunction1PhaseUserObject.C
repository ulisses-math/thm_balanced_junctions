//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* computeFluxesAndResiduals() is adapted from ADVolumeJunction1PhaseUserObject of the MOOSE
//* thermal_hydraulics module (Battelle Energy Alliance, LLC, LGPL 2.1).

#include "ADHydrostaticVolumeJunction1PhaseUserObject.h"
#include "SinglePhaseFluidProperties.h"
#include "THMIndicesVACE.h"
#include "HydrostaticFace.h"
#include "ADNumericalFlux3EqnBase.h"

registerMooseObject("ThmBalancedJunctionsApp", ADHydrostaticVolumeJunction1PhaseUserObject);

InputParameters
ADHydrostaticVolumeJunction1PhaseUserObject::validParams()
{
  InputParameters params = ADVolumeJunction1PhaseUserObject::validParams();

  params.addParam<RealVectorValue>(
      "gravity_vector", RealVectorValue(0, 0, -9.81), "Gravitational acceleration vector [m/s^2]");
  params.addParam<bool>("hydrostatic_correction",
                        true,
                        "Extrapolate the channel cell state hydrostatically to the junction face");

  params.addClassDescription("Computes and caches flux and residual vectors for a 1-phase volume "
                             "junction, with the channel state extrapolated hydrostatically to "
                             "the junction face");
  return params;
}

ADHydrostaticVolumeJunction1PhaseUserObject::ADHydrostaticVolumeJunction1PhaseUserObject(
    const InputParameters & params)
  : ADVolumeJunction1PhaseUserObject(params),
    _gravity_vector(getParam<RealVectorValue>("gravity_vector")),
    _hydrostatic_correction(getParam<bool>("hydrostatic_correction"))
{
}

void
ADHydrostaticVolumeJunction1PhaseUserObject::computeFluxesAndResiduals(const unsigned int & c)
{
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
}
