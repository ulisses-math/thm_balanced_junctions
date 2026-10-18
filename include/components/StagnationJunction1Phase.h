//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FlowJunction1Phase.h"

/**
 * Junction between 1-phase flow channels that conserves mass and energy only, with no 0-D momentum
 * equation, similar to a RELAP5 branch volume.
 *
 * Notes:
 * - the junction variables are rho*V and rho*E*V; the momentum of each connected channel is coupled
 * to the junction through the junction pressure, and gravity acts only in the 1-D channels;
 * - the junction state used in the numerical flux moves with the velocity of the connected cell, so
 * the junction adds no acoustic pressure jump rho * c * du (and also no acoustic brake);
 * - the channel state at the junction face is extrapolated hydrostatically (HydrostaticFace.h);
 * without this correction a fluid at rest gets a spurious loop flow when the cells on the two sides
 * of a junction have different heights;
 * - K is accepted for input compatibility with the official junctions, but it is not applied,
 * because there is no momentum equation; form losses go in the connected channels;
 * - we use it in our reactor model, where the 0-D momentum equation of the official junctions gave
 * convergence problems during the flow reversal.
 */
class StagnationJunction1Phase : public FlowJunction1Phase
{
public:
  StagnationJunction1Phase(const InputParameters & params);

  virtual void addVariables() override;
  virtual void addMooseObjects() override;

  /// Indices of the junction variables and equations (mass and energy only)
  enum StagnationJunction1PhaseIndices
  {
    RHOV_INDEX = 0,
    RHOEV_INDEX = 1
  };
  /// Number of equations of the junction
  static const unsigned int N_EQ;

  UserObjectName getStagnationJunctionUserObjectName() const { return _junction_uo_name; }

protected:
  virtual void setupMesh() override;
  virtual void check() const override;

  /// Builds the user object that computes and stores the fluxes and residuals
  virtual void buildStagnationJunctionUserObject();

  /// Adds a junction variable (constant monomial on the junction subdomain)
  void addJunctionVariable(bool is_nonlinear, const VariableName & var, Real scaling_factor = 1.0);

  /// Adds a VolumeJunction1PhaseIC to initialize a junction quantity from p and T
  void addStagnationJunctionIC(const VariableName & var, const std::string & quantity);

  /// Volume of the junction
  const Real _volume;
  /// Spatial position of the center of the junction
  const Point & _position;

  /// Scaling factor for rho*V
  const Real & _scaling_factor_rhoV;
  /// Scaling factor for rho*E*V
  const Real & _scaling_factor_rhoEV;

  /// rho*V variable name for the junction
  const VariableName _rhoV_var_name;
  /// rho*E*V variable name for the junction
  const VariableName _rhoEV_var_name;

  /// Form loss coefficient (not applied)
  const Real & _K;
  /// Reference area (not used)
  const Real & _A_ref;

public:
  static InputParameters validParams();
};
