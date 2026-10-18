//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADJunctionParallelChannels1PhaseUserObject.h"

/**
 * Computes and caches flux and residual vectors for a 1-phase junction that connects parallel flow
 * channels, with the channel state extrapolated hydrostatically to the junction face.
 *
 * This is ADJunctionParallelChannels1PhaseUserObject with three changes:
 * - change 1: the channel state at the junction face is extrapolated hydrostatically from the cell
 * center to the face (HydrostaticFace.h);
 * - change 2: the wall pressure of the abrupt area change term is the face pressure, not the
 * cell-center pressure;
 * - change 3: the inlet/outlet classification of the connections uses one flow direction for the
 * whole junction (the direction of connection 0, which the official wall term already uses), held
 * at its previous value while |rho*u*A| of connection 0 is below mdot_tolerance.
 *
 * Notes:
 * - change 3 is needed for a fluid at rest: with change 1 the small standing velocities go away,
 * the sign of rho*u*A of each connection becomes round-off noise, and with the official
 * classification all connections can end up on the same side, so the wall term jumps between
 * p * |A_out - A_in| and p * (A_out + A_in) from one Newton iteration to the next;
 * - hydrostatic_correction = false reverts changes 1 and 2 only, change 3 stays active;
 * - everything else (0-D vector momentum, velocity scaling, form loss, the wall term itself) is the
 * official formulation, evaluated with the face state;
 * - as in the official class, the channel ends connected to the junction are assumed to be at the
 * same elevation (one junction pressure serves all faces).
 */
class ADHydrostaticJunctionParallelChannels1PhaseUserObject
  : public ADJunctionParallelChannels1PhaseUserObject
{
public:
  ADHydrostaticJunctionParallelChannels1PhaseUserObject(const InputParameters & params);

  virtual void threadJoin(const UserObject & uo) override;
  virtual void finalize() override;

protected:
  virtual void computeFluxesAndResiduals(const unsigned int & c) override;

  /// Gravitational acceleration vector
  const RealVectorValue _gravity_vector;
  /// Extrapolate the cell state hydrostatically to the junction face?
  const bool _hydrostatic_correction;
  /// Band on rho*u*A of connection 0 below which the flow direction is held [kg/s]
  const Real _mdot_tol;

  /// Channel direction at each connection (stored in execute, broadcast in finalize)
  std::vector<RealVectorValue> _dirs;
  /// rho*u*A at each connection (stored in execute, broadcast in finalize)
  std::vector<Real> _stored_rhouA;
  /// Flow direction held from the previous evaluation
  RealVectorValue _d_flow_held;
  /// Has _d_flow_held been initialized?
  bool _d_flow_initialized;

public:
  static InputParameters validParams();
};
