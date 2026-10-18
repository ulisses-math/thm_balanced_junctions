//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADVolumeJunction1PhaseUserObject.h"

/**
 * Computes and caches flux and residual vectors for a 1-phase volume junction, with the channel
 * state extrapolated hydrostatically to the junction face.
 *
 * This is ADVolumeJunction1PhaseUserObject with one change:
 * - change 1: the channel state at the junction face is extrapolated hydrostatically from the cell
 * center to the face (HydrostaticFace.h).
 *
 * Notes:
 * - the official class has no wall term and no inlet/outlet classification, so changes 2 and 3 of
 * ADHydrostaticJunctionParallelChannels1PhaseUserObject do not apply here;
 * - hydrostatic_correction = false gives the official face state;
 * - everything else (0-D vector momentum, velocity scaling, form loss) is the official formulation,
 * evaluated with the face state;
 * - the channel ends connected to the junction are assumed to be at the same elevation (one
 * junction pressure serves all faces), and HydrostaticVolumeJunction1Phase checks it.
 */
class ADHydrostaticVolumeJunction1PhaseUserObject : public ADVolumeJunction1PhaseUserObject
{
public:
  ADHydrostaticVolumeJunction1PhaseUserObject(const InputParameters & params);

protected:
  virtual void computeFluxesAndResiduals(const unsigned int & c) override;

  /// Gravitational acceleration vector
  const RealVectorValue _gravity_vector;
  /// Extrapolate the cell state hydrostatically to the junction face?
  const bool _hydrostatic_correction;

public:
  static InputParameters validParams();
};
