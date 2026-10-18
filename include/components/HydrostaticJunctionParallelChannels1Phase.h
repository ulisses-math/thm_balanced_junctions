//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "JunctionParallelChannels1Phase.h"

/**
 * Junction between parallel 1-phase flow channels, hydrostatically balanced.
 *
 * Notes:
 * - this is JunctionParallelChannels1Phase with the user object
 * ADHydrostaticJunctionParallelChannels1PhaseUserObject, which extrapolates the channel state
 * hydrostatically from the cell center to the junction face, uses the face pressure in the wall
 * term and one held flow direction for the inlet/outlet classification;
 * - the input syntax is the same as for JunctionParallelChannels1Phase, plus gravity_vector
 * (usually taken from GlobalParams), hydrostatic_correction and mdot_tolerance;
 * - hydrostatic_correction = false gives the cell-center face state and wall pressure of the
 * official component, but the held flow direction stays active.
 */
class HydrostaticJunctionParallelChannels1Phase : public JunctionParallelChannels1Phase
{
public:
  HydrostaticJunctionParallelChannels1Phase(const InputParameters & params);

protected:
  virtual void check() const override;
  virtual void buildVolumeJunctionUserObject() override;

public:
  static InputParameters validParams();
};
