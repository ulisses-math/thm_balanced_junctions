//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VolumeJunction1Phase.h"

/**
 * Junction connecting one or more 1-phase flow channels, hydrostatically balanced.
 *
 * Notes:
 * - this is VolumeJunction1Phase with the user object ADHydrostaticVolumeJunction1PhaseUserObject,
 * which extrapolates the channel state hydrostatically from the cell center to the junction face;
 * - the input syntax is the same as for VolumeJunction1Phase, plus gravity_vector (usually taken
 * from GlobalParams) and hydrostatic_correction (default true; false gives the official
 * component).
 */
class HydrostaticVolumeJunction1Phase : public VolumeJunction1Phase
{
public:
  HydrostaticVolumeJunction1Phase(const InputParameters & params);

protected:
  virtual void check() const override;
  virtual void buildVolumeJunctionUserObject() override;

public:
  static InputParameters validParams();
};
