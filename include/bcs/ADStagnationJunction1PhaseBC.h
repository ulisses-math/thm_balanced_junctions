//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADOneDIntegratedBC.h"

class ADVolumeJunctionBaseUserObject;

/**
 * Adds the boundary fluxes of a stagnation junction to a connected flow channel.
 *
 * Notes:
 * - the form is the same as ADVolumeJunction1PhaseBC;
 * - the junction user object is taken with its base type, ADVolumeJunctionBaseUserObject, so that
 * the BC works with ADStagnationJunction1PhaseUserObject.
 */
class ADStagnationJunction1PhaseBC : public ADOneDIntegratedBC
{
public:
  ADStagnationJunction1PhaseBC(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;

  /// Maps the variable numbers of the connected flow channel to flux equation indices
  std::map<unsigned int, unsigned int> getFlowChannelIndexMapping() const;

  /// Connection index of the connected flow channel
  const unsigned int _connection_index;
  /// Junction user object (base type, which provides getFlux)
  const ADVolumeJunctionBaseUserObject & _junction_uo;

  /// Cross-sectional area, elemental
  const ADVariableValue & _A_elem;
  /// Cross-sectional area, linear
  const ADVariableValue & _A_linear;

  /// Variable number of rho*A
  const unsigned int _rhoA_jvar;
  /// Variable number of rho*u*A
  const unsigned int _rhouA_jvar;
  /// Variable number of rho*E*A
  const unsigned int _rhoEA_jvar;

  /// Map of coupled flow variable index to flux equation index
  const std::map<unsigned int, unsigned int> _flow_channel_jvar_map;
  /// Flux equation index for the variable of this BC
  const unsigned int _equation_index;

public:
  static InputParameters validParams();
};
