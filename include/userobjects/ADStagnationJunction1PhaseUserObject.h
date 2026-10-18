//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADVolumeJunctionBaseUserObject.h"

class SinglePhaseFluidProperties;
class ADNumericalFlux3EqnBase;

/**
 * Computes and caches flux and residual vectors for a 1-phase stagnation junction.
 *
 * Notes:
 * - only the mass (rho*V) and energy (rho*E*V) residuals are assembled, there is no 0-D momentum;
 * - the junction state in the numerical flux has the junction pressure and internal energy, but it
 * moves with the velocity of the connected cell, so the junction gives no acoustic pressure jump
 * rho * c * du;
 * - the channel state at the junction face is extrapolated hydrostatically from the cell center to
 * the face (HydrostaticFace.h); with hydrostatic_correction = false the cell-center state is used,
 * as in the official VolumeJunction1Phase.
 */
class ADStagnationJunction1PhaseUserObject : public ADVolumeJunctionBaseUserObject
{
public:
  ADStagnationJunction1PhaseUserObject(const InputParameters & params);

  virtual void finalize() override;

protected:
  virtual void computeFluxesAndResiduals(const unsigned int & c) override;

  virtual std::vector<const MooseVariableBase *> getFlowChannelVariables() const override;
  virtual std::vector<const MooseVariableBase *> getJunctionVariables() const override;

  /**
   * Computes the 3D numerical flux between the junction state and the connected channel state
   */
  std::vector<ADReal> compute3DFlux(const ADNumericalFlux3EqnBase & numerical_flux,
                                    const std::vector<ADReal> & Ui,
                                    const RealVectorValue & ni) const;

  /// Cross-sectional area of connected flow channels
  const ADVariableValue & _A;
  /// rho*A of the connected flow channel
  const ADVariableValue & _rhoA;
  /// rho*u*A of the connected flow channel
  const ADVariableValue & _rhouA;
  /// rho*E*A of the connected flow channel
  const ADVariableValue & _rhoEA;

  /// Form loss coefficient (accepted for input compatibility, not used)
  const Real & _K;
  /// Reference area (accepted for input compatibility, not used)
  const Real & _A_ref;

  /// Gravitational acceleration vector
  const RealVectorValue _gravity_vector;
  /// Extrapolate the cell state hydrostatically to the junction face?
  const bool _hydrostatic_correction;

  /// Single-phase fluid properties
  const SinglePhaseFluidProperties & _fp;

  /// Numerical flux user object for each connection
  std::vector<const ADNumericalFlux3EqnBase *> _numerical_flux_uo;

public:
  static InputParameters validParams();
};
