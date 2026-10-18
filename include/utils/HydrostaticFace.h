//* This file is part of thm_balanced_junctions, an application based on the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADReal.h"
#include "SinglePhaseFluidProperties.h"

/**
 * Hydrostatic extrapolation of a flow channel cell state to a junction face.
 *
 * The 0-D junctions of the thermal_hydraulics module build the channel state at the junction face
 * from the cell average of the adjacent cell, that is, from the pressure at the cell center, half a
 * cell above or below the junction. The junction itself has no elevation and no gravity term. Under
 * gravity, half a cell of hydrostatic head is left without a counterpart, and a fluid at rest gets
 * a pressure error of rho * g * dz / 4 in each cell next to a junction.
 *
 * The function in this namespace moves the cell state to the face:
 * - p_face = p_c + rho_c * g . (x_face - x_c), with the density of the cell;
 * - rho and u are kept, and the specific internal energy is changed at constant specific volume,
 * e_face = e_c + (p_face - p_c) / (dp/de)_v.
 *
 * Notes:
 * - the result is exact for the stiffened gas equation of state (p is linear in e) and first order
 * in (p_face - p_c) for a general equation of state;
 * - no inverse equation of state call is needed, so no round-off cancellation is added;
 * - for a linear hydrostatic profile the face pressure is exact, so a fluid at rest stays at rest.
 */
namespace HydrostaticFace
{

/**
 * Replaces rho*E*A of a cell by the value of the hydrostatically extrapolated face state
 *
 * @param fp Single-phase fluid properties
 * @param A Cell-average area of the junction-adjacent cell (the area used in the flux)
 * @param g_dot_dx Gravity vector times (x_face - x_c) [m^2/s^2], positive when the face is below
 * the cell center
 * @param rhoA rho*A of the cell (the face has the same density)
 * @param rhouA rho*u*A of the cell (the face has the same velocity)
 * @param rhoEA rho*E*A of the cell on input, rho*E*A of the face on output
 */
inline void
extrapolate(const SinglePhaseFluidProperties & fp,
            const ADReal & A,
            const Real g_dot_dx,
            const ADReal & rhoA,
            const ADReal & rhouA,
            ADReal & rhoEA)
{
  // horizontal connection or no gravity: the face state is the cell state
  if (g_dot_dx == 0.0)
    return;

  const ADReal rho = rhoA / A;
  const ADReal u = rhouA / rhoA;
  const ADReal v = 1.0 / rho;
  const ADReal e = rhoEA / rhoA - 0.5 * u * u;

  ADReal p, dp_dv, dp_de;
  fp.p_from_v_e(v, e, p, dp_dv, dp_de);

  const ADReal dp = rho * g_dot_dx;
  const ADReal e_face = e + dp / dp_de;

  rhoEA = rho * (e_face + 0.5 * u * u) * A;
}

}
