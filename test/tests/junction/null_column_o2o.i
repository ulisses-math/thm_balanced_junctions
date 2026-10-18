# Column of water at rest with a JunctionOneToOne1Phase at mid-height
#
# The column of null_column.i (7 m, open at the top, closed at the bottom, two vertical channels
# joined at mid-height) with the official JunctionOneToOne1Phase, which has no volume and
# reconstructs the state across the junction.
#
# Notes:
# - the fluid is at rest under gravity, so any flow is a numerical error;
# - the slope reconstruction across the junction (ADJunctionOneToOne1PhaseUserObject::finalize(),
# through SlopeReconstruction1DInterface::getBoundaryElementSlopes()) places the cell on the other
# side at the mirrored position of the interior neighbor, not at its real position;
# - this application has no corrected version of JunctionOneToOne1Phase;
# - JunctionOneToOne1Phase needs the same area in the two channels, so there is no AB parameter;
# - nA and nB are the numbers of cells above and below the junction;
# - scripts/check_hydro_profile.py compares the final profile with the hydrostatic profile, cell by
# cell.

nA = 50
nB = 50
AA = 0.5
g = -9.81
rdg = minmod
H = 7.0
T0 = 313.15
p_top = 101325

[GlobalParams]
  gravity_vector = '0 0 ${g}'
  rdg_slope_reconstruction = ${rdg}
  closures = thm_closure
  fp = water
  initial_T = ${T0}
  initial_p = p_init_fn
  initial_vel = 0
  scaling_factor_1phase = '1 1 1'
[]

[FluidProperties]
  [water]
    type = StiffenedGasFluidProperties
    gamma = 2.35
    cv = 1778.3
    p_inf = 7.456e8
    q = -1.800e5
    q_prime = 0
    mu = 6.53e-4
    k = 0.631
  []
[]

[Closures]
  [thm_closure]
    type = Closures1PhaseTHM
  []
[]

[Functions]
  # hydrostatic initial pressure with rho = 991.3 kg/m^3 (this stiffened gas at 313 K and 1 atm)
  [p_init_fn]
    type = ParsedFunction
    expression = '${p_top} + ${fparse -g} * 991.3 * (${H} - z)'
  []
[]

[Components]
  [top_bc]
    type = Outlet1Phase
    input = 'pipeA:in'
    p = ${p_top}
  []
  [pipeA]
    type = FlowChannel1Phase
    position = '0 0 ${H}'
    orientation = '0 0 -1'
    length = ${fparse H / 2}
    n_elems = ${nA}
    A = ${AA}
    D_h = 0.8
  []
  [jct]
    type = JunctionOneToOne1Phase
    connections = 'pipeA:out pipeB:in'
  []
  [pipeB]
    type = FlowChannel1Phase
    position = '0 0 ${fparse H / 2}'
    orientation = '0 0 -1'
    length = ${fparse H / 2}
    n_elems = ${nB}
    A = ${AA}
    D_h = 0.8
  []
  [bot_bc]
    type = SolidWall1Phase
    input = 'pipeB:out'
  []
[]

[Postprocessors]
  [mdot_A]
    type = SideAverageValue
    variable = rhouA
    boundary = 'pipeA:out'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [mdot_B]
    type = SideAverageValue
    variable = rhouA
    boundary = 'pipeB:in'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [p_A_out]
    type = SideAverageValue
    variable = p
    boundary = 'pipeA:out'
  []
  [p_B_in]
    type = SideAverageValue
    variable = p
    boundary = 'pipeB:in'
  []
  [offset]
    type = ParsedPostprocessor
    expression = 'p_A_out - p_B_in'
    pp_names = 'p_A_out p_B_in'
  []
  [p_bottom]
    type = SideAverageValue
    variable = p
    boundary = 'pipeB:out'
  []
  # pressure at the bottom minus p_top + rho * g * H, with the fixed rho of p_init_fn; this is only
  # a reference, because the density of the solution changes with depth
  [p_bottom_err]
    type = ParsedPostprocessor
    expression = 'p_bottom - (${p_top} + ${fparse -g} * 991.3 * ${H})'
    pp_names = 'p_bottom'
  []
[]

[VectorPostprocessors]
  # final p, rho and rho*u*A of each cell, sorted by elevation
  [prof]
    type = ElementValueSampler
    variable = 'p rho rhouA'
    block = 'pipeA pipeB'
    sort_by = z
    execute_on = FINAL
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  start_time = 0
  end_time = 300
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-2
    optimal_iterations = 12
    iteration_window = 3
    growth_factor = 1.5
    cutback_factor = 0.5
  []
  dtmin = 1e-6
  dtmax = 2.0
  solve_type = NEWTON
  line_search = basic
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-2
  nl_max_its = 30
  l_tol = 1e-4
  l_max_its = 100
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  [csv]
    type = CSV
    execute_vector_postprocessors_on = FINAL
  []
  [console]
    type = Console
    max_rows = 5
    execute_on = 'TIMESTEP_END'
  []
[]
