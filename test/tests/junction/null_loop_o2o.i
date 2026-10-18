# Closed loop of water at rest with two JunctionOneToOne1Phase
#
# The loop of null_loop.i (two vertical legs of 7 m, joined at the top and at the bottom) with the
# official JunctionOneToOne1Phase, which has no volume and reconstructs the state across the
# junction.
#
# Notes:
# - the fluid is at rest under gravity, so any flow is a numerical error;
# - the slope reconstruction across the junction (ADJunctionOneToOne1PhaseUserObject::finalize(),
# through SlopeReconstruction1DInterface::getBoundaryElementSlopes()) places the cell on the other
# side at the mirrored position of the interior neighbor, not at its real position;
# - this application has no corrected version of JunctionOneToOne1Phase;
# - JunctionOneToOne1Phase needs the same area in the two channels, so there is no AB parameter;
# - nA and nB are the numbers of cells of each leg (dz = 7 m / n);
# - g is the vertical component of gravity [m/s^2], and rdg is the slope reconstruction (minmod or
# none);
# - mdot_A is rho*u*A at the top of leg A (positive: up through leg A and down through leg B).

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
  # hydrostatic initial pressure with rho = 991.3 kg/m^3 (this stiffened gas at 313 K and 1 atm);
  # with g = 0 the pressure is uniform
  [p_init_fn]
    type = ParsedFunction
    expression = '${p_top} + ${fparse -g} * 991.3 * (${H} - z)'
  []
[]

[Components]
  [legA]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '0 0 1'
    length = ${H}
    n_elems = ${nA}
    A = ${AA}
    D_h = 0.8
  []
  [legB]
    type = FlowChannel1Phase
    position = '1 0 ${H}'
    orientation = '0 0 -1'
    length = ${H}
    n_elems = ${nB}
    A = ${AA}
    D_h = 0.8
  []
  [top_jct]
    type = JunctionOneToOne1Phase
    connections = 'legA:out legB:in'
  []
  [bot_jct]
    type = JunctionOneToOne1Phase
    connections = 'legB:out legA:in'
  []
[]

[Postprocessors]
  [mdot_A]
    type = SideAverageValue
    variable = rhouA
    boundary = 'legA:out'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [mdot_B]
    type = SideAverageValue
    variable = rhouA
    boundary = 'legB:out'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [p_A_top]
    type = SideAverageValue
    variable = p
    boundary = 'legA:out'
  []
  [p_B_top]
    type = SideAverageValue
    variable = p
    boundary = 'legB:in'
  []
  # pressure jump between the cell centers next to the top junction; the exact hydrostatic value
  # is -rho * g * (dz_B - dz_A) / 2
  [offset_top]
    type = ParsedPostprocessor
    expression = 'p_A_top - p_B_top'
    pp_names = 'p_A_top p_B_top'
  []
  [T_A_top]
    type = SideAverageValue
    variable = T
    boundary = 'legA:out'
  []
[]

[VectorPostprocessors]
  # final p, rho and rho*u*A of each cell of each leg, sorted by elevation
  [prof_A]
    type = ElementValueSampler
    variable = 'p rho rhouA'
    block = 'legA'
    sort_by = z
    execute_on = FINAL
  []
  [prof_B]
    type = ElementValueSampler
    variable = 'p rho rhouA'
    block = 'legB'
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
  # the tolerances of our reactor model; with nl_rel_tol = 1e-6 and nl_abs_tol = 1e-4 the Newton
  # iteration reaches the residual floor at the junctions and the run fails
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
