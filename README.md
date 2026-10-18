thm_balanced_junctions
======================

A small MOOSE-based application with 0-D junctions for the `thermal_hydraulics` module that keep a
fluid at rest under gravity in hydrostatic balance, together with minimal reproducers and regression
tests.

We wrote these classes during a validation study of the `thermal_hydraulics` module against
loss-of-flow data of the IEA-R1 pool-type research reactor. The framework and the module are not
modified: the classes derive from the official components or follow them closely.

## The issue in short

- The channel state that `VolumeJunction1Phase` and `JunctionParallelChannels1Phase` receive from a
connected `FlowChannel1Phase` is the cell average of the adjacent cell
(`ADVolumeJunction1PhaseUserObject.C`, lines 104-110), that is, the pressure at the cell center,
half a cell above or below the junction.
- The junction itself has no elevation and no gravity term, so half a cell of hydrostatic head is
left without a counterpart.
- In a fluid at rest, each cell next to a junction gets a steady pressure error of rho * g * dz / 4,
where dz is the height of that cell (with the default rDG reconstruction; with
`rdg_slope_reconstruction = none` the error is twice as large).
- The 0-D momentum equation of the official junctions absorbs the imbalance with a small standing
velocity in the cells next to the junction, so the error is a fixed pressure offset and not a
circulation. It only matters when the physical driving pressures are of the same order, as in
buoyancy-driven flow reversal.
- `JunctionOneToOne1Phase` is not balanced either, for a different reason: the slope reconstruction
places the cell on the other side of the junction at the mirrored position of the interior neighbor
(`SlopeReconstruction1DInterface.h`, "assumed to be the same distance away"), not at the real
distance between the cell centers. With unequal cells, the slope across the junction is too large in
the smaller cell and too small in the larger one, and the minmod limiter keeps the smaller slope.
The reconstruction is done in `ADJunctionOneToOne1PhaseUserObject::finalize()`, through
`getBoundaryElementSlopes()`; a corrected version would pass the real position of the cell on the
other side. This application has reproducers for it, but no corrected version.
- The documentation of `JunctionOneToOne1Phase` says that the two channels must have the same flow
area, but `check()` only checks the number of connections, so channels with different areas are
accepted without a warning.

## The correction

Before the channel-side Riemann state is formed, the cell-center state is extrapolated
hydrostatically to the junction face with the density of the cell:

- p_face = p_c + rho_c * g . (x_face - x_c);
- rho and u are kept, and e_face = e_c + (p_face - p_c) / (dp/de)_v.

Notes:

- this is exact for the stiffened gas equation of state and first order for a general equation of
state, and it needs no inverse equation of state call;
- it is the junction-face counterpart of the hydrostatic reconstruction of Audusse et al. (2004),
SIAM J. Sci. Comput. 25(6), 2050-2065;
- the function is in `include/utils/HydrostaticFace.h`.

## Contents

| Object | Derived from | Changes |
| :- | :- | :- |
| `HydrostaticVolumeJunction1Phase` | `VolumeJunction1Phase` | face state only |
| `HydrostaticJunctionParallelChannels1Phase` | `JunctionParallelChannels1Phase` | face state; wall pressure of the area change term at the face; one held flow direction for the inlet/outlet classification |
| `StagnationJunction1Phase` | `FlowJunction1Phase` | a junction with mass and energy only (no 0-D momentum), with the face state |

Notes:

- each class has the parameter `hydrostatic_correction` (default `true`); with `false`,
`HydrostaticVolumeJunction1Phase` gives the official results to all digits in the reproducers,
and `HydrostaticJunctionParallelChannels1Phase` keeps only change 3;
- the changes are marked in the source with `change 1`, `change 2` and `change 3`, and the class
documentation in the headers explains each one;
- the correction assumes that the channel ends connected to a junction are at the same elevation
(the junction has one pressure for all faces), and the components check it;
- `StagnationJunction1Phase` is the junction that we use in our reactor model; it also shows the
imbalance without the acoustic brake of the 0-D momentum equation.

Reproducers, in `test/tests/junction` (all at rest, isothermal, stiffened gas water at 313 K):

- `null_column.i`: a 7 m column, open at the top and closed at the bottom, with a junction at
mid-height (channels on opposite sides of the junction);
- `null_tee.i`: the same column with a third, horizontal channel at the junction (opposite sides and
a vertical to horizontal connection);
- `null_loop.i`: a closed loop of two 7 m vertical legs and two junctions (both legs on the same side
of each junction);
- `null_column_o2o.i` and `null_loop_o2o.i`: the column and the loop with the official
`JunctionOneToOne1Phase`, which takes different parameters (no volume, the same area on both sides).

The junction type and the meshes are set on the command line, for example:

```
./thm_balanced_junctions-opt -i test/tests/junction/null_column.i jt=VolumeJunction1Phase nA=100 nB=25
./thm_balanced_junctions-opt -i test/tests/junction/null_column.i jt=HydrostaticVolumeJunction1Phase nA=100 nB=25
```

## Build and test

Notes:

- we use MOOSE commit `13877661c5` (2026-05-08) and the `moose-dev` conda environment of that
commit;
- the `Makefile` looks for MOOSE in `../moose` (or set `MOOSE_DIR`).

```
make -j 8
./run_tests -j 8
scripts/reproduce.sh
```

`run_tests` runs seven CSVDiff tests (absolute tolerance 1e-4 kg/s on the mass flow rate) with the
three classes. `scripts/reproduce.sh` runs the reproducers with the official junctions and with the
classes of this application (14 runs), and prints the comparison (about one minute on 16 cores).
`scripts/check_hydro_profile.py` compares a final profile with the hydrostatic profile, cell by cell,
using the density of the solution.

## Results

Output of `scripts/reproduce.sh`, 35 mm cells on one side of the junction and 140 mm on the other,
at t = 300 s.

| Case | Official junction | This application |
| :- | :- | :- |
| column, `JunctionParallelChannels1Phase` | +85.1 Pa (cell above), -340.6 Pa (cell below) | 0.000 Pa in every cell |
| column, `VolumeJunction1Phase` | +85.1 Pa, -340.6 Pa | 0.000 Pa in every cell |
| column with a 4.6:1 contraction, `JunctionParallelChannels1Phase` | +85.1 Pa (cell above), -34.1 Pa (cell below), +306.5 Pa in the rest of the lower channel | 0.000 Pa in every cell |
| T junction, `VolumeJunction1Phase` | +85.1 Pa, -340.6 Pa, -85.1 Pa in the horizontal channel | 0.000 Pa in every cell of the three channels |
| loop, `VolumeJunction1Phase` | jump between the cell centers next to the top junction -255.4 Pa (exact: about -510.6 Pa); rho*u*A -0.032 kg/s | -510.86 Pa; below 1e-5 kg/s |
| loop, `StagnationJunction1Phase` | with `hydrostatic_correction = false`: loop flow of 1214 kg/s | below 1e-5 kg/s |
| column, `JunctionOneToOne1Phase` | +21.3 Pa (cell above), -242.7 Pa (cell below) | no corrected version |
| loop, `JunctionOneToOne1Phase` | jump between the cell centers next to the top junction -306.5 Pa (exact: about -510.6 Pa); rho*u*A -0.034 kg/s | no corrected version |

Notes:

- the deviations are p - p_ref, with p_ref integrated from the top cell with the density of the
solution;
- with the official junctions, the cells next to the junction keep a standing rho*u*A of 0.03 to
0.13 kg/s; it is not a flow through the column, which is closed at the bottom;
- the column error follows the law rho * g * dz / 4: 9725 * 0.035 / 4 = 85 Pa and
9725 * 0.140 / 4 = 340 Pa;
- in the loop, the error of the jump at the top junction is rho * g * (0.140 - 0.035) / 4 = 255 Pa;
- the T junction gives the same result on 1 and 2 MPI ranks.

## Status

Notes:

- this is work in progress, shared with the MOOSE developers before we submit a paper on it;
- the inputs of our reactor model are available on request;
- we would be glad to contribute the face correction to the `thermal_hydraulics` module if it is
useful.

## License

LGPL 2.1, the same as MOOSE (see `LICENSE`). Parts of the source are adapted from the
`thermal_hydraulics` module of MOOSE (Battelle Energy Alliance, LLC), and each file says which ones.

## Contact

Ulisses Alves Maciel Neto (Faculdade Ibmec São Paulo and Centro Tecnológico da Marinha em São Paulo,
CTMSP), ulisses.neto@professores.ibmec.edu.br.
