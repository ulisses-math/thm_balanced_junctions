#!/usr/bin/env python3
"""
Hydrostatic balance, cell by cell, of a column at rest (null_column.i or null_tee.i).

Notes:
- the script reads the final ElementValueSampler profile of the vertical channels (p, rho, rhouA,
sorted by z) and, if given, the profile of the horizontal channel of null_tee.i (sorted by x);
- the reference is the hydrostatic profile integrated from the top cell with the density of the
solution (trapezoidal rule between cell centers), so it does not depend on a fixed density:
p_ref(z_k) = p_top_cell + g * sum_{j > k} (rho_j + rho_{j+1}) / 2 * (z_{j+1} - z_j);
- for the horizontal channel the reference is the pressure at the junction elevation,
extrapolated from the first cell above the junction with its own density.

Usage: check_hydro_profile.py <vertical_prof.csv> [<horizontal_prof.csv>] [--zj 3.5] [--g 9.81]
"""
import argparse
import csv


def load(path, key):
    rows = list(csv.DictReader(open(path)))
    rows.sort(key=lambda r: float(r[key]))
    return [(float(r[key]), float(r["p"]), float(r["rho"]), float(r["rhouA"])) for r in rows]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vertical")
    ap.add_argument("horizontal", nargs="?")
    ap.add_argument("--zj", type=float, default=3.5, help="junction elevation [m]")
    ap.add_argument("--g", type=float, default=9.81, help="magnitude of gravity [m/s^2]")
    a = ap.parse_args()

    prof = load(a.vertical, "z")
    n = len(prof)
    ref = [0.0] * n
    ref[-1] = prof[-1][1]
    for k in range(n - 2, -1, -1):
        z0, _, r0, _ = prof[k]
        z1, _, r1, _ = prof[k + 1]
        ref[k] = ref[k + 1] + a.g * 0.5 * (r0 + r1) * (z1 - z0)
    dev = [prof[k][1] - ref[k] for k in range(n)]

    # first cell below and first cell above the junction
    kb = max(k for k in range(n) if prof[k][0] < a.zj)
    ka = kb + 1
    print(f"{n} cells; junction at z = {a.zj} m between cells {kb} and {ka}")
    print(f"  {'z [m]':>8s} {'dz [mm]':>8s} {'p - p_ref [Pa]':>15s} {'rhouA [kg/s]':>13s}")
    for k in range(max(kb - 2, 0), min(ka + 3, n)):
        # cell size from the neighbor on the same side of the junction
        same = (k + 1 if k < n - 1 else k - 1) if k >= ka else (k - 1 if k > 0 else k + 1)
        dz = abs(prof[k][0] - prof[same][0])
        tag = "  <- above junction" if k == ka else ("  <- below junction" if k == kb else "")
        print(f"  {prof[k][0]:8.4f} {dz * 1e3:8.1f} {dev[k]:15.3f} {prof[k][3]:13.6f}{tag}")
    print(f"  max |p - p_ref| in the vertical channels = {max(abs(d) for d in dev):.3f} Pa")
    print(f"  max |rhouA| in the vertical channels     = {max(abs(r[3]) for r in prof):.3e} kg/s")

    if a.horizontal:
        hc = load(a.horizontal, "x")
        za, pa, ra, _ = prof[ka]
        p_j = pa + a.g * ra * (za - a.zj)
        devs = [p - p_j for _, p, _, _ in hc]
        print(f"  horizontal channel, {len(hc)} cells: max |p - p_ref| = "
              f"{max(abs(d) for d in devs):.3f} Pa, max |rhouA| = "
              f"{max(abs(r[3]) for r in hc):.3e} kg/s")


if __name__ == "__main__":
    main()
