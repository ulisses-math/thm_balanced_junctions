#!/bin/bash
# Runs the three reproducers with the official junctions and with the junctions of this
# application, and prints the comparison.
#
# Notes:
# - run it from the root of the repository, after make;
# - the runs go in parallel (about 1 minute each) and the results go to reproduce_out/;
# - set MOOSE_EXEC to use another executable.

cd "$(dirname "$0")/.." || exit 1
EXEC=${MOOSE_EXEC:-./thm_balanced_junctions-opt}
D=test/tests/junction
O=reproduce_out
mkdir -p $O

run() # run <tag> <input> <args...>
{
  local tag=$1 inp=$2
  shift 2
  $EXEC -i $D/$inp "$@" Outputs/file_base=$O/$tag > $O/$tag.log 2>&1
  echo "$tag: exit code $?"
}

# column (channels on opposite sides of the junction), 35 mm cells above and 140 mm below
run col_jpc null_column.i jt=JunctionParallelChannels1Phase nA=100 nB=25 &
run col_hjpc null_column.i jt=HydrostaticJunctionParallelChannels1Phase nA=100 nB=25 &
run col_vj null_column.i jt=VolumeJunction1Phase nA=100 nB=25 &
run col_hvj null_column.i jt=HydrostaticVolumeJunction1Phase nA=100 nB=25 &
# column with a 4.6:1 area contraction at the junction
run colc_jpc null_column.i jt=JunctionParallelChannels1Phase nA=100 nB=25 AB=0.1087 &
run colc_hjpc null_column.i jt=HydrostaticJunctionParallelChannels1Phase nA=100 nB=25 AB=0.1087 &
# T junction (opposite sides and a horizontal channel)
run tee_vj null_tee.i jt=VolumeJunction1Phase nA=100 nB=25 &
run tee_hvj null_tee.i jt=HydrostaticVolumeJunction1Phase nA=100 nB=25 &
# closed loop (both legs on the same side of each junction), 35 mm cells in leg A, 140 mm in leg B
run loop_vj null_loop.i jt=VolumeJunction1Phase nA=200 nB=50 &
run loop_hvj null_loop.i jt=HydrostaticVolumeJunction1Phase nA=200 nB=50 &
run loop_stag_off null_loop.i jt=StagnationJunction1Phase nA=200 nB=50 \
  Components/top_jct/hydrostatic_correction=false Components/bot_jct/hydrostatic_correction=false &
run loop_stag null_loop.i jt=StagnationJunction1Phase nA=200 nB=50 &
# official JunctionOneToOne1Phase (no corrected version in this application)
run col_o2o null_column_o2o.i nA=100 nB=25 &
run loop_o2o null_loop_o2o.i nA=200 nB=50 &
wait

for tag in col_jpc col_hjpc col_vj col_hvj colc_jpc colc_hjpc col_o2o; do
  echo
  echo "== $tag"
  python3 scripts/check_hydro_profile.py $O/${tag}_prof_[0-9]*.csv | tail -7
done
for tag in tee_vj tee_hvj; do
  echo
  echo "== $tag"
  python3 scripts/check_hydro_profile.py $O/${tag}_prof_[0-9]*.csv $O/${tag}_prof_C_*.csv | tail -8
done
echo
python3 - $O <<'EOF'
import csv
import sys

import glob

o = sys.argv[1]


def top_cell(path):
    rows = sorted(csv.DictReader(open(glob.glob(path)[0])), key=lambda r: float(r["z"]))
    return float(rows[-1]["z"]), float(rows[-1]["rho"])


print("== closed loop, at t = 300 s")
print(f"  {'case':14s} {'mdot_A [kg/s]':>14s} {'offset_top [Pa]':>16s} {'exact [Pa]':>11s}")
for tag in ["loop_vj", "loop_hvj", "loop_stag_off", "loop_stag", "loop_o2o"]:
    r = list(csv.DictReader(open(f"{o}/{tag}.csv")))[-1]
    # exact jump between the two cell centers next to the top junction, with their own density
    zA, rA = top_cell(f"{o}/{tag}_prof_A_*.csv")
    zB, rB = top_cell(f"{o}/{tag}_prof_B_*.csv")
    exact = -9.81 * 0.5 * (rA + rB) * (zA - zB)
    print(f"  {tag:14s} {float(r['mdot_A']):14.5f} {float(r['offset_top']):16.2f} {exact:11.2f}")
EOF
