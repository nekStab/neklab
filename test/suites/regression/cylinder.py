"""Regression analyses of the 2D cylinder (fast and full levels).

   test/regression/cylinder_steady/     newton_fp  (later: eigs_dir, eigs_adj)
   test/regression/cylinder_periodic/   later: newton_upo, floquet_dir, floquet_adj

Each Analysis below becomes Cylinder<name>Fast and Cylinder<name>Full.
"""
from lib.neklabRegression import Analysis, Case, make_tests

CYLINDER_STEADY = Case(
   folder="regression/cylinder_steady",
   name="1cyl",
)


NEWTON_FP = Analysis(
   case=CYLINDER_STEADY,
   name="NewtonFixedPoinshashareshar",
   reference="cylinder/newton_fp",
   # BF.fld is read twice by the .usr: as the start (fast) and to measure the
   # error of the converged base flow (both variants)
   inputs={"BF.fld": "cylinder/newton_fp/BF_ref.fld"},
   # ADAPT: the name outpost_dnek(bf, "BF_") gives the file (check after one run)
   saves_fields={"BF_1cyl0.f00001": "BF_ref.fld"},
   checks=["newton_fp.started_off_fixed_point"],
   compare={"bf_norm": 1e-6},
   bounds={"bf_rel_error": (0.0, 1e-5)},
   fast=dict(
      # the defaults of the .usr: load BF.fld, add noise
      # max={"newton_iterations": 4},    # once the .usr reports it
   ),
   full=dict(
      settings={"if_load_BF": False, "noise_level": 0.0},
   ),
)


# Next, in the same folder (group = userParam01 once the .usr switches on it):
# EIGS_DIR = Analysis(case=CYLINDER_STEADY, name="EigsDir", group=2,
#                     reference="cylinder/eigs_dir", needs=(NEWTON_FP,),
#                     inputs={"BF.fld": "cylinder/newton_fp/BF_ref.fld",
#                             "mode1.fld": "cylinder/eigs_dir/mode1_ref.fld"}, ...)


globals().update(make_tests(NEWTON_FP, module=__name__, prefix="Cylinder"))