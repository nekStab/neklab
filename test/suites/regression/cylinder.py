"""Regression analyses of the 2D cylinder (fast and full levels).

One case folder per analysis, each with its own .usr:

   test/regression/cylinder_steady_newton/     NewtonFp   Newton for the steady base flow
   test/regression/cylinder_steady_stability/  Eigs       direct + adjoint eigenpairs of it

Each Analysis below becomes Cylinder<name>Fast and Cylinder<name>Full.
"""
from lib.neklabRegression import Analysis, Case, make_tests
from lib.neklabUnit import merged

# --- Newton for the steady base flow ---------------------------------------------
NEWTON_FP = Analysis(
   case=Case(folder="regression/cylinder_steady_newton", name="1cyl"),
   name="NewtonFp",
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


# --- Eigenpairs of the base flow (direct and adjoint, one run) -------------------
EIGS_REF = "cylinder/eigs"
NEV = 4                    # leading modes checked by the full test


def eigenvalues(nev, tol):
   """{result: tol} for Re and |Im| of the first nev eigenvalues, both sides."""
   return {f"{side}_lambda_{i}_{part}": tol
           for side in ("dir", "adj") for i in range(1, nev + 1) for part in ("re", "abs_im")}


def mode_errors(nev, high):
   """{result: high} for the Re and Im errors of the first nev modes, both sides."""
   return {f"{side}_mode_{i}_{part}_err": high
           for side in ("dir", "adj") for i in range(1, nev + 1) for part in ("re", "im")}


def mode_files(side, nev):
   """{file the run writes: name in references/}. The .usr outposts Re then Im
   of each mode with prefix `side`: <side>1cyl0.f00001 = Re(mode 1), f00002 = Im(mode 1), ..."""
   files = {}
   for i in range(1, nev + 1):
      files[f"{side}1cyl0.f{2 * i - 1:05d}"] = f"{side}{i}_re.fld"
      files[f"{side}1cyl0.f{2 * i:05d}"] = f"{side}{i}_im.fld"
   return files


MODE_FILES = merged(mode_files("dir", NEV), mode_files("adj", NEV))

EIGS = Analysis(
   case=Case(folder="regression/cylinder_steady_stability", name="1cyl"),
   name="Eigs",
   reference=EIGS_REF,
   needs=(NEWTON_FP,),
   inputs={"BF.fld": "cylinder/newton_fp/BF_ref.fld"},
   saves_fields=MODE_FILES,
   # the .usr compares mode i with ref_<side><i>_re.fld / _im.fld when present
   reference_fields={f"ref_{name}": f"{EIGS_REF}/{name}" for name in MODE_FILES.values()},
   compare=eigenvalues(NEV, tol=1e-6),
   field_errors=mode_errors(NEV, high=1e-4),
   full=dict(
      settings={"nev": NEV, "kdim": 64, "load_x0": False},
   ),
   fast=dict(
      # start from the leading reference modes + noise, converge only those
      settings={"nev": 1, "kdim": 16, "load_x0": True},
      inputs={"x0r.fld": f"{EIGS_REF}/dir1_re.fld", "x0i.fld": f"{EIGS_REF}/dir1_im.fld",
              "y0r.fld": f"{EIGS_REF}/adj1_re.fld", "y0i.fld": f"{EIGS_REF}/adj1_im.fld"},
      compare=eigenvalues(1, tol=1e-6),          # x tol_factor (10)
      field_errors=mode_errors(1, high=1e-3),
   ),
)


globals().update(make_tests(NEWTON_FP, EIGS, module=__name__, prefix="Cylinder"))