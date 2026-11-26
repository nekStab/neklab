"""Regression analyses of the 2D cylinder (fast and full levels).

One case folder per analysis, each with its own .usr:

   test/regression/cylinder_steady_newton/     NewtonFp   Newton for the steady base flow
   test/regression/cylinder_steady_stability/  Eigs       direct + adjoint eigenpairs of it

Each Analysis below becomes Cylinder<name>Fast and Cylinder<name>Full.
"""
from lib.neklabRegression import Analysis, Case, make_tests
from lib.neklabUnit import merged

# --- Newton for the steady base flow ---------------------------------------------
NEWTON_REF = "cylinder/newton_fp"

NEWTON_FP = Analysis(
   case=Case(folder="regression/cylinder_steady_newton", name="1cyl"),
   name="NewtonFp",
   reference=NEWTON_REF,
   # the .usr compares the converged base flow with ref_BF.fld when it is there
   # and reports bf_error (not on the first --update-references: no reference yet)
   reference_fields={"ref_BF.fld": f"{NEWTON_REF}/BF_ref.fld"},
   saves_fields={"BF_1cyl0.f00001": "BF_ref.fld"},
   checks=["newton_fp.started_off_fixed_point"],
   compare={"bf_norm": 1e-6},
   field_errors={"bf_error": 1e-6},
   fast=dict(
      inputs={"BF.fld": f"{NEWTON_REF}/BF_ref.fld"},
   ),
   full=dict(
      settings={"if_load_BF": False, "noise_level": 0.0},     # start from useric
   ),
)


# --- Eigenpairs of the base flow (direct and adjoint, one run) -------------------
# The .usr is still a stub: it does not report its results yet, so both tests fail.
EIGS_REF = "cylinder/eigs"
NEV = 2                    # leading modes checked by the full test
MODES = ("dir", "adj")
PARTS = ("re", "im")


def eigenvalues(nev, tol):
   """{result: tol} for Re and |Im| of the first nev eigenvalues, both modes."""
   return {f"{mode}_lambda_{i}_{part}": tol
           for mode in MODES for i in range(1, nev + 1) for part in PARTS}


def mode_errors(nev, high):
   """{result: high} for the Re and Im errors of the first nev modes, both modes."""
   return {f"{mode}_mode_{i}_{part}_err": high
           for mode in ("dir", "adj") for i in range(1, nev + 1) for part in ("re", "im")}


def mode_files(mode, nev):
   """{file the run writes: name in references/}. The .usr outposts Re then Im
   of each mode with prefix `mode`: <mode>1cyl0.f00001 = Re(mode 1), f00002 = Im(mode 1), ..."""
   files = {}
   for i in range(1, nev + 1):
      files[f"{mode}1cyl0.f{2 * i - 1:05d}"] = f"{mode}{i}_re.fld"
      files[f"{mode}1cyl0.f{2 * i:05d}"] = f"{mode}{i}_im.fld"
   return files


MODE_FILES = merged(mode_files("dir", NEV), mode_files("adj", NEV))

EIGS = Analysis(
   case=Case(folder="regression/cylinder_steady_stability", name="1cyl"),
   name="Eigs",
   reference=EIGS_REF,
   needs=(NEWTON_FP,),
   inputs={"BF.fld": f"{NEWTON_REF}/BF_ref.fld"},
   saves_fields=MODE_FILES,
   # the .usr compares mode i with ref_<side><i>_re.fld / _im.fld when present
   reference_fields={f"ref_{name}": f"{EIGS_REF}/{name}" for name in MODE_FILES.values()},
   compare=eigenvalues(NEV, tol=1e-7),
   field_errors=mode_errors(NEV, high=1e-4),
   full=dict(
      settings={"nev": NEV, "kdim": 128, "load_x0": False},
   ),
   fast=dict(
      # start from the leading reference modes + noise, converge only those
      settings={"nev": 2, "kdim": 128, "load_x0": True},
      inputs={"x0.fld": f"{EIGS_REF}/dir1_re.fld",
              "y0.fld": f"{EIGS_REF}/adj1_re.fld"},
      compare=eigenvalues(1, tol=1e-6),          # x tol_factor (10)
      field_errors=mode_errors(1, high=1e-3),
   ),
)


globals().update(make_tests(NEWTON_FP, EIGS, module=__name__, prefix="Cylinder"))
