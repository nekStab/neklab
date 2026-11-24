"""Level 3 -- full non-regression tests.

Each test runs a case from scratch, exactly as a user would, and compares
the result with test/references/<reference>/values.json.

The reference values are produced by running, on the trusted version:

   ./neklabTests.py --levels full --update-references

which records the compared values instead of comparing them, and copies the
converged fields (saves_fields) that the fast tests start from.
"""
import json
import shutil

from lib.neklabFiles import read_eigenvalues
from lib.neklabTestCase import NeklabTestCase, say
from suites.fast import FastTest      # same machinery: inputs, compare, bounds


class FullTest(FastTest):
   """Base of the full tests of the regression cases (shared with fast).

   saves_fields  {file written by the run: name in references/<reference>/},
                 copied with --update-references
   """

   abstract = True
   level = "full"
   saves_references = True
   saves_fields = {}

   def _save_reference_values(self):
      """Only the compared values (not errors, start norms, ...), plus the fields."""
      values = self._reference_values()
      values.update({key: self.measured[key] for key in self.compare})
      path = self._reference_file()
      path.parent.mkdir(parents=True, exist_ok=True)
      path.write_text(json.dumps(values, indent=2, sort_keys=True) + "\n")
      say("refs", f"{', '.join(self.compare) or 'no values'} -> {path}")
      for produced, saved in self.saves_fields.items():
         if not (self.dir / produced).exists():
            raise FileNotFoundError(f"{produced} was not written by the run, cannot save "
                                    f"it as {path.parent / saved} (check saves_fields)")
         shutil.copy2(self.dir / produced, path.parent / saved)
         say("refs", f"{produced} -> {path.parent / saved}")


# =================================================================================
#  2D cylinder, steady: test/regression/cylinder_steady/ (shared with the fast tests)
# =================================================================================

class CylinderNewtonFull(FullTest):
   """Newton-Krylov for the steady base flow, from the initial condition (useric).

   With --update-references: writes bf_norm to values.json and the converged
   base flow as BF_ref.fld, which CylinderNewtonFast starts from.
   """
   case_dir = "regression/cylinder_steady_newton"
   case_name = "1cyl"
   settings_source = "1cyl.usr"
   reference = "cylinder/newton_fp"

   inputs = {"BF.fld": "cylinder/newton_fp/BF_ref.fld"}
   settings = {"if_load_BF": False, "noise_level": 0.0}
   checks = ["newton_fp.started_off_fixed_point"]
   compare = {"bf_norm": 1.0e-6}
   bounds = {"bf_error": (0.0, 1.0e-5)}
   saves_fields = {"BF_1cyl0.f00001": "BF_ref.fld"}