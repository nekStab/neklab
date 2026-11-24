"""Level 2 -- fast tests.

Each fast test runs the same case folder and .usr as its full test, but starts
next to the answer: a converged reference field + a little noise. It passes if
the solver gets back to the reference answer:

   * every NEKLAB_CHECK the run prints passes and the run reaches its end,
   * `compare`: the results agree with test/references/<reference>/values.json
     (written by the full test with --update-references) to `tol`,
   * `bounds`: the results lie in absolute bounds (e.g. the error of the
     converged field against the reference field).

   ./neklabTests.py --levels fast
   ./neklabTests.py --levels fast -k Newton


ADDING A FAST TEST
------------------
   class CylinderEigsDirFast(FastTest):
      \"\"\"What it starts from and what it checks.\"\"\"
      case_dir = "regression/cylinder_steady"
      case_name = "1cyl"
      reference = "cylinder/eigs_dir"
      inputs = {"BF.fld": "cylinder/newton_fp/BF_ref.fld"}
      settings = {...}                     # variables of the .usr, see below
      checks = ["eigs_dir.started_off_eigenvector"]
      compare = {"lambda_1_re": 1e-5}
      bounds = {"mode_overlap": (0.999, 1.001)}

`settings` edits initialised variables of the .usr (`real(dp) :: noise_level
= 1.0e-5_dp`), never `parameter`s; their names must be unique in the file.
"""
import shutil

from lib.neklabTestCase import NeklabTestCase, TEST_DIR, say

REFERENCES = TEST_DIR / "references"


class FastTest(NeklabTestCase):
   """Base of the fast tests: copy the reference inputs, run, check the results."""

   abstract = True
   level = "fast"

   # --- what a fast test sets ------------------------------------------------------
   inputs = {}          # {name in the scratch case: path under test/references/}
   checks = ()          # NEKLAB_CHECK names the run must report
   compare = {}         # {result: tol}, against references/<reference>/values.json
   bounds = {}          # {result: (low, high)}
   completed = "tests completed"     # printed by tests_completed() at the end of a run

   def test_run(self):
      self.copy_inputs()
      self.prepare()
      self.build()
      self.run_nek()
      self.check_driver(expected=self.checks, completed=self.completed)
      self.check_results()
      self.finish()

   def copy_inputs(self):
      """Put the reference fields into the scratch case, before the run."""
      for target, source in self.inputs.items():
         path = REFERENCES / source
         if path.exists():
            shutil.copy2(path, self.dir / target)
            say("input", f"{target} <- references/{source}")
         elif (self.dir / target).exists():
            say("input", f"{target}: references/{source} does not exist yet, "
                         f"using the file of {self.case_dir}/")
         else:
            raise FileNotFoundError(
               f"input {target}: neither references/{source} nor {self.case_dir}/{target} "
               f"exists. Run the full test with --update-references first.")

   def check_results(self):
      for key, tol in self.compare.items():
         value = self._result(key)
         if value is not None:
            self.check_reference(key, value, tol)
      for key, (low, high) in self.bounds.items():
         value = self._result(key)
         if value is not None:
            self.check(low <= value <= high, key, f"= {value!r} (bounds [{low}, {high}])")

   def _result(self, key):
      value = self.measured.get(key)
      if not isinstance(value, (int, float)):
         self.check(False, key, "not reported by the run (NEKLAB_RESULT line missing)")
         return None
      return value


# =================================================================================
#  2D cylinder, steady: test/regression/cylinder_steady/ (shared with the full tests)
# =================================================================================

class CylinderNewtonFast(FastTest):
   """Newton-Krylov for the steady base flow, from the reference base flow + noise.

   The .usr loads BF.fld, adds noise_level*||BF|| of random noise, reconverges,
   and loads BF.fld again to print the error of the result against it.
   """
   case_dir = "regression/cylinder_steady"
   case_name = "1cyl"
   settings_source = "1cyl.usr"
   reference = "cylinder/newton_fp"

   inputs = {"BF.fld": "cylinder/newton_fp/BF_ref.fld"}
   settings = {"if_load_BF": True, "noise_level": 5.0e-6}
   checks = ["newton_fp.started_off_fixed_point"]
   compare = {"bf_norm": 1.0e-5}
   bounds = {"bf_error": (0.0, 1.0e-5)}