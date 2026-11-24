"""What a regression analysis is built from (used by suites/regression/*.py).

An Analysis is one solver run on one case (Newton for the steady base flow,
direct eigenvalues, ...). It is defined once and make_tests() turns it into
two tests that run the SAME .usr and differ only in their settings:

   full   from scratch, as a user would; compares with references/<reference>/
          values.json and, with --update-references, writes the answer and the
          converged fields there.
   fast   starts next to the answer (a reference field + a little noise); must
          get back to the same answer (looser tolerance) in few iterations.

See suites/regression/__init__.py for how to add one.
"""
import json
import shutil
import unittest
from dataclasses import dataclass, field

from lib.neklabTestCase import NeklabTestCase, TEST_DIR, say
from lib.neklabUnit import Case, merged

__all__ = ["Analysis", "Case", "RegressionTest", "make_tests"]

REFERENCES = TEST_DIR / "references"
VARIANTS = ("fast", "full")


@dataclass(frozen=True)
class Analysis:
   """One analysis of a case folder, and how its fast and full tests differ.

   case             Case (folder, name, .par/settings shared by its analyses)
   name             class name part: "NewtonFp" -> <prefix>NewtonFpFast / ...Full
   reference        folder under test/references/ for values.json and fields
   group            value of case.group_key in the .par (None: not set)
   settings_source  file whose settings `settings` edits [<case.name>.usr]
   inputs           {name in scratch: path under references/}, copied before the
                    run in both variants (if the reference does not exist yet,
                    the file of the case folder is used)
   saves_fields     {file written by the run: name in references/<reference>/},
                    copied by full --update-references
   checks           NEKLAB_CHECK names the run must report
   compare          {result: tol}, compared with values.json (fast: tol * tol_factor)
   bounds           {result: (low, high)}, absolute, both variants
   needs            analyses whose references this one reads (sets the order of
                    --update-references)
   fast, full       per-variant changes: settings, par, size, checks, skip,
                    and for fast also tol_factor [10] and
                    max {result: upper bound}  (e.g. iteration counts)
   """

   case: Case
   name: str
   reference: str
   group: int = None
   settings_source: str = None
   inputs: dict = field(default_factory=dict)
   saves_fields: dict = field(default_factory=dict)
   checks: tuple = ()
   compare: dict = field(default_factory=dict)
   bounds: dict = field(default_factory=dict)
   needs: tuple = ()
   fast: dict = field(default_factory=dict)
   full: dict = field(default_factory=dict)

   def variant(self, which):
      return getattr(self, which)

   @property
   def depth(self):
      """0 for an analysis that needs nothing, 1 + the deepest need otherwise."""
      return 1 + max(a.depth for a in self.needs) if self.needs else 0


class RegressionTest(NeklabTestCase):
   """Base of the generated fast/full tests."""

   abstract = True
   saves_references = True

   analysis = None        # Analysis
   variant = None         # "fast" or "full"
   completed = "tests completed"     # printed by tests_completed() at the end of a run

   # --- shortcuts --------------------------------------------------------------
   @property
   def opts(self):
      return self.analysis.variant(self.variant)

   # --- the test -----------------------------------------------------------------
   def test_analysis(self):
      self.prepare()
      self.build()
      self.run_nek()
      a, opts = self.analysis, self.opts
      self.check_driver(expected=opts.get("checks", a.checks), completed=self.completed)
      self.check_results()
      self.finish()

   def prepare(self):
      a, opts, case = self.analysis, self.opts, self.analysis.case
      group = {"GENERAL": {case.group_key: a.group}} if a.group is not None else {}
      self.par_values = merged(case.par, group, opts.get("par", {}))
      self.settings = merged(case.settings, opts.get("settings", {}))
      self.size_params = opts.get("size", {})
      self.copy_inputs()
      super().prepare()

   def copy_inputs(self):
      for target, source in self.analysis.inputs.items():
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
               f"exists. Run the full tests of {[n.name for n in self.analysis.needs]} "
               f"with --update-references first.")

   def check_results(self):
      a, opts = self.analysis, self.opts
      factor = opts.get("tol_factor", 10.0) if self.variant == "fast" else 1.0

      def value(key):
         v = self.measured.get(key)
         if not isinstance(v, (int, float)):
            self.check(False, key, "not reported by the run (NEKLAB_RESULT line missing)")
            return None
         return v

      for key, tol in a.compare.items():
         v = value(key)
         if v is not None:
            self.check_reference(key, v, float(f"{tol * factor:.6g}"))
      for key, (low, high) in a.bounds.items():
         v = value(key)
         if v is not None:
            self.check(low <= v <= high, key, f"= {v!r} (bounds [{low}, {high}])")
      for key, high in opts.get("max", {}).items():
         v = value(key)
         if v is not None:
            self.check(v <= high, key, f"= {v!r} (at most {high})")

   # --- references ---------------------------------------------------------------
   def _save_reference_values(self):
      """full: the compared values and the fields; fast: its `max` values, as fast.<key>.

      The rest of measured_values.json (errors, norms of the start, ...) is not
      a reference and is not saved.
      """
      a, opts = self.analysis, self.opts
      if self.variant == "full":
         new = {key: self.measured[key] for key in a.compare}
      else:
         new = {f"fast.{key}": self.measured[key] for key in opts.get("max", {})}
      values = self._reference_values()
      values.update(new)
      path = self._reference_file()
      path.parent.mkdir(parents=True, exist_ok=True)
      path.write_text(json.dumps(values, indent=2, sort_keys=True, default=str) + "\n")
      say("refs", f"{', '.join(new) or 'nothing'} -> {path.relative_to(TEST_DIR)}")

      if self.variant == "full":
         for produced, saved in a.saves_fields.items():
            src = self.dir / produced
            if not src.exists():
               self.check(False, produced, f"not written by the run, cannot save it "
                                           f"as references/{a.reference}/{saved}")
               continue
            shutil.copy2(src, path.parent / saved)
            say("refs", f"{produced} -> references/{a.reference}/{saved}")


def make_tests(*analyses, module, prefix=""):
   """Two RegressionTest classes per analysis: {name: class}, for globals().update().

   module  __name__ of the suite file (the runner only picks up classes that
           belong to the file it imports)
   """
   tests = {}
   for a in analyses:
      for variant in VARIANTS:
         opts = a.variant(variant)
         name = f"{prefix}{a.name}{variant.capitalize()}"
         attrs = dict(
            __module__=module,
            __doc__=f"{a.name} ({variant}), references/{a.reference}",
            analysis=a, variant=variant, level=variant, depth=a.depth,
            case_dir=a.case.folder, case_name=a.case.name,
            settings_source=a.settings_source or f"{a.case.name}.usr",
            reference=a.reference,
         )
         cls = type(name, (RegressionTest,), attrs)
         if opts.get("skip"):
            cls.test_analysis = unittest.skip(opts["skip"])(RegressionTest.test_analysis)
         tests[name] = cls
   return tests