"""What a regression analysis is built from (used by suites/regression/*.py).

An Analysis is one solver run on one case (Newton for the steady base flow,
direct eigenvalues, ...). It is defined once and make_tests() turns it into
two tests that run the SAME .usr and differ only in their settings:

   full   from scratch, as a user would; compares with references/<reference>/
          values.json and, with --update-references, writes the answer and the
          converged fields there.
   fast   starts next to the answer (a reference field + a little noise); must
          get back to the same answer (looser tolerance).

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

# options of Analysis that fast/full may override (inputs and reference_fields
# are merged on top of the shared ones, the others replace them)
MERGED_OPTIONS = ("inputs", "reference_fields")
REPLACED_OPTIONS = ("checks", "compare", "bounds", "field_errors")


@dataclass(frozen=True)
class Analysis:
   """One analysis of a case folder, and how its fast and full tests differ.

   case              Case (folder, name, .par/settings shared by its analyses)
   name              class name part: "NewtonFp" -> <prefix>NewtonFpFast / ...Full
   reference         folder under test/references/ for values.json and fields
   group             value of case.group_key in the .par (None: not set)
   settings_source   file whose settings `settings` edits [<case.name>.usr]
   inputs            {name in scratch: path under references/}, REQUIRED inputs,
                     copied before the run. If the reference does not exist yet,
                     the file of the case folder is used, except with
                     --update-references: then the test fails (the references
                     must not be built on anything but references).
   reference_fields  {name in scratch: path under references/}, OPTIONAL inputs:
                     copied when they exist (never with --update-references).
                     The .usr compares its results with them when present and
                     reports the errors, checked with field_errors.
   saves_fields      {file written by the run: name in references/<reference>/},
                     copied by full --update-references
   checks            NEKLAB_CHECK names the run must report
   compare           {result: tol}, compared with values.json (fast: tol * tol_factor)
   bounds            {result: (low, high)}, absolute
   field_errors      {result: upper bound}, errors against reference_fields;
                     checked only when all the reference fields were copied
   needs             analyses whose references this one reads (sets the order:
                     an analysis runs after the ones it needs)
   fast, full        per-variant changes: settings, par, size, skip, and any of
                     inputs, reference_fields (merged on top of the shared ones),
                     checks, compare, bounds, field_errors (replace them);
                     fast also tol_factor [10] and
                     max {result: upper bound}  (e.g. iteration counts)
   """

   case: Case
   name: str
   reference: str
   group: int = None
   settings_source: str = None
   inputs: dict = field(default_factory=dict)
   reference_fields: dict = field(default_factory=dict)
   saves_fields: dict = field(default_factory=dict)
   checks: tuple = ()
   compare: dict = field(default_factory=dict)
   bounds: dict = field(default_factory=dict)
   field_errors: dict = field(default_factory=dict)
   needs: tuple = ()
   fast: dict = field(default_factory=dict)
   full: dict = field(default_factory=dict)

   def __post_init__(self):
      known = set(MERGED_OPTIONS + REPLACED_OPTIONS) | {"settings", "par", "size", "skip"}
      for which in VARIANTS:
         extra = set(self.variant(which)) - known
         if which == "fast":
            extra -= {"tol_factor", "max"}
         if extra:
            raise TypeError(f"Analysis {self.name!r}: unknown {which} option(s) {sorted(extra)}")

   def variant(self, which):
      return getattr(self, which)

   def option(self, name, which):
      """The value of option `name` for variant `which` (fast/full overrides applied)."""
      opts = self.variant(which)
      if name in MERGED_OPTIONS:
         return merged(getattr(self, name), opts.get(name, {}))
      return opts.get(name, getattr(self, name))

   @property
   def depth(self):
      """0 for an analysis that needs nothing, 1 + the deepest need otherwise."""
      return 1 + max(a.depth for a in self.needs) if self.needs else 0


class RegressionTest(NeklabTestCase):
   """Base of the generated fast/full tests."""

   abstract = True

   analysis = None        # Analysis
   variant = None         # "fast" or "full"
   completed = "tests completed"     # printed by tests_completed() at the end of a run

   # --- shortcuts --------------------------------------------------------------
   @property
   def opts(self):
      return self.analysis.variant(self.variant)

   def option(self, name):
      return self.analysis.option(name, self.variant)

   # --- the test -----------------------------------------------------------------
   def test_analysis(self):
      self.prepare()
      self.build()
      self.run_nek()
      self.check_driver(expected=self.option("checks"), completed=self.completed)
      self.check_results()
      self.finish()

   def prepare(self):
      a, opts, case = self.analysis, self.opts, self.analysis.case
      group = {"GENERAL": {case.group_key: a.group}} if a.group is not None else {}
      self.par_values = merged(case.par, group, opts.get("par", {}))
      self.settings = merged(case.settings, opts.get("settings", {}))
      self.size_params = opts.get("size", {})
      self.copy_inputs()
      self.copy_reference_fields()
      super().prepare()

   def _producers(self, source):
      """Which full tests write references/<source>: for the error messages."""
      a = self.analysis
      producers = [f"{n.name} (full)" for n in a.needs]
      if source.startswith(a.reference.rstrip("/") + "/"):
         producers.append(f"{a.name} (full)")
      return producers or ["the analysis that writes it"]

   def copy_inputs(self):
      for target, source in self.option("inputs").items():
         path = REFERENCES / source
         if path.exists():
            shutil.copy2(path, self.dir / target)
            say("input", f"{target} <- references/{source}")
         elif self.update_references:
            raise FileNotFoundError(
               f"input {target}: references/{source} does not exist. With "
               f"--update-references every input must come from the references, "
               f"written earlier in the same run by {self._producers(source)}: "
               f"see why that test failed.")
         elif (self.dir / target).exists():
            say("input", f"{target}: references/{source} does not exist yet, "
                         f"using the file of {self.case_dir}/")
         else:
            raise FileNotFoundError(
               f"input {target}: neither references/{source} nor {self.case_dir}/{target} "
               f"exists. Run {self._producers(source)} with --update-references first.")

   def copy_reference_fields(self):
      """Copy the reference fields that exist; remember whether all of them did.

      Not with --update-references: the references are being rebuilt, and the
      old fields would only be compared with what replaces them.
      """
      wanted = self.option("reference_fields")
      self.have_reference_fields = False
      if not wanted:
         return
      if self.update_references:
         say("refs", f"{len(wanted)} reference fields not copied (--update-references): "
                     f"field errors not checked")
         return
      missing = []
      for target, source in wanted.items():
         path = REFERENCES / source
         if path.exists():
            shutil.copy2(path, self.dir / target)
         else:
            missing.append(source)
      self.have_reference_fields = not missing
      if missing:
         say("refs", f"{len(missing)}/{len(wanted)} reference fields missing "
                     f"(e.g. references/{missing[0]}): field errors not checked")
      else:
         say("refs", f"{len(wanted)} reference fields copied")

   def check_results(self):
      factor = self.opts.get("tol_factor", 10.0) if self.variant == "fast" else 1.0

      def value(key):
         v = self.measured.get(key)
         if not isinstance(v, (int, float)):
            self.check(False, key, "not reported by the run (NEKLAB_RESULT line missing)")
            return None
         return v

      for key, tol in self.option("compare").items():
         v = value(key)
         if v is not None:
            self.check_reference(key, v, float(f"{tol * factor:.6g}"))
      for key, (low, high) in self.option("bounds").items():
         v = value(key)
         if v is not None:
            self.check(low <= v <= high, key, f"err {v!r} (bounds [{low}, {high}])")
      if self.have_reference_fields:
         for key, high in self.option("field_errors").items():
            v = value(key)
            if v is not None:
               self.check(v <= high, key, f"err {v!r}  tol {high}")
      for key, high in self.opts.get("max", {}).items():
         v = value(key)
         if v is not None:
            self.check(v <= high, key, f"err {v!r}  tol {high}")

   # --- references (full only, see NeklabTestCase.finish) -------------------------
   def check_reference_outputs(self):
      """Before anything is written: every field to save must exist."""
      a = self.analysis
      for produced, saved in a.saves_fields.items():
         self.check((self.dir / produced).exists(), produced,
                    f"not written by the run, cannot save it as "
                    f"references/{a.reference}/{saved}", quiet=True)

   def _save_reference_values(self):
      """The compared values and the fields.

      The rest of measured_values.json (errors, norms of the start, ...) is not
      a reference and is not saved.
      """
      a = self.analysis
      new = {key: self.measured[key] for key in self.option("compare")}
      values = self._reference_values()
      values.update(new)
      path = self._reference_file()
      path.parent.mkdir(parents=True, exist_ok=True)
      path.write_text(json.dumps(values, indent=2, sort_keys=True, default=str) + "\n")
      say("refs", f"{', '.join(new) or 'nothing'} -> {path.relative_to(TEST_DIR)}")
      for produced, saved in a.saves_fields.items():
         shutil.copy2(self.dir / produced, path.parent / saved)
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
