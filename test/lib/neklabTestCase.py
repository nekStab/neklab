"""Base class of all neklab tests.

A test is a complete case folder (SIZE, .par, .usr, mesh, input fields) plus
a few values to change in it. The test:

   1. copies the folder to test/scratch/<TestClass>/ (the original is never
      touched; symlinks are followed, so the copy is self-contained),
   2. changes size_params (SIZE), par_values (.par) and settings (initialised
      variables of the .usr) in that copy,
   3. prints all parameters, marking those that differ from the case folder,
   4. builds the case from clean with makeneklab and runs it,
   5. reads the NEKLAB_CHECK / NEKLAB_RESULT lines of the log and checks them.

Test definitions do not use this class directly but one of
   lib/neklabUnit.py        UnitTest        (level unit)
   lib/neklabRegression.py  Analysis        (levels fast and full)
"""
import sys
import json
import os
import shutil
import time
import unittest
from fnmatch import fnmatch
from pathlib import Path

from lib import neklabBuildRun as nek
from lib import neklabFiles as files

TEST_DIR = Path(__file__).resolve().parent.parent      # the test/ folder
LEVELS = ("unit", "fast", "full")
NAME_WIDTH = 60                         # width of check and result names

# Never copied into the scratch case: build products and outputs of earlier runs
IGNORE = ("obj", "nek5000", "makefile", "makefile_usr.inc", "*.log.*", "*.out",
          "SESSION.NAME", "ioinfo", "eigs_output.txt", "nek_dummy.usr",
          ".neklab_stash", "*.mod", "*.smod")


if sys.stdout.isatty():
   GREEN, RED, RESET = "\033[32m", "\033[31m", "\033[0m"
else:
   GREEN = RED = RESET = ""


def status(ok):
   return f"{GREEN}PASSED{RESET}" if ok else f"{RED}FAILED{RESET}"


def say(step, message=""):
   """One line of the runner's narration: '[step]    message'."""
   print(f"  {'[' + step + ']':10} {message}", flush=True)


class NeklabTestCase(unittest.TestCase):

   # --- what a test class sets -------------------------------------------------
   level = "full"        # "unit", "fast" or "full"
   case_dir = ""         # complete case folder, relative to test/
   case_name = ""        # name of the .usr/.par/.re2 files
   size_params = {}      # {"lx1": 8}                         -> SIZE (full rebuild)
   par_values = {}       # {"GENERAL": {"endTime": 0.05}}     -> <case>.par
   settings = {}         # {"tau": 0.05}                      -> settings of the .usr
   reference = None      # folder in test/references with values.json
   settings_source = None    # file whose settings `settings` edits [<case_name>.usr]

   # SIZE parameters shown in the summary (all changed ones are shown as well)
   SIZE_SHOWN = ("ldim", "lx1", "lxd", "lx2", "lelg", "lpmin", "ldimt", "lpert",
                 "lorder", "mxprev")

   # Safety net: lines with these markers make the test fail, whatever the
   # checks say (neklab_test.f90 `report` prints NEKLAB_TEST_FAIL)
   FAIL_MARKERS = ("NEKLAB_TEST_FAIL",)

   # check() in neklab_test_helpers.f prints at most this many characters of a name
   CHECK_NAME_MAX = 45

   # --- set up -----------------------------------------------------------------

   def setUp(self):
      env = os.environ
      self.neklab_root = Path(env["NEKLAB_SOURCE_ROOT"])
      self.nek_root = Path(env.get("NEK_SOURCE_ROOT", self.neklab_root / "Nek5000"))
      self.nprocs = int(env.get("NPROCS", "4"))
      self.mpiexec = env.get("MPIEXEC", "mpiexec")
      self.verbose = env.get("VERBOSE_TESTS") == "true"
      self.update_references = env.get("UPDATE_REFERENCES") == "true"

      # Fresh copy of the case: <scratch>/<TestClass>/
      source = (TEST_DIR / self.case_dir).resolve()
      self.source = source
      self.dir = Path(env.get("SCRATCH_ROOT", TEST_DIR / "scratch")) / type(self).__name__
      shutil.rmtree(self.dir, ignore_errors=True)
      shutil.copytree(source, self.dir, ignore=lambda d, names:
                      [n for n in names if any(fnmatch(n, p) for p in IGNORE)])
      print(f"\n=== {type(self).__name__}  ({self.level}: {self.case_dir})", flush=True)
      say("case", f"copied to {self.dir}")

      self.failures = []    # collected by the check* methods
      self.measured = {}    # every value checked; written to measured_values.json
      self.log = None

   @property
   def settings_file(self):
      return self.settings_source or f"{self.case_name}.usr"

   # --- prepare, build, run ----------------------------------------------------

   def prepare(self):
      """Apply size_params, par_values and settings; make the case runnable."""
      size, par = self.dir / "SIZE", self.dir / f"{self.case_name}.par"
      if self.size_params:
         files.set_size_parameters(size, self.size_params)
      # lpmin = minimum number of ranks the SIZE file is dimensioned for
      lpmin = files.read_size(size).get("lpmin", "")
      if lpmin.isdigit() and int(lpmin) > self.nprocs:
         files.set_size_parameters(size, {"lpmin": self.nprocs})
      if self.par_values:
         files.set_par_values(par, self.par_values)
      if self.settings:
         files.set_settings(self.dir / self.settings_file, self.settings)
      if not (self.dir / f"{self.case_name}.ma2").exists():
         say("genmap", f"{self.case_name}.re2 -> {self.case_name}.ma2")
         nek.run_genmap(self.dir, self.case_name, self.nek_root, verbose=self.verbose)
      self.summarize_parameters()

   # --- summary of all parameters ------------------------------------------------

   def summarize_parameters(self):
      """Print SIZE, .par and settings; '*' marks values that differ from the case folder.

      With PARAMS=changed (runner option --params changed) only the changed
      values are listed.
      """
      only_changed = os.environ.get("PARAMS") == "changed"
      say("params", f"* = differs from file in {self.case_dir}/*  (value there in brackets)")

      def block(title, rows, indent=""):
         shown = [r for r in rows if r[0] or not only_changed]
         if not shown and only_changed:
            return
         print(f"   {title}")
         for changed, name, value, original in shown:
            mark = "*" if changed else " "
            extra = f"   ({original})" if changed else ""
            print(f"   {mark}   {indent}{name:24} {value}{extra}")

      # SIZE (changing it recompiles)
      orig = files.read_size(self.source / "SIZE")
      now = files.read_size(self.dir / "SIZE")
      rows = []
      for name, value in now.items():
         changed = not files.same_value(value, orig.get(name))
         if name in self.SIZE_SHOWN or changed:
            rows.append((changed, name, value, orig.get(name, "not set")))
      block("SIZE", rows)

      # .par (read at run time)
      par = f"{self.case_name}.par"
      orig = files.read_par(self.source / par)
      for section, entries in files.read_par(self.dir / par).items():
         before = {k.lower(): v for k, v in orig.get(section, {}).items()}
         rows = []
         for key, value in entries.items():
            old = before.get(key.lower())
            rows.append((old is None or not files.same_value(value, old), key, value,
                         "not set" if old is None else old))
         block(f"{par} [{section}]", rows)

      # settings (initialised variables of the .usr)
      if (self.dir / self.settings_file).exists():
         orig = files.read_settings(self.source / self.settings_file)
         rows = [(not files.same_value(value, orig.get(name)), name, value, orig.get(name))
                 for name, value in files.read_settings(self.dir / self.settings_file).items()]
         block(f"{self.settings_file} settings", rows)

   def build(self):
      """Clean build of the case with makeneklab (no reuse of earlier builds)."""
      say("build", f"makeneklab (clean) in {self.dir}  (log: build.out) ...")
      start = time.time()
      nek.build_neklab(self.dir, self.neklab_root, self.nek_root, verbose=self.verbose)
      say("build", f"done in {time.time() - start:.0f} s")

   def run_nek(self):
      for old in self.dir.glob("eigs_output.txt"):
         old.unlink()
      say("run", f"{self.mpiexec} -np {self.nprocs} ./nek5000  "
                 f"(log: {self.case_name}.log.{self.nprocs}) ...")
      start = time.time()
      self.log = nek.run_nek(self.dir, self.case_name, self.nprocs, self.mpiexec,
                             verbose=self.verbose)
      say("run", f"finished in {time.time() - start:.0f} s")
      return self.log

   # --- checks (all collected, reported together by finish) --------------------

   def check(self, ok, name, info="", quiet=False):
      if not quiet or not ok:
         print(f"    {name:{NAME_WIDTH}}         {status(ok)}  {info}".rstrip(), flush=True)
      if not ok:
         self.failures.append(f"{name}  {info}".rstrip())

   def check_close(self, label, value, target, tol):
      self.measured[label] = value
      self.check(abs(value - target) <= tol, label, f"= {value!r} (target {target!r} +/- {tol})")

   def check_driver(self, expected=(), completed=None):
      """Read the NEKLAB_CHECK / NEKLAB_RESULT lines printed by the case's .usr.

      expected   check names that must be reported
      completed  text the driver prints at its end (a run that stops early fails)

      Any line with one of FAIL_MARKERS also fails the test (safety net for
      failures neklab reports itself, outside the checks of the .usr).
      """
      checks = files.read_checks(self.log)
      results = files.read_results(self.log)
      self.measured.update(results)

      say("checks", f"{len(checks)} reported in {Path(self.log).name}:")
      for name, (status, message) in checks.items():
         self.check(status == "PASS", name, message)
      self.check(bool(checks), "the driver reported its checks", quiet=True)
      for name in expected:
         if name not in checks:
            hint = (f" -- longer than the {self.CHECK_NAME_MAX} characters check() prints"
                    if len(name) > self.CHECK_NAME_MAX else "")
            self.check(False, name, f"(expected, never reported){hint}")
      for line in files.find_lines(self.log, self.FAIL_MARKERS):
         self.check(False, "failure reported by neklab", line)
      if completed:
         reached = files.log_contains(self.log, completed)
         self.check(reached, f'run reached "{completed}"', quiet=True)
         if reached:
            say("checks", f'driver reached "{completed}"')

      if results:
         say("results", f"{len(results)} values (written to measured_values.json):")
         for key, value in results.items():
            shown = f"{value: .5e}" if isinstance(value, float) else value
            print(f"    {key:{NAME_WIDTH}}   {shown}")

   def check_reference(self, key, value, tol):
      """Compare with test/references/<reference>/values.json (or record it in update mode)."""
      self.measured[key] = value
      if self.update_references:
         return
      refs = self._reference_values()
      if key not in refs:
         self.check(False, key, "no reference value yet (run --update-references)")
      else:
         self.check_close(key, value, refs[key], tol)

   def finish(self):
      """Write measured_values.json, save references in update mode, report failures.

      References are written only if every check passed and every file to be
      saved exists (check_reference_outputs), so a failing run never leaves
      half-updated references behind.
      """
      with open(self.dir / "measured_values.json", "w") as f:
         json.dump(self.measured, f, indent=2, sort_keys=True, default=str)
      if self.update_references and self.level == "full":
         if not self.failures:
            self.check_reference_outputs()
         if self.failures:
            say("refs", "not written: the test failed")
         else:
            self._save_reference_values()
      if self.failures:
         say("result", f"{len(self.failures)} check(s) failed")
         self.fail(f"{len(self.failures)} check(s) failed:\n  " + "\n  ".join(self.failures))
      say("result", "all checks passed")

   # --- reference values ---------------------------------------------------------

   def check_reference_outputs(self):
      """Before saving references: add a failure for anything that cannot be saved."""

   def _reference_file(self):
      return TEST_DIR / "references" / (self.reference or self.case_name) / "values.json"

   def _reference_values(self):
      path = self._reference_file()
      return json.loads(path.read_text()) if path.exists() else {}

   def _save_reference_values(self):
      values = self._reference_values()
      values.update(self.measured)
      path = self._reference_file()
      path.parent.mkdir(parents=True, exist_ok=True)
      path.write_text(json.dumps(values, indent=2, sort_keys=True, default=str) + "\n")
      say("refs", f"values -> {path.relative_to(TEST_DIR)}")
