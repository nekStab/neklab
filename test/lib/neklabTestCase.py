"""Base class of all neklab tests.

A test is a complete case folder (SIZE, .par, .usr, mesh, input fields) plus
a few values to change in it. The test:

   1. copies the folder to a scratch directory (the original is never touched),
   2. changes the values listed in size_params / par_values / usr_params,
   3. builds and runs the case,
   4. checks the results.

Example:

   class UnitCylinder(NeklabTestCase):
      level = "unit"
      case_dir = "unit/cylinder"            # relative to test/
      case_name = "1cyl"
      par_values = {"GENERAL": {"endTime": 0.05}}
      usr_params = {"ntrials": 2}

      def test_axioms(self):
         self.prepare()
         self.build()
         self.run_nek()
         self.check_driver(expected=["vector.nek_dvector.axioms"])
         self.finish()
"""
import hashlib
import json
import os
import shutil
import unittest
from fnmatch import fnmatch
from pathlib import Path

from lib import neklabBuildRun as nek
from lib import neklabFiles as files

TEST_DIR = Path(__file__).resolve().parent.parent      # the test/ folder
LEVELS = ("unit", "fast", "full")

# Never copied into the scratch case: build products and outputs of earlier runs
IGNORE = ("obj", "nek5000", "makefile", "makefile_usr.inc", "*.log.*", "*.out",
          "SESSION.NAME", "ioinfo", "eigs_output.txt", "nek_dummy.usr",
          ".neklab_stash", "*.mod", "*.smod")


class NeklabTestCase(unittest.TestCase):

   # nek5000 executables built in this session, by build key (see build)
   _builds = {}

   # --- what a test class sets -------------------------------------------------
   level = "full"        # "unit", "fast" or "full"
   case_dir = ""         # complete case folder, relative to test/
   case_name = ""        # name of the .usr/.par/.re2 files
   size_params = {}      # {"lx1": 8}                         -> SIZE
   par_values = {}       # {"GENERAL": {"endTime": 0.05}}     -> <case>.par
   usr_params = {}       # {"tau": 0.05}  (Fortran parameters) -> <case>.usr
   fortran_params = {}   # {"unit_tests.f": {"tau": 0.05}}   -> any other source file
   reference = None      # folder in test/references with values.json

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
      self.dir = Path(env.get("SCRATCH_ROOT", TEST_DIR / "scratch")) / type(self).__name__
      shutil.rmtree(self.dir, ignore_errors=True)
      shutil.copytree(source, self.dir, ignore=lambda d, names:
                      [n for n in names if any(fnmatch(n, p) for p in IGNORE)])
      print(f"Case copied to {self.dir}")

      self.failures = []    # collected by the check* methods
      self.measured = {}    # every value checked; written to measured_values.json
      self.log = None

   # --- prepare, build, run ----------------------------------------------------

   def prepare(self):
      """Apply size_params, par_values and usr_params; make the case runnable."""
      size = self.dir / "SIZE"
      if self.size_params:
         files.set_fortran_parameters(size, self.size_params)
      # lpmin = minimum number of ranks the SIZE file is dimensioned for
      lpmin = files.get_fortran_parameter(size, "lpmin")
      if lpmin and lpmin.isdigit() and int(lpmin) > self.nprocs:
         files.set_fortran_parameters(size, {"lpmin": self.nprocs})
         print(f"SIZE: lpmin {lpmin} -> {self.nprocs} (ranks used by this test)")
      if self.par_values:
         files.set_par_values(self.dir / f"{self.case_name}.par", self.par_values)
      if self.usr_params:
         files.set_fortran_parameters(self.dir / f"{self.case_name}.usr", self.usr_params)
      for name, params in self.fortran_params.items():
         if params:
            files.set_fortran_parameters(self.dir / name, params)
      if not (self.dir / f"{self.case_name}.ma2").exists():
         nek.run_genmap(self.dir, self.case_name, self.nek_root, verbose=self.verbose)

   def build(self):
      """Build the case, or reuse the executable of an identical earlier build.

      Two tests share a build when all their Fortran sources (SIZE, .usr and
      the files it includes) are identical after prepare() -- e.g. several
      groups of one case selected by a .par value.
      """
      key = self._build_key()
      earlier = NeklabTestCase._builds.get(key)
      if earlier is not None and (earlier / "nek5000").exists():
         shutil.copy2(earlier / "nek5000", self.dir / "nek5000")
         print(f"Reusing nek5000 built in {earlier}")
         return
      nek.build_neklab(self.dir, self.neklab_root, self.nek_root, verbose=self.verbose)
      NeklabTestCase._builds[key] = self.dir

   SOURCES = ("SIZE", "*.usr", "*.f", "*.f90", "*.F", "*.F90", "*.inc")

   def _build_key(self):
      h = hashlib.sha1(str(self.neklab_root).encode())
      for path in sorted(self.dir.iterdir()):
         if path.is_file() and any(fnmatch(path.name, p) for p in self.SOURCES):
            h.update(path.name.encode())
            h.update(path.read_bytes())
      return h.hexdigest()

   def run_nek(self):
      for old in self.dir.glob("eigs_output.txt"):
         old.unlink()
      self.log = nek.run_nek(self.dir, self.case_name, self.nprocs, self.mpiexec,
                             verbose=self.verbose)
      return self.log

   # --- checks (all collected, reported together by finish) --------------------

   def check(self, ok, message):
      print(("  ok    " if ok else "  FAIL  ") + message)
      if not ok:
         self.failures.append(message)

   def check_close(self, label, value, target, tol):
      self.measured[label] = value
      self.check(abs(value - target) <= tol, f"{label} = {value!r} (target {target!r} +/- {tol})")

   def check_driver(self, expected=(), completed=None):
      """Read the NEKLAB_CHECK / NEKLAB_RESULT lines printed by the case's .usr.

      expected   check names that must be reported
      completed  text the driver prints at its end (a run that stops early fails)
      """
      checks = files.read_checks(self.log)
      self.measured.update(files.read_results(self.log))
      self.check(bool(checks), "the driver reported its checks")
      if completed:
         self.check(files.log_contains(self.log, completed), f'the driver reached "{completed}"')
      for name, (status, message) in sorted(checks.items()):
         self.check(status == "PASS", f"{name} {message}".strip())
      for name in expected:
         self.check(name in checks, f"{name} was reported")

   def check_reference(self, key, value, tol):
      """Compare with test/references/<reference>/values.json (or record it in update mode)."""
      self.measured[key] = value
      if self.update_references:
         return
      refs = self._reference_values()
      if key not in refs:
         self.check(False, f"{key}: no reference value yet (run --update-references)")
      else:
         self.check_close(key, value, refs[key], tol)

   def finish(self):
      """Write measured_values.json, save references in update mode, report failures."""
      with open(self.dir / "measured_values.json", "w") as f:
         json.dump(self.measured, f, indent=2, sort_keys=True, default=str)
      if self.update_references and self.level == "full" and not self.failures:
         self._save_reference_values()
      if self.failures:
         self.fail(f"{len(self.failures)} check(s) failed:\n  " + "\n  ".join(self.failures))

   # --- reference values ---------------------------------------------------------

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
      print(f"Reference values written to {path}")
