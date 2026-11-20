"""Level 1 -- unit tests of the neklab types.

Each test is a complete case folder in test/unit/<case>/ (SIZE, .par, .usr,
mesh, input fields). Its .usr runs LightKrylov's checks in userchk:

   verify_vector_axioms          vectors
   verify_linop_axioms           linear operators, including Jacobians
   verify_system_linearization   systems (Taylor test)

and prints one line per check:  NEKLAB_CHECK PASS|FAIL <name> [message].
The test fails if a check fails, if the run crashes, or if a check listed in
`expected` is never reported.

Changing parameters from here (the case files hold the defaults):

   size_params     -> SIZE                       {"lx1": 8}
   par_values      -> <case>.par                 {"GENERAL": {"endTime": 0.1}}
   usr_params      -> Fortran parameters in <case>.usr
   fortran_params  -> declarations in other files of the case
                      {"unit_tests.f": {"ntrials_lin": 2}}

Any name declared with a value can be changed (`parameter`s and initialised
variables such as the settings of unit_tests.f); a misspelled name stops the
test with a KeyError. Booleans become .true./.false. Floats are written with the _dp kind, strings
exactly as given (so a file name needs its quotes: "'BF.fld'").

.par values are read at run time; everything else is compiled in. Tests whose
sources end up identical share one build, so prefer .par values for what
differs between groups, and change compiled settings for all groups at once.
"""
import unittest

from lib.neklabTestCase import NeklabTestCase


class UnitTest(NeklabTestCase):
   """What every unit test does: build the case, run it once, read its checks."""

   level = "unit"
   expected = ()
   completed = None      # text printed at the end of the driver, if any

   def test_checks(self):
      self.prepare()
      self.build()
      self.run_nek()
      self.check_driver(expected=self.expected, completed=self.completed)
      self.finish()


def merged(*dicts):
   """Merge nested dicts left to right: merged({"A": {"x": 1}}, {"A": {"y": 2}})."""
   out = {}
   for d in dicts:
      for key, value in d.items():
         if isinstance(value, dict):
            out[key] = merged(out.get(key, {}), value)
         else:
            out[key] = value
   return out


# --- Cylinder: unit/cylinder (1cyl.par, 1cyl.usr + unit_tests.f) -------------
#
# One case, three groups selected with userParam01 (run time, one build).
# The settings below apply to every group; each group can add its own.

CYLINDER_PAR = {
   "GENERAL":  {"endTime": 0.1},
   # Residual projection makes each solve depend on the previous ones, so the
   # operators are no longer exactly linear (fails the linearity axioms).
   "PRESSURE": {"residualProj": "no"},
   "VELOCITY": {"residualProj": "no"},
}

CYLINDER_SETTINGS = {
   # Settings block of unit_tests.f (defaults are in the file). Examples:
   # "ntrials_lin": 2,
   # "ptol":        1e-8,
   # "lin_tol":     1e-6,
   # "adj_tol":     1e-2,
   # "eps0":        1.0,
   # "verb":        False,
   # "base_flow":   "'BF.fld'",
}


class CylinderGroup(UnitTest):
   """Common part of the cylinder groups (not a test itself)."""

   abstract = True
   case_dir = "unit/cylinder"
   case_name = "1cyl"
   completed = "tests completed"       # printed by tests_completed()
   group = None      # userParam01 in 1cyl.par
   par = {}          # .par changes for this group, on top of CYLINDER_PAR
   settings = {}     # unit_tests.f changes for this group (forces its own build)

   def prepare(self):
      self.par_values = merged(CYLINDER_PAR, {"GENERAL": {"userParam01": self.group}},
                               self.par)
      self.fortran_params = {"unit_tests.f": merged(CYLINDER_SETTINGS, self.settings)}
      super().prepare()


def linop_checks(name, adjoint=True):
   """Names printed by check_linop in unit_tests.f."""
   return ((f"{name}.linearity",) + ((f"{name}.adjoint",) if adjoint else ())
           + (f"GLOBAL : {name}",))


class UnitCylinderFixedPoint(CylinderGroup):
   """Group 1: fixed point -- vectors, exptA, nek_system and its Jacobian (BF.fld)."""

   group = 1
   expected = (
      "nek_dvector.axioms",
      "nek_zvector.axioms",
      *linop_checks("exptA_linop.axioms"),
      "nek_system.linearization",
      *linop_checks("nek_system%jacobian.axioms"),
   )


class UnitCylinderResolvent(CylinderGroup):
   """Group 2: resolvent operator (BF.fld)."""

   group = 2
   expected = (
      "nek_zvector.axioms",
      *linop_checks("resolvent_linop.axioms"),
   )

   @unittest.skip("test_resolvent in unit_tests.f exits at once ('Currently not available')")
   def test_checks(self):
      super().test_checks()


class UnitCylinderPeriodicOrbit(CylinderGroup):
   """Group 3: periodic orbit -- extended vector, UPO system and its Jacobian (snapshot.fld).

   The adjoint is off: off the orbit jac_adjoint_map is not the transpose of
   jac_direct_map (base flow integrated forward, border not transposed), so
   that check would fail by construction. Remove the line to test it anyway.
   Changing a compiled setting gives this group its own build.
   """

   group = 3
   settings = {"test_adjoint": False}
   expected = (
      "nek_ext_dvector.axioms",
      "nek_upo_system.linearization",
      *linop_checks("nek_upo_system%jacobian.axioms", adjoint=False),
   )


class UnitRayleighBenard(UnitTest):
   """Temperature: vectors with theta, exptA_temp, temperature fixed-point system."""

   case_dir = "unit/rayBen"
   case_name = "rayBen"
   par_values = {"GENERAL": {"endTime": 0.05}}
   usr_params = {"tau": 0.05, "ntrials_lin": 2}
   expected = (
      "setup.temperature_active",
      "vector.nek_dvector_temp.axioms",
      "setup.base_flow",
      "linop.exptA_temp.linearity", "linop.exptA_temp.adjoint",
      "system.nek_temp_system.taylor",
      "linop.nek_temp_jacobian.linearity", "linop.nek_temp_jacobian.adjoint",
   )
