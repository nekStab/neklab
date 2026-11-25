"""What a unit test definition is built from (used by suites/unit/*.py).

A unit test is one complete case folder test/unit/<folder>/ with its own
.usr: the checks are written directly in userchk, and the settings are the
initialised variables at its top. The test passes if every check the run
reports passes, every name in `checks` is reported and the run reaches
tests_completed(). See suites/unit/__init__.py for how to add one.
"""
import unittest
from dataclasses import dataclass, field

from lib.neklabTestCase import NeklabTestCase


@dataclass(frozen=True)
class Case:
   """A case folder and what all of its analyses share.

   Only used by lib/neklabRegression.py; the unit tests have one folder each.
   """

   folder: str                                   # relative to test/, e.g. "regression/cylinder_steady"
   name: str                                     # basename of the .usr/.par/.re2 files
   par: dict = field(default_factory=dict)       # .par changes for every analysis of the case
   settings: dict = field(default_factory=dict)  # settings for every analysis
   group_key: str = "userParam01"                # .par key for Analysis.group


def linop(name, adjoint=True):
   """The check names that check_linop prints for one linear operator."""
   return ([f"{name}.linearity"] + ([f"{name}.adjoint"] if adjoint else []))


def merged(*dicts):
   """Merge nested dicts left to right: merged({"A": {"x": 1}}, {"A": {"y": 2}})."""
   out = {}
   for d in dicts:
      for key, value in d.items():
         out[key] = merged(out.get(key, {}), value) if isinstance(value, dict) else value
   return out


class UnitTest(NeklabTestCase):
   """Base of the definitions: a test is a folder, its .usr and the expected checks.

   A class without `folder` is not a test; use one to share `name` and `par`
   between the tests of a case (see suites/unit/cylinder.py).
   """

   level = "unit"

   # --- what a definition sets ----------------------------------------------------
   folder = None        # case folder in test/unit/, e.g. "cylinder_fixed_point"
   name = None          # basename of its .usr/.par/.re2 files, e.g. "1cyl"
   checks = ()          # check names the run must report
   par = {}             # .par changes
   settings = {}        # changes to the settings at the top of userchk in <name>.usr
   size = {}            # SIZE changes
   skip = None          # reason to skip this test, or None

   completed = "tests completed"     # printed by tests_completed() at the end of a run

   def __init_subclass__(cls, **kwargs):
      super().__init_subclass__(**kwargs)
      if cls.folder:
         cls.case_dir = f"unit/{cls.folder}"
         cls.case_name = cls.name
         cls.settings_source = f"{cls.name}.usr"
      if cls.skip:
         cls.test_checks = unittest.skip(cls.skip)(cls.test_checks)

   def prepare(self):
      self.par_values = self.par
      self.size_params = self.size
      super().prepare()

   def test_checks(self):
      self.prepare()
      self.build()
      self.run_nek()
      self.check_driver(expected=self.checks, completed=self.completed)
      self.finish()
