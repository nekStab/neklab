"""What a unit test definition is built from (used by suites/unit/*.py).

A unit test runs one group of checks of a case folder and passes if every
check the driver reports passes and the driver gets to its end. See
suites/unit/__init__.py for how to add one.
"""
import unittest
from dataclasses import dataclass, field

from lib.neklabTestCase import NeklabTestCase


@dataclass(frozen=True)
class Case:
   """A case folder in test/unit/ and what all of its tests share."""

   folder: str                                   # relative to test/, e.g. "unit/cylinder"
   name: str                                     # basename of the .usr/.par/.re2 files
   par: dict = field(default_factory=dict)       # .par changes for every test of the case
   settings: dict = field(default_factory=dict)  # unit_tests.f settings for every test
   group_key: str = "userParam01"                # .par key read by run_tests(group)


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
   """Base of the definitions: a test is a Case, a group and the expected checks."""

   abstract = True
   level = "unit"

   # --- what a definition sets ----------------------------------------------------
   case = None          # Case
   group = None         # number passed to run_tests() in the case's unit_tests.f
   checks = ()          # check names the driver must report
   par = {}             # .par changes for this test (on top of case.par)
   settings = {}        # unit_tests.f settings for this test (on top of case.settings)
   size = {}            # SIZE changes (full rebuild)
   skip = None          # reason to skip this test, or None

   completed = "tests completed"     # printed by tests_completed() at the end of a run

   def __init_subclass__(cls, **kwargs):
      super().__init_subclass__(**kwargs)
      if cls.case is not None:
         cls.case_dir, cls.case_name = cls.case.folder, cls.case.name
      if cls.skip:
         cls.test_checks = unittest.skip(cls.skip)(cls.test_checks)

   def prepare(self):
      group = {"GENERAL": {self.case.group_key: self.group}} if self.group is not None else {}
      self.par_values = merged(self.case.par, group, self.par)
      self.settings = merged(self.case.settings, type(self).settings)
      self.size_params = self.size
      super().prepare()

   def test_checks(self):
      self.prepare()
      self.build()
      self.run_nek()
      self.check_driver(expected=self.checks, completed=self.completed)
      self.finish()
