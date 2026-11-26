"""Unit tests on Rayleigh-Benard convection: one case folder per test in test/unit/.

   rayBen_fixed_point/   on u = 0, T = 1 - y   vectors with temperature,
                                               exptA_temp, nek_temp_system + Jacobian

The mesh and neklab_test_helpers.f are symlinks to test/common/.
"""
from lib.neklabUnit import UnitTest, linop


class RayleighBenardUnitTest(UnitTest):
   """What the Rayleigh-Benard tests share (not a test itself: it has no folder)."""
   name = "rayBen"
   par = {
      "GENERAL":     {"endTime": 0.1},     # = tau in rayBen.usr
      "VELOCITY":    {"residualProj": "no"},
      "PRESSURE":    {"residualProj": "no"},
      "TEMPERATURE": {"residualProj": "no"},
   }


class RayleighBenardFixedPoint(RayleighBenardUnitTest):
   """Vectors with temperature, exptA_temp, nek_temp_system and its Jacobian,
   on the conduction state u = 0, T = 1 - y (built in rayBen.usr)."""
   folder = "rayBen_fixed_point"
   checks = [
      "nek_dvector_temp.axioms",
      *linop("exptA_temp_linop.axioms"),
      "nek_temp_system.linearization",
      *linop("nek_temp_system%jacobian.axioms"),
   ]
