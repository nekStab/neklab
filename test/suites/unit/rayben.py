"""Unit tests on Rayleigh-Benard convection: test/unit/rayBen/ (rayBen.par, rayBen.usr, unit_tests.f)."""
from lib.neklabUnit import Case, UnitTest, linop

RAYBEN = Case(
   folder="unit/rayBen",
   name="rayBen",
   par={
      "GENERAL":     {"endTime": 0.1},     # = tau in unit_tests.f
      # residual projection makes the operators history-dependent
      "VELOCITY":    {"residualProj": "no"},
      "PRESSURE":    {"residualProj": "no"},
      "TEMPERATURE": {"residualProj": "no"},
   },
)


class RayleighBenardFixedPoint(UnitTest):
   """Vectors with temperature, exptA_temp, nek_temp_system and its Jacobian,
   on the conduction state u = 0, T = 1 - y (built in unit_tests.f)."""
   case = RAYBEN
   group = 1
   checks = [
      "nek_dvector_temp.axioms",
      *linop("exptA_temp_linop.axioms"),
      "nek_temp_system.linearization",
      *linop("nek_temp_system%jacobian.axioms"),
   ]
