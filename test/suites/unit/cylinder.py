"""Unit tests on the 2D cylinder: test/unit/cylinder/ (1cyl.par, 1cyl.usr, unit_tests.f)."""
from lib.neklabUnit import Case, UnitTest, linop

CYLINDER = Case(
   folder="unit/cylinder",
   name="1cyl",
   par={
      "GENERAL":  {"endTime": 0.1},
      # residual projection makes each solve depend on the previous ones,
      # so the operators would not be exactly linear
      "PRESSURE": {"residualProj": "yes"},
      "VELOCITY": {"residualProj": "yes"},
   },
)


class CylinderFixedPoint(UnitTest):
   """Vectors, exptA, nek_system and its Jacobian, on the steady base flow BF.fld."""
   case = CYLINDER
   group = 1
   checks = [
      "nek_dvector.axioms",
      "nek_zvector.axioms",
      *linop("exptA_linop.axioms"),
      "nek_system.linearization",
      *linop("nek_system%jacobian.axioms"),
   ]


class CylinderResolvent(UnitTest):
   """Resolvent operator, on BF.fld."""
   case = CYLINDER
   group = 2
   checks = [
      "nek_zvector.axioms",
      *linop("resolvent_linop.axioms"),
   ]
   skip = "test_resolvent in unit_tests.f exits at once ('Currently not available')"


class CylinderPeriodicOrbit(UnitTest):
   """Extended vector, UPO system and its Jacobian, on the unsteady snapshot.fld.

   No adjoint check: off a converged orbit jac_adjoint_map is not the transpose
   of jac_direct_map, so it would fail by construction.
   """
   case = CYLINDER
   group = 3
   settings = {"test_adjoint": True}
   checks = [
      "nek_ext_dvector.axioms",
      "nek_upo_system.linearization",
      *linop("nek_upo_system%jacobian.axioms", adjoint=True),
   ]
