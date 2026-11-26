"""Unit tests on the 2D cylinder: one case folder per test in test/unit/.

   cylinder_fixed_point/      on BF.fld        vectors, exptA, nek_system + Jacobian
   cylinder_resolvent/        on BF.fld        resolvent (skipped: not available yet)
   cylinder_periodic_orbit/   on snapshot.fld  extended vector, UPO system + Jacobian

Each folder has its own 1cyl.usr with the checks written in userchk; the
mesh, the fields and neklab_test_helpers.f are symlinks to test/common/.
"""
from lib.neklabUnit import UnitTest, linop


class CylinderUnitTest(UnitTest):
   """What the cylinder tests share (not a test itself: it has no folder)."""
   name = "1cyl"
   par = {
      "GENERAL":  {"endTime": 0.1},
      "PRESSURE": {"residualProj": "yes"},
      "VELOCITY": {"residualProj": "yes"},
   }


class CylinderFixedPoint(CylinderUnitTest):
   """Vectors, exptA, nek_system and its Jacobian, on the steady base flow BF.fld."""
   folder = "cylinder_fixed_point"
   checks = [
      "nek_dvector.axioms",
      "nek_zvector.axioms",
      *linop("exptA_linop.axioms"),
      "nek_system.linearization",
      *linop("nek_system%jacobian.axioms"),
   ]


class CylinderResolvent(CylinderUnitTest):
   """Resolvent operator, on BF.fld."""
   folder = "cylinder_resolvent"
   checks = [
      "nek_zvector.axioms",
      *linop("resolvent_linop.axioms"),
   ]
   skip = "resolvent not available yet: the .usr stops before its checks"


class CylinderPeriodicOrbit(CylinderUnitTest):
   """Extended vector, UPO system and its Jacobian, on the unsteady snapshot.fld.

   The snapshot is not on a periodic orbit, so jac_adjoint_map is only close
   to the transpose of jac_direct_map; with the short guessed period the
   adjoint check passes within adj_tol (see the header of the .usr).
   """
   folder = "cylinder_periodic_orbit"
   checks = [
      "nek_ext_dvector.axioms",
      "nek_upo_system.linearization",
      *linop("nek_upo_system%jacobian.axioms"),
   ]
