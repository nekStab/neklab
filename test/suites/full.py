"""Level 3 -- full non-regression tests.

Each test runs a case from scratch, exactly as a user would, and compares
the result with test/references/<reference>/values.json.

The reference values are produced by running, on the trusted version:

   ./neklabTests.py --levels full --update-references

which records every checked value instead of comparing it.
"""
from lib.neklabFiles import read_eigenvalues
from lib.neklabTestCase import NeklabTestCase


class CylEigsDir(NeklabTestCase):
   """Leading direct eigenvalue of the 2D cylinder wake, from scratch."""

   level = "full"
   case_dir = "../examples/cylinder/stability/direct"   # the example itself
   case_name = "1cyl"
   reference = "cylinder_stability_direct"
   size_params = {}      # e.g. {"lx1": 8} to run at another resolution

   def test_leading_eigenvalue(self):
      self.prepare()
      self.build()
      self.run_nek()

      eigs = read_eigenvalues(self.dir / "eigs_output.txt")
      self.check(len(eigs) > 0, "at least one converged eigenvalue")
      if eigs:
         # |mu_1| of the exponential propagator, and |Im mu_1| (the sign of the
         # imaginary part depends on which of a conjugate pair comes first)
         self.check_reference("lambda_1_modulus", eigs[0]["modulus"], tol=1e-4)
         self.check_reference("lambda_1_abs_imag", abs(eigs[0]["im"]), tol=1e-4)
      self.finish()
