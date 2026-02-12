#!/usr/bin/env python3

import sys
from shutil import copyfile
import os

if sys.version_info < (3, 6):
   print("Sorry, requires Python > 3.6")
   sys.exit(1)

from lib.neklabTestCase import (
   NeklabTestCase,
   pn_pn_2_parallel,
)

class CylEigsDir(NeklabTestCase):
    
   example_subdir = os.path.join('cylinder','stability', 'direct')
   case_name = "1cyl"

   def setUp(self):
      self.size_params = dict(
         ldim="2",
         lpmin="12",
         lx1="6",
         lxd="9",
         lx2="lx1-2",
         lelg="1996",
         lpelt="lelt",
      )

      self.build_tools(["genmap"])
      self.run_genmap()

   @pn_pn_2_parallel
   def test_PnPn2_Parallel(self):
      self.size_params["lx2"] = "lx1-2"
      self.config_size()
      self.build_neklab(usr_file="1cyl")
      self.run_nek(step_limit=None)

      eigs = self.get_converged_eigs_data()
      self.assertAlmostEqualDelayed(
         eigs['lambda_1']['modulus'], target_val=1.0156, delta=1e-04, label="growth rate"
      )

      self.assertDelayedFailures()

if __name__ == "__main__":
   import unittest
   import argparse

   # Get arguments from command line
   parser = argparse.ArgumentParser()
   parser.add_argument(
      "--f77",
      default="mpif77",
      help="The Fortran 77 compiler to use [default: mpif77]",
   )
   parser.add_argument(
      "--cc", default="mpicc", help="The C compiler to use [default: mpicc]"
   )
   parser.add_argument(
      "--ifmpi",
      default="true",
      choices=["true", "false"],
      help="Enable/disable parallel tests with MPI [default: true]",
   )
   parser.add_argument(
      "--nprocs",
      default="4",
      help="Number of processes to use for MPI tests [default: 4]",
   )
   parser.add_argument(
      "-v", "--verbose", action="store_true", help="Enable verbose output"
   )

   args = parser.parse_args()

   # # Set environment
   os.environ["CC"] = args.cc
   os.environ["FC"] = args.f77
   os.environ["IFMPI"] = args.ifmpi
   os.environ["PARALLEL_PROCS"] = 12 #args.nprocs
   if args.verbose:
      os.environ["VERBOSE_TESTS"] = "true"
      ut_verbose = 2
   else:
      os.environ["VERBOSE_TESTS"] = "false"
      ut_verbose = 1

   testList = (
      CylEigsDir,
   )

   suite = unittest.TestSuite(
      [unittest.TestLoader().loadTestsFromTestCase(t) for t in testList]
   )
   unittest.TextTestRunner(verbosity=ut_verbose, buffer=True).run(suite)
