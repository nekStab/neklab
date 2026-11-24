'''Level 1 -- unit tests of the neklab types: one file per case folder.

   suites/unit/cylinder.py   test/unit/cylinder/   2D cylinder (isothermal)
   suites/unit/rayben.py     test/unit/rayBen/     Rayleigh-Benard (temperature)

Each test runs one group of checks of a case: LightKrylov's verify_vector_axioms,
verify_linop_axioms and verify_system_linearization on neklab's types. A test
passes if every NEKLAB_CHECK the run prints passes, every check in `checks` is
reported, and the run reaches tests_completed().


ADDING A TEST FOR YOUR TYPE
---------------------------
1. Fortran, in test/unit/<case>/unit_tests.f:
     - write   subroutine test_mytype()   using check / check_linop / check_taylor
       (see the existing drivers), and make it public;
     - add one line to run_tests:   case (4); call test_mytype()

2. Python, in suites/unit/<case>.py:

     class CylinderMyType(UnitTest):
        """What it checks, on which state."""
        case = CYLINDER
        group = 4
        checks = ["mytype.axioms", *linop("mytype_linop.axioms")]

   Optional: par = {...}, settings = {...}, size = {...}, skip = "reason".

The test is found automatically (./neklabTests.py --list).


ADDING A CASE
-------------
A complete case folder test/unit/<case>/ (SIZE, .par, .usr, mesh, fields),
whose .usr includes 'unit_tests.f' and calls in userchk

     call log_settings()
     call run_tests(nint(uparam(1)))
     call tests_completed()

plus test/unit/<case>/neklab_test_helpers.f -> ../common/neklab_test_helpers.f
(symlink) and a new suites/unit/<case>.py with a Case(...) and its tests.


CHANGING PARAMETERS
-------------------
   par      -> <case>.par                       {"GENERAL": {"endTime": 0.1}}
   settings -> Settings block of unit_tests.f   {"ntrials_lin": 2, "base_flow": "BF.fld"}
   size     -> SIZE (full rebuild)              {"lx1": 8}

Case(...) holds what all tests of a case share; a test adds its own on top.
Every test prints all values before it builds ('*' = differs from the case
folder), and its log repeats the settings and Nek5000 setup the run used.
'''
