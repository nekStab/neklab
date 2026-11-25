'''Level 1 -- unit tests of the neklab types: one case folder per test.

   suites/unit/cylinder.py   test/unit/cylinder_*/   2D cylinder (isothermal)
   suites/unit/rayben.py     test/unit/rayBen_*/     Rayleigh-Benard (temperature)

Each test is a complete, independent case folder whose .usr runs LightKrylov's
verify_vector_axioms, verify_linop_axioms and verify_system_linearization on
neklab's types, directly in userchk. A test passes if every NEKLAB_CHECK the
run prints passes, every check in `checks` is reported, and the run reaches
tests_completed().

   ./neklabTests.py --list                    show the tests
   ./neklabTests.py -k CylinderFixedPoint     run one test

A folder also runs by hand, without the harness (makeneklab, then mpiexec).


FOLDERS
-------
   test/unit/common/                  shared files: neklab_test_helpers.f, meshes,
                                      SIZE, converged base flows, snapshots
   test/unit/cylinder_fixed_point/    1cyl.usr, 1cyl.par + symlinks to common/
   ...

Every test folder holds its own <name>.usr (and .par); everything else is a
symlink into common/. The harness copies the folder to test/scratch/<Test>/
and follows the symlinks, so the scratch copy is self-contained.


ADDING A TEST
-------------
1. Fortran: a new folder test/unit/<case>_<what>/ with
     - <name>.usr: start from an existing one. It includes
       'neklab_test_helpers.f'; userchk declares its settings as initialised
       variables (`real(dp) :: tau = 0.1_dp`, never `parameter`, names unique
       in the file), logs them with log_setting, runs the checks with
       check / check_linop / check_taylor, and ends with tests_completed();
     - <name>.par;
     - symlinks to the files it needs in ../common/.

2. Python, in suites/unit/<case>.py:

     class CylinderMyType(CylinderUnitTest):
        """What it checks, on which state."""
        folder = "cylinder_my_type"
        checks = ["mytype.axioms", *linop("mytype_linop.axioms")]

   Optional: par = {...}, settings = {...}, size = {...}, skip = "reason".

The test is found automatically (./neklabTests.py --list).


CHANGING PARAMETERS
-------------------
   par      -> <name>.par                       {"GENERAL": {"endTime": 0.1}}
   settings -> settings at the top of userchk   {"ntrials_lin": 2, "tau": 0.05}
   size     -> SIZE (full rebuild)              {"lx1": 8}

The shared class of a case (CylinderUnitTest, ...) holds `name` and `par`;
a test can override them. Every test prints all values before it builds
('*' = differs from the case folder), and its log repeats the settings the
run used.
'''
