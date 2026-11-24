'''Levels 2 and 3 -- fast and full regression tests: one file per case.

   suites/regression/cylinder.py   test/regression/cylinder_steady/  (and _periodic/)
   suites/regression/examples.py   the examples/ themselves, unchanged

Each Analysis (lib/neklabRegression.py) is one solver run of a case folder and
gives two tests that run the same .usr with different settings:

   full   from scratch; compares with test/references/<reference>/values.json
   fast   from a reference field + noise; same answer (looser), few iterations

   ./neklabTests.py --levels fast            quick check during the port
   ./neklabTests.py --levels full            before merging
   ./neklabTests.py --levels full,fast --update-references     trusted version only


ADDING AN ANALYSIS
------------------
1. Fortran, in test/regression/<case>/<case>.usr (a copy of the example's .usr):
     - the settings fast and full differ in are initialised variables
       (`logical :: if_load_BF = .true.`), never `parameter`, with names
       unique in the file;
     - userchk prints NEKLAB_RESULT <key> <value> for every compared value,
       NEKLAB_CHECK lines with check(), and ends with tests_completed().

2. Python, in suites/regression/<case>.py:

     MY_ANALYSIS = Analysis(
        case=CYLINDER_STEADY, name="MyAnalysis", reference="cylinder/my_analysis",
        inputs={"BF.fld": "cylinder/newton_fp/BF_ref.fld"}, needs=(NEWTON_FP,),
        saves_fields={"<file the run writes>": "mode1_ref.fld"},
        compare={"lambda_1_re": 1e-6}, bounds={"mode_overlap": (0.999, 1.001)},
        fast=dict(settings={...}, max={"iterations": 3}),
        full=dict(settings={...}),
     )

   and add it to make_tests(...) at the end of the file.


REFERENCES
----------
On the trusted version, run --update-references for full first (it writes the
compared values and the converged fields that the fast tests start from),
then for fast (it records the `max` values as fast.<key>). The runner orders
the tests itself: full before fast, and an analysis after its `needs`.
'''