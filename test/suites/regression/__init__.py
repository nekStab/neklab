'''Levels 2 and 3 -- fast and full regression tests: one file per case.

   suites/regression/cylinder.py   test/regression/cylinder_steady_*/

Each Analysis (lib/neklabRegression.py) is one solver run of a case folder and
gives two tests that run the same .usr with different settings:

   full   from scratch; compares with test/references/<reference>/values.json
   fast   from a reference field + noise; same answer (looser tolerance)

   ./neklabTests.py --levels fast            quick check during the port
   ./neklabTests.py --levels full            before merging
   ./neklabTests.py --levels full --update-references     trusted version only


ADDING AN ANALYSIS
------------------
1. Fortran: a case folder test/regression/<case>_<what>/ whose <name>.usr
     - includes 'neklab_test_helpers.f' (first line; symlink to test/common/);
     - declares the settings fast and full differ in as initialised variables
       (`logical :: if_load_BF = .true.`), never `parameter`, with names
       unique in the file;
     - prints NEKLAB_RESULT lines with report() for every compared value and
       NEKLAB_CHECK lines with check(), and ends with tests_completed();
     - reads its reference fields (ref_*.fld) only if they exist: they are
       missing on the first --update-references.

2. Python, in suites/regression/<case>.py:

     MY_ANALYSIS = Analysis(
        case=Case(folder="regression/<case>_<what>", name="1cyl"),
        name="MyAnalysis", reference="cylinder/my_analysis",
        needs=(NEWTON_FP,),                                    # runs after it
        inputs={"BF.fld": "cylinder/newton_fp/BF_ref.fld"},    # required
        reference_fields={"ref_mode1.fld": "cylinder/my_analysis/mode1.fld"},
        saves_fields={"<file the run writes>": "mode1.fld"},
        compare={"lambda_1_re": 1e-6},                         # vs values.json
        field_errors={"mode_1_err": 1e-4},                     # vs reference_fields
        fast=dict(settings={...}, max={"iterations": 3}),
        full=dict(settings={...}),
     )

   and add it to make_tests(...) at the end of the file.


REFERENCES
----------
On the trusted version:

   ./neklabTests.py --levels full --update-references

Every full test that passes writes its `compare` values to values.json and
copies its saves_fields; a test that fails writes nothing. The runner orders
the tests itself (an analysis after its `needs`), and in this mode an input
that is missing from the references is an error, never replaced by the file
of the case folder. Reference fields are not compared in this mode.
The fast tests have no references of their own: their `max` bounds are set
in the suite file.
'''
