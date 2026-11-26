#!/usr/bin/env python3
"""Run the neklab tests.

   ./neklabTests.py --list                         show the tests
   ./neklabTests.py --levels unit --nprocs 4       run the unit tests
   ./neklabTests.py -k Cylinder                    only tests whose name contains "Cylinder"
   ./neklabTests.py --levels full --update-references   (trusted version only)

Levels: unit (suites/unit/, one file per case), fast and full (both from
suites/regression/: every Analysis gives a fast and a full test).
Each test runs in a copy of its case under test/scratch/; the copies are
deleted when every test passed (keep them with --keep-scratch).
The exit code is 0 only if every selected test passed.

The compiler is chosen by app/makeneklab, not by this script.
"""
import argparse
import importlib
import inspect
import os
import pkgutil
import shutil
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from lib.neklabTestCase import NeklabTestCase, LEVELS

# Where the tests of each level are defined: a package, every file of which is read.
SUITES = {"unit": "suites.unit", "fast": "suites.regression", "full": "suites.regression"}


def level_modules(level):
   """The package of the level and every file in it."""
   module = importlib.import_module(SUITES[level])
   yield module
   for info in pkgutil.iter_modules(getattr(module, "__path__", [])):
      yield importlib.import_module(f"{SUITES[level]}.{info.name}")


def find_tests(levels, patterns):
   """(class, method) of every test of the given levels whose name matches a pattern."""
   tests = []
   for level in levels:
      for module in level_modules(level):
         for _, cls in inspect.getmembers(module, inspect.isclass):
            if issubclass(cls, NeklabTestCase) and cls.__module__ == module.__name__ \
                  and cls.level == level and cls.case_dir \
                  and not cls.__dict__.get("abstract", False):
               for method in unittest.TestLoader().getTestCaseNames(cls):
                  name = f"{cls.__name__}.{method}"
                  if not patterns or any(p.lower() in name.lower() for p in patterns):
                     tests.append((cls, method))
   # an analysis after the ones it needs, full before fast (fast starts from full's fields)
   order = {"unit": 0, "full": 1, "fast": 2}
   return sorted(tests, key=lambda t: (getattr(t[0], "depth", 0), order.get(t[0].level, 9)))


def nothing_selected(levels, patterns):
   where = ", ".join(sorted({SUITES[l].replace(".", "/") + "/" for l in levels}))
   msg = f"No test selected (levels {','.join(levels)}; looked in {where})."
   if patterns:
      msg += f" Check the -k filter {patterns}: ./neklabTests.py --levels unit,fast,full --list"
   return msg


def raise_stack_limit():
   """neklab keeps large work vectors on the stack; the usual 8 MB is not enough
   when few ranks are used. Raise the limit for the nek5000 runs started here."""
   try:
      import resource
      _, hard = resource.getrlimit(resource.RLIMIT_STACK)
      resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))
   except (ImportError, ValueError, OSError):
      pass


def main():
   p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
   p.add_argument("--levels", "--level", default="unit", help="comma-separated: unit,fast,full [unit]")
   p.add_argument("-k", dest="patterns", action="append", default=[], help="name filter (repeatable)")
   p.add_argument("--list", action="store_true", help="list the selected tests and exit")
   p.add_argument("--nprocs", type=int, default=4, help="MPI ranks per run [4]")
   p.add_argument("--mpiexec", default="mpiexec", help="MPI launcher [mpiexec]")
   p.add_argument("--update-references", action="store_true",
                  help="with --levels full only: store the results and fields as the new references")
   p.add_argument("--keep-scratch", action="store_true", help="keep the scratch copies")
   p.add_argument("--params", choices=("all", "changed"), default="all",
                  help="parameters listed before each build: all, or only those that differ from the case folder [all]")
   p.add_argument("-v", "--verbose", action="store_true", help="also echo the full makeneklab and nek5000 output")
   args = p.parse_args()

   levels = [l.strip() for l in args.levels.split(",")]
   if any(l not in LEVELS for l in levels):
      sys.exit(f"Unknown level in {levels}; choose from {LEVELS}.")
   if args.update_references and levels != ["full"]:
      sys.exit("--update-references only works with --levels full "
               "(the fast tests have no references of their own).")

   tests = find_tests(levels, args.patterns)
   if args.list:
      for cls, method in tests:
         skip = getattr(getattr(cls, method), "__unittest_skip_why__", "")
         print(f"{cls.level:5}  {cls.__name__}.{method}" + (f"   [skipped: {skip}]" if skip else ""))
      if not tests:
         print(nothing_selected(levels, args.patterns))
      return 0
   if not tests:
      sys.exit(nothing_selected(levels, args.patterns))
   if "NEKLAB_SOURCE_ROOT" not in os.environ:
      sys.exit("Set NEKLAB_SOURCE_ROOT (and NEK_SOURCE_ROOT if Nek5000 is not in neklab/Nek5000).")

   scratch = HERE / "scratch"
   os.environ.update(
      NPROCS=str(args.nprocs),
      MPIEXEC=args.mpiexec,
      SCRATCH_ROOT=str(scratch),
      VERBOSE_TESTS="true" if args.verbose else "false",
      PARAMS=args.params,
      UPDATE_REFERENCES="true" if args.update_references else "false",
   )
   raise_stack_limit()

   suite = unittest.TestSuite(cls(method) for cls, method in tests)
   # Output is never buffered: every test narrates its steps as they happen.
   # -v additionally echoes the makeneklab and nek5000 output.
   result = unittest.TextTestRunner(verbosity=2, buffer=False).run(suite)

   if result.wasSuccessful() and not args.keep_scratch:
      shutil.rmtree(scratch, ignore_errors=True)
   else:
      print(f"\nScratch copies kept in {scratch}")
   return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
   sys.exit(main())
