import unittest
import inspect
import os, sys
from functools import wraps

def pn_pn_2_parallel(method):
   @wraps(method)
   def wrapper(self, *args, **kwargs):
      self.mpi_procs = self.parallel_procs
      if not self.ifmpi:
         self.skipTest(f'Skipping "{self.id()}"; MPI is not enabled.')
      else:
         # Set number of mpi procs
         self.log_suffix = ".pn_pn_2"
         if self.ifmpi:
               self.log_suffix += ".parallel"
         else:
               self.log_suffix += ".serial"
         method(self, *args, **kwargs)

   return wrapper

###############################################################################
#  BASE TEST CASE
###############################################################################


class NeklabTestCase(unittest.TestCase):
   """Base class for Neklab unittests

   This defines a setUpClass method to:
      (a) get the relevant environment variables for compilers, directories
      (b) add env vars to maketools, makenek
      (b) build tools
   All subclassed TestCases will need to do these things.

   Class attributes:
      f77 (str):                 The Fortran 77 compiler to use     [default: 'gfortran']
      cc (str):                  The C compiler to use              [default: 'gcc']
      ifmpi (bool):              Perform compilation/tests with MPI [default: False]
      nek_source_root (str):     Path to Nek source directory;overridden by $NEK_SOURCE_ROOT env variable
      neklab_source_root (str):  Path to Neklab source directory;overridden by $NEKLAB_SOURCE_ROOT env variable
      tools_root (str):          Path to Nek tools directory; overridden by $TOOLS_ROOT env variable
      examples_root (str):       Path to Nek examples directory; overridden by $EXAMPLES_ROOT env variable
      makeneklab (str):          Path to makeneklab                 [default: neklab_source_root/app/makeneklab]
      makenek (str):             Path to makenek                    [default: nek_source_root/makenek]
      tools_bin (str):           Directory to place compiled tools  [default: tools_root/bin]

   Subclass attributes:
      These aren't meaningful in the base class.  They're intended for a subclass that represents
      a particular example problem.
      example_subdir (str):      The subdirectory for the subclass' example.  Assumed that it's in example_root
      rea_file (str):            The .rea file for the subclass' example, minus the '.rea' extension.  Assumed
                                 that it's in example_root/example_dir
      size_file (str):           The SIZE file for the subclass' example.  Assumed that it's in
                                 example_root/example_subdir
   """

   # Defined in subclasses only; declared here to make syntax checker happy
   example_subdir = ""
   case_name = ""

   def __init__(self, *args, **kwargs):
      # These can be overridden by self.get_opts
      self.f77 = ""
      self.cc = ""
      self.pplist = ""
      self.usr_lflags = ""
      self.ifmpi = True

      self.neklab_source_root = os.path.dirname(os.path.dirname(inspect.getabsfile(self.__class__)))
      self.examples_root = os.path.join(self.neklab_source_root, "examples")
      self.makeneklab = os.path.join(self.neklab_source_root, "app", "makeneklab")
      self.nek_source_root = os.path.join(self.neklab_source_root, "Nek5000")
      self.tools_root = os.path.join(self.nek_source_root, "tools")
      self.tools_bin = os.path.join(self.nek_source_root, "bin")
      self.log_root = ""
      self.verbose = True
      self.serial_procs = 1
      self.parallel_procs = 12
      self.size_params = {}

      # These are overridden by method decorators (pn_pn_serial, pn_pn_parallel,
      # pn_pn_2_serial, and pn_pn_2_parallel)
      self.log_suffix = ""
      self.mpi_procs = None

      # Empy list of delayed fails
      self._delayed_failures = []

      self.get_opts()

      unittest.TestCase.__init__(self, *args, **kwargs)

   def assertAlmostEqualDelayed(self, test_val, target_val, delta, label):
      if abs(test_val - target_val) <= delta:
         msg = "    SUCCESS: {}: Test value {} equals target value {} +/- {}".format(
               label, test_val, target_val, delta
         )
      else:
         msg = (
               "    FAILURE: {}: Test value {} exceeds target value {} +/- {}".format(
                  label, test_val, target_val, delta
               )
         )
         self._delayed_failures.append(msg)
      print(msg)

   def assertIsNotNullDelayed(self, test_val, label):
      if test_val:
         msg = f'SUCCESS: Found phrase "{label}" in logfile.'
      else:
         msg = 'FAILURE: Unexpectedly did not find phrase "{}" in logfile'.format(
               label
         )
         self._delayed_failures.append(msg)
      print(msg)

   def assertIsNullDelayed(self, test_val, label):
      if test_val:
         msg = f'FAILURE: Found phrase "{label}" in logfile.'
         self._delayed_failures.append(msg)
      else:
         msg = f'SUCCESS: Did not find phrase "{label}" in logfile'
      print(msg)

   def assertDelayedFailures(self):
      if self._delayed_failures:
         report = ["\n\nFailed assertions:{}\n".format(len(self._delayed_failures))]
         for i, failure in enumerate(self._delayed_failures, start=1):
               report.append(f"{i}: {failure}")
         # self._delayed_failures = []
         self.fail("\n".join(report))

   def get_opts(self):

      print("Getting setup options...")

      # Get compiler options from env
      self.f77 = os.environ.get("FC", self.f77)
      self.cc = os.environ.get("CC", self.cc)
      self.pplist = os.environ.get("PPLIST", self.pplist)
      self.usr_lflags = os.environ.get("USR_LFLAGS", self.usr_lflags)
      self.ifmpi = os.environ.get("MPI", self.ifmpi)

      # Get paths from env
      try:
         self.neklab_source_root = os.path.abspath(os.environ["NEKLAB_SOURCE_ROOT"])
      except KeyError:
         print('NEKLAB_SOURCE_ROOT nto defined')
         sys.exit(1)
      else:
         self.makeneklab = os.path.join(self.neklab_source_root, "app", "makeneklab")

      try:
         self.nek_source_root = os.path.abspath(os.environ["NEK_SOURCE_ROOT"])
      except KeyError:
         pass
      else:
         self.makenek = os.path.join(self.nek_source_root, "bin", "makenek")
         self.tools_root = os.path.join(self.nek_source_root, "tools")
         self.tools_bin = os.path.join(self.nek_source_root, "bin")

      self.examples_root = os.path.abspath(
         os.environ.get("EXAMPLES_ROOT", self.examples_root)
      )
      self.tools_root = os.path.abspath(os.environ.get("TOOLS_ROOT", self.tools_root))
      self.tools_bin = os.path.abspath(os.environ.get("TOOLS_BIN", self.tools_bin))

      try:
         self.log_root = os.path.abspath(os.environ["LOG_ROOT"])
      except KeyError:
         pass

      self.verbose = (
         str(os.environ.get("VERBOSE_TESTS", self.verbose)).lower() == "true"
      )
      self.parallel_procs = int(os.environ.get("PARALLEL_PROCS", self.parallel_procs))

      # Print everything out
      for varname, varval in (
         ("FC", self.f77),
         ("CC", self.cc),
         ("PPLIST", self.pplist),
         ("USR_LFLAGS", self.usr_lflags),
         ("IFMPI", self.ifmpi),
         ("NEKLAB_SOURCE_ROOT", self.neklab_source_root),
         ("NEK_SOURCE_ROOT", self.nek_source_root),
         ("EXAMPLES_ROOT", self.examples_root),
         ("LOG_ROOT", self.log_root),
         ("TOOLS_ROOT", self.tools_root),
         ("TOOLS_BIN", self.tools_bin),
         ("VERBOSE_TESTS", self.verbose),
         ("PARALLEL_PROCS", self.parallel_procs),
      ):
         if varval:
               print(f'    Using {varname:14} = "{varval}"')

      # Verify that pathnames are valid
      for varname, varval in (
         ("NEKLAB_SOURCE_ROOT", self.neklab_source_root),
         ("NEK_SOURCE_ROOT", self.nek_source_root),
         ("EXAMPLES_ROOT", self.examples_root),
         ("LOG_ROOT", self.log_root),
         ("TOOLS_ROOT", self.tools_root),
         ("TOOLS_BIN", self.tools_bin),
      ):
         if varval and not os.path.isdir(varval):
               raise OSError(
                  'The {0} directory "{1}" does not exist. Please the env variable ${0} to a valid directory.'.format(
                     varname, varval
                  )
               )

      print("Finished getting setup options!")

   def build_tools(
         self,
         targets=None,
         tools_root=None,
         tools_bin=None,
         f77=None,
         cc=None,
         bigmem=None,
         verbose=None,
      ):
         from lib.neklabBinBuild import build_tools

         build_tools(
            targets=targets if targets else ("clean", "genmap"),
            tools_root=tools_root if tools_root else self.tools_root,
            tools_bin=tools_bin if tools_bin else self.tools_bin,
            f77=f77,
            cc=cc,
            bigmem=bigmem if bigmem else "false",
            verbose=verbose if verbose else self.verbose,
         )

   def config_size(self, params=None, infile=None, outfile=None):
      from lib.neklabFileConfig import config_size

      cls = self.__class__

      if not infile:
         infile = os.path.join(self.nek_source_root, "core", "SIZE.template")
      if not outfile:
         outfile = os.path.join(self.examples_root, cls.example_subdir, "SIZE")
      if not params:
         params = self.size_params

      config_size(params=params, infile=infile, outfile=outfile)

   def mkSIZE(self, case=None):
      cls = self.__class__

      if not case:
         case = cls.case_name

      workdir = os.path.join(self.examples_root, cls.example_subdir)
      os.system("cd " + workdir + " ; " + self.source_root + "/bin/mkSIZE " + case)

   def config_parfile(self, opts=None, infile=None, outfile=None):
      from lib.neklabFileConfig import config_parfile

      cls = self.__class__

      if not infile:
         infile = os.path.join(
               self.examples_root, cls.example_subdir, cls.case_name + ".par"
         )
      if not outfile:
         outfile = infile
      if not opts:
         opts = {}

      config_parfile(opts=opts, infile=infile, outfile=outfile)

   def run_genmap(self, rea_file=None, tol="0.5"):

      from lib.neklabBinRun import run_meshgen

      cls = self.__class__

      if not rea_file:
         rea_file = cls.case_name

      run_meshgen(
         command=os.path.join(self.tools_bin, "genmap"),
         stdin=[rea_file, tol],
         cwd=os.path.join(self.examples_root, cls.example_subdir),
         verbose=self.verbose,
      )

   def run_genbox(self, box_file=None):
      from lib.neklabBinRun import run_meshgen

      if not box_file:
         box_file = self.__class__.case_name

      # Fix extension, in case user doesn't provide it
      root, ext = os.path.splitext(box_file)
      if ext != ".box":
         box_file = root + ext + ".box"

      run_meshgen(
         command=os.path.join(self.tools_bin, "genbox"),
         stdin=[box_file],
         cwd=os.path.join(self.examples_root, self.__class__.example_subdir),
         verbose=self.verbose,
      )

   def run_n2to3(self, stdin):
      from lib.neklabBinRun import run_meshgen

      run_meshgen(
         command=os.path.join(self.tools_bin, "n2to3"),
         stdin=stdin,
         cwd=os.path.join(self.examples_root, self.__class__.example_subdir),
      )

   def run_gmsh2nek(self, dim="3", fmsh_file=None, smsh_file=None, out_file=None, fP_list=None, sP_list=None):
      from lib.neklabBinRun import run_meshgen

      if not fmsh_file:
         fmsh_file = self.__class__.case_name

      stdin = [dim, fmsh_file]

      if not smsh_file:
         ifCHT = '0'
         stdin.append(ifCHT)
      else:
         ifCHT = '1'
         stdin.append(ifCHT)
         stdin.append(smsh_file)

      if not fP_list:
         P_pairs = '0'
         stdin.append(P_pairs)
      else:
         P_pairs = str(len(fP_list))
         stdin.append(P_pairs)
         stdin.extend(fP_list)

      if sP_list:
         P_pairs = str(len(sP_list))
         stdin.append(P_pairs)
         stdin.extend(sP_list)

      if not out_file:
         out_file = self.__class__.case_name

      stdin.append(out_file) 

      run_meshgen(
         command=os.path.join(self.tools_bin, "gmsh2nek"),
         stdin=stdin,
         cwd=os.path.join(self.examples_root, self.__class__.example_subdir),
         verbose=self.verbose,
      )

   def build_neklab(self, opts=None, usr_file=None):
      from lib.neklabBinBuild import build_neklab

      cls = self.__class__

      if not usr_file:
         usr_file = cls.case_name

      all_opts = dict(
         FC=self.f77,
         CC=self.cc,
         PPLIST=self.pplist,
         USR_LFLAGS=self.usr_lflags,
         MPI=int(self.ifmpi),
      )
      if opts:
         all_opts.update(opts)

      print(self.neklab_source_root)
      build_neklab(
         source_root=self.neklab_source_root,
         usr_file=usr_file,
         cwd=os.path.join(self.examples_root, cls.example_subdir),
         opts=all_opts,
         verbose=self.verbose,
      )

   def run_nek(self, rea_file=None, step_limit=None):
      from lib.neklabBinRun import run_nek

      cls = self.__class__
      run_nek(
         cwd=os.path.join(self.examples_root, cls.example_subdir),
         rea_file=cls.case_name if not rea_file else rea_file,
         ifmpi=self.ifmpi,
         log_suffix=self.log_suffix,
         n_procs=self.mpi_procs,
         step_limit=step_limit,
         verbose=self.verbose,
      )

   def mvn(self, src_prefix, dest_prefix):
      from lib.neklabBinRun import mvn

      cls = self.__class__
      mvn(
         src_prefix,
         dest_prefix,
         cwd=os.path.join(self.examples_root, cls.example_subdir),
      )

   def get_converged_eigs_data(self, logfile=None):
      cls = self.__class__
      if not logfile:
         logfile = os.path.join(
               self.examples_root,
               cls.example_subdir,
               f"eigs_output.txt",
         )
      eigs = {}

      with open(logfile) as file:
        for i, line in enumerate(file):
            parts = line.split()

            # Skip header or malformed lines
            if len(parts) < 6:
               continue
            try:
               re_part = float(parts[1])
               im_part = float(parts[2])
               modulus = float(parts[3])
               residual = float(parts[4])
               conv_flag = parts[5]
            except ValueError:
               # Skip header or bad lines
               continue

            if conv_flag == "T":
               label = f'lambda_{len(eigs)+1:d}'
               eigs[label] = {
                  "Re": re_part,
                  "Im": im_part,
                  "eigenvalue": complex(re_part, im_part),
                  "modulus": modulus,
                  "residual": residual,
                  "converged": True,
               }

      if not eigs:
         raise ValueError(
               f"No converged eigenvalues found in logfile '{logfile}'."
         )
      return eigs         

   def get_value_from_log(self, label, column, row=0, logfile=None):
      cls = self.__class__
      if not logfile:
         logfile = os.path.join(
               self.examples_root,
               cls.example_subdir,
               f"{cls.case_name}.log.{self.mpi_procs}{self.log_suffix}",
         )
      # Get all lines with label
      with open(logfile) as file:
         line_list = [l for l in file if label in l]
      if not line_list:
         raise ValueError(
               f'Could not find label "{label}" in logfile "{logfile}".  '
               "The run may have failed."
         )
      try:
         value = float(line_list[row].split()[column])
      except ValueError:
         raise ValueError(
               f'Attempted to parse non-numerical value in logfile, "{logfile}".'
               "  Logfile may be malformatted or run may have failed"
         )
      except IndexError:
         raise IndexError(
               f'Fewer rows and/or columns than expected in logfile, "{logfile}".'
               "  Logfile may be malformmated or run may have failed."
         )
      else:
         return value

   def get_phrase_from_log(self, label, logfile=None, row=0):
      cls = self.__class__
      if not logfile:
         logfile = os.path.join(
               self.examples_root,
               cls.example_subdir,
               f"{cls.case_name}.log.{self.mpi_procs}{self.log_suffix}",
         )

      with open(logfile) as file:
         line_list = [l for l in file if label in l]

      try:
         line = line_list[row]
      except IndexError:
         return None
      else:
         return line