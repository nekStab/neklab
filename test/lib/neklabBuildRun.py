"""Building and running a Nek5000/neklab case.

Every function raises an exception on failure, with the end of the relevant
log in the message, so a broken case makes its own test fail and the rest of
the suite continues.
"""
import os
import subprocess
from pathlib import Path


class BuildError(RuntimeError):
   pass


class RunError(RuntimeError):
   pass


def tail(path, n=30):
   try:
      with open(path) as f:
         return "".join(f.readlines()[-n:])
   except OSError:
      return f"(no file {path})"


def _run(command, cwd, env, logfile, stdin="", verbose=False, append=False):
   """Run a command, send its output to logfile (and to the screen if verbose)."""
   with open(logfile, "a" if append else "w") as log:
      proc = subprocess.Popen([str(c) for c in command], cwd=cwd, env=env, text=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
      proc.stdin.write(stdin)
      proc.stdin.close()
      for line in proc.stdout:
         log.write(line)
         if verbose:
            print(line, end="")
      return proc.wait()


_TOOLS_BUILT = set()


def build_tool(tool, nek_source_root, verbose=False):
   """Build a Nek5000 tool (e.g. genmap) once per test session."""
   if tool in _TOOLS_BUILT:
      return
   tools = Path(nek_source_root) / "tools"
   env = dict(os.environ, bin_nek_tools=str(Path(nek_source_root) / "bin"))
   rc = _run([tools / "maketools", tool], tools, env, tools / f"{tool}.build.out",
             verbose=verbose)
   if rc != 0:
      raise BuildError(f"maketools {tool} failed:\n{tail(tools / tool / 'build.log')}")
   _TOOLS_BUILT.add(tool)


def run_genmap(case_dir, case_name, nek_source_root, tolerance="0.01", verbose=False):
   """Create <case>.ma2 from <case>.re2."""
   build_tool("genmap", nek_source_root, verbose)
   genmap = Path(nek_source_root) / "bin" / "genmap"
   log = Path(case_dir) / "genmap.out"
   rc = _run([genmap], case_dir, os.environ.copy(), log, stdin=f"{case_name}\n{tolerance}\n",
             verbose=verbose)
   if rc != 0 or not (Path(case_dir) / f"{case_name}.ma2").exists():
      raise BuildError(f"genmap failed:\n{tail(log)}")


def build_neklab(case_dir, neklab_source_root, nek_source_root, verbose=False):
   """Clean build of the case with makeneklab (the script finds the .usr itself)."""
   env = dict(os.environ, NEKLAB_SOURCE_ROOT=str(neklab_source_root),
              NEK_SOURCE_ROOT=str(nek_source_root))
   makeneklab = Path(neklab_source_root) / "app" / "makeneklab"
   log = Path(case_dir) / "build.out"
   _run([makeneklab, "clean"], case_dir, env, log, stdin="Y\n", verbose=verbose)
   rc = _run([makeneklab], case_dir, env, log, verbose=verbose, append=True)
   if rc != 0 or not (Path(case_dir) / "nek5000").exists():
      raise BuildError(f"makeneklab failed (code {rc}):\n{tail(log)}")


def run_nek(case_dir, case_name, nprocs, mpiexec="mpiexec", verbose=False):
   """Run nek5000 on nprocs MPI ranks; returns the logfile path."""
   case_dir = Path(case_dir)
   (case_dir / "SESSION.NAME").write_text(f"1\n{case_name}\n{case_dir}/\n")
   (case_dir / "ioinfo").unlink(missing_ok=True)    # a leftover would stop the run early
   log = case_dir / f"{case_name}.log.{nprocs}"
   rc = _run([mpiexec, "-np", nprocs, case_dir / "nek5000"], case_dir, os.environ.copy(),
             log, verbose=verbose)
   if rc != 0:
      raise RunError(f"nek5000 exited with code {rc}. End of {log}:\n{tail(log)}")
   return log
