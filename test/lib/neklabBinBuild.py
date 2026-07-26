import os, re
from subprocess import Popen, PIPE, STDOUT, DEVNULL
from pathlib import Path

def build_tools(
    tools_root,
    tools_bin,
    f77=None,
    cc=None,
    bigmem=None,
    targets=("clean", "all"),
    verbose=False,
):

    tools_root = Path(tools_root)

    print("Compiling tools... ")
    print(f'    Using output directory "{tools_bin}"')
    print(f'    Using FC "{f77}"')
    print(f'    Using CC "{cc}"')

    maketools_in = tools_root / "maketools"

    my_env = os.environ.copy()
    if f77:
        my_env["FC"] = f77
    if cc:
        my_env["CC"] = cc
    my_env["bin_nek_tools"] = tools_bin

    if targets[0] == "all":
        targets = [t for t in os.listdir(tools_root) if "maketools" not in t]
        print("Targets:", targets)

    for t in targets:
        proc = Popen([maketools_in, t], env=my_env, cwd=tools_root, stderr=STDOUT)
        proc.wait()
        logfile = tools_root / t / "build.log"
        if proc.returncode != 0:
            with open(logfile, "r") as file:
                text = file.read()
            print(text)
            exit(-1)

_OBJ = re.compile(r'-o\s+(\S+\.o)\b')
_ERR = re.compile(r'[Ee]rror|FATAL|make(\[|\:)')

def _stream_filtered(proc, logfile):
    with open(logfile, "w") as log:
        for line in proc.stdout:
            log.write(line)          # raw log: everything
            log.flush()
            m = _OBJ.search(line)    # stdout: filtered only
            if m:
                print(os.path.basename(m.group(1)))
            elif _ERR.search(line):
                print(line.rstrip())
    proc.wait()
    return None

def build_nek(source_root, usr_file, cwd=None, opts=None, verbose=False):

    if not opts:
        _opts = {}
    else:
        _opts = opts.copy()
    _opts.update(NEK_SOURCE_ROOT=source_root)

    print("Compiling nek5000...")
    print(f'    Using working directory "{cwd}"')
    print(f'    Using .usr file "{usr_file}"')
    for key, val in list(_opts.items()):
        print(f'    Using {key}="{val}"')

    my_env = os.environ.copy()
    if source_root:
        my_env["NEK_SOURCE_ROOT"] = source_root
    if _opts.get("F77"):
        my_env["FC"] = _opts.get("F77")
    if _opts.get("CC"):
        my_env["CC"] = _opts.get("CC")
    if _opts.get("PPLIST"):
        my_env["PPLIST"] = _opts.get("PPLIST")

    makenek_in = Path(source_root) / "bin" / "makenek"
    logfile = Path(cwd) / "build.log"

    proc = Popen([makenek_in, "clean"], cwd=cwd, env=my_env, stdin=PIPE, stdout=PIPE, stderr=STDOUT, text=True)
    proc.communicate(input="Y\n")
    proc.wait()

    proc = Popen([makenek_in, usr_file], cwd=cwd, env=my_env, stdin=PIPE, stdout=PIPE, stderr=STDOUT, text=True)
    out, _ = proc.communicate(input="\n")

    proc.wait()

    if proc.returncode != 0:
        print(out)
        with open(logfile, "w") as f:
            f.write(out)
        exit(-1)

def build_neklab(source_root, usr_file, cwd=None, opts=None, verbose=False):

    if not opts:
        _opts = {}
    else:
        _opts = opts.copy()
    _opts.update(NEKLAB_SOURCE_ROOT=source_root)

    print("Compiling neklab...")
    print(f'    Using working directory "{cwd}"')
    print(f'    Using .usr file "{usr_file}"')
    for key, val in list(_opts.items()):
        print(f'    Using {key}="{val}"')

    my_env = os.environ.copy()
    if source_root:
        my_env["NEKLAB_SOURCE_ROOT"] = source_root
    if _opts.get("F77"):
        my_env["FC"] = _opts.get("F77")
    if _opts.get("CC"):
        my_env["CC"] = _opts.get("CC")
    if _opts.get("PPLIST"):
        my_env["PPLIST"] = _opts.get("PPLIST")
    if "MPI" in _opts:
        my_env["MPI"] = str(int(_opts["MPI"]))   # force 0/1, not "true"

    makeneklab_in = Path(source_root) / "app" / "makeneklab"
    logfile = Path(cwd) / "build.log"

    cwd_p = Path(cwd)
    # Recover from a crashed prior run of makeneklab, which does:
    #   mv <case>.usr <case>.usr_tmp ; cp dummy ; ... ; cp <case>.usr_tmp <case>.usr
    # A crash can leave <case>.usr_tmp without <case>.usr, and/or a stray nek_dummy.usr.
    for tmp in cwd_p.glob("*.usr_tmp"):
        real = tmp.with_suffix("")          # cav.usr_tmp -> cav.usr
        real = tmp.parent / tmp.name[:-4]   # strip "_tmp" -> cav.usr
        if not real.exists():
            # swap was interrupted; the tmp IS the real case file — restore it
            tmp.rename(real)
        else:
            # both exist: real is intact, tmp is a stale leftover — safe to drop
            tmp.unlink()

    # Only ever a build artifact; never the real case file.
    dummy = cwd_p / "nek_dummy.usr"
    if dummy.exists():
        dummy.unlink()

    proc = Popen([makeneklab_in, "clean"], cwd=cwd, env=my_env,
             stdin=PIPE, text=True)
    out, _ = proc.communicate(input="Y\n")

    proc = Popen([makeneklab_in], cwd=cwd, env=my_env,
                stdin=DEVNULL, text=True)
    proc.wait()

    if proc.returncode != 0:
        exit(-1)
    