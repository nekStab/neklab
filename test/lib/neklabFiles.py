"""Small, targeted edits of case files, and reading results from a log.

Every edit changes only the values it is asked to change; all other lines,
including comments, stay exactly as they are.
"""
import re


def _fortran_value(value):
   """Python value -> Fortran literal. Floats get the _dp kind; strings are used as written."""
   if isinstance(value, bool):
      return ".true." if value else ".false."
   if isinstance(value, float):
      return f"{value!r}_dp"
   return str(value)


def _value_end(text, i):
   """End of the Fortran value starting at text[i]: the first top-level ',' or
   ')', a '!' comment or the end of the line (strings and brackets skipped)."""
   depth, quote = 0, None
   while i < len(text) and text[i] != "\n":
      c = text[i]
      if quote:
         if c == quote:
            quote = None
      elif c in "'\"":
         quote = c
      elif c == "(":
         depth += 1
      elif c == ")":
         if depth == 0:
            break
         depth -= 1
      elif c == "," and depth == 0 or c == "!":
         break
      i += 1
   return i


def set_fortran_parameters(path, params):
   """Change the value of Fortran parameters or initialised variables in a file.

   Works for the declaration forms
         parameter (lpmin=8)                                 ! SIZE style
         real(dp), parameter :: tau = 0.05_dp                ! parameter
         real(dp), public    :: tau = 0.1_dp                 ! module variable
         real(dp), public    :: omega = 8.0_dp*atan(1.0_dp)  ! expressions too
   Only the value is replaced; comments and spacing are kept. Assignments and
   keyword arguments (no `::` or `parameter` before the name) are not touched.
   Raises KeyError if a name is not found, so a typo cannot go unnoticed.
   """
   with open(path) as f:
      text = f.read()
   for name, value in params.items():
      pattern = re.compile(
         rf"^[^!\n]*(?:\bparameter\b|::)[^!\n]*?\b{re.escape(name)}\s*=(?!=)\s*",
         re.IGNORECASE | re.MULTILINE,
      )
      new, pos, n = [], 0, 0
      for m in pattern.finditer(text):
         end = _value_end(text, m.end())
         old = text[m.end():end]
         trail = old[len(old.rstrip()):]          # keep the spaces before a comment
         new += [text[pos:m.end()], _fortran_value(value), trail]
         pos, n = end, n + 1
      if n == 0:
         raise KeyError(f'"{name}" is not declared with a value in {path}')
      text = "".join(new) + text[pos:]
   with open(path, "w") as f:
      f.write(text)


def set_par_values(path, values):
   """Change `key = value` lines in a .par file: {"GENERAL": {"endTime": 0.05}}.

   Keys and sections are matched case-insensitively, inline comments are kept.
   A key that is missing is added at the end of its section (the section is
   created if needed).
   """
   with open(path) as f:
      lines = f.read().splitlines()

   for section, entries in values.items():
      for key, value in entries.items():
         # locate the section
         start = next((i for i, l in enumerate(lines)
                       if l.strip().lower() == f"[{section.lower()}]"), None)
         if start is None:
            lines += ["", f"[{section.upper()}]"]
            start = len(lines) - 1
         end = next((i for i in range(start + 1, len(lines))
                     if lines[i].strip().startswith("[")), len(lines))
         # replace the key inside the section, or add it
         key_re = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)([^#]*?)(\s*#.*)?$", re.IGNORECASE)
         for i in range(start + 1, end):
            m = key_re.match(lines[i])
            if m:
               lines[i] = m.group(1) + str(value) + (m.group(3) or "")
               break
         else:
            insert_at = end
            while insert_at > start + 1 and not lines[insert_at - 1].strip():
               insert_at -= 1
            lines.insert(insert_at, f"{key} = {value}")

   with open(path, "w") as f:
      f.write("\n".join(lines) + "\n")


def get_fortran_parameter(path, name):
   """Value of a Fortran parameter as written in the file (string), or None."""
   with open(path) as f:
      m = re.search(rf"^[^!\n]*\bparameter\b[^!\n]*?\b{name}\s*=\s*([^,)\s!]+)",
                    f.read(), re.IGNORECASE | re.MULTILINE)
   return m.group(1) if m else None


_ANSI = re.compile(r"\x1b\[[0-9;]*m")
_CHECK_NEW = re.compile(                     # NEKLAB_CHECK <name> [tol = x] PASSED|FAILED [msg]
   r"NEKLAB_CHECK\s+(?P<name>.*?)\s+(?:tol\s*=\s*\S+\s+)?(?P<status>PASSED|FAILED)\b\s*(?P<msg>.*)$")
_CHECK_OLD = re.compile(                     # NEKLAB_CHECK PASS|FAIL <name> [msg]
   r"NEKLAB_CHECK\s+(?P<status>PASS|FAIL)\s+(?P<name>\S+)\s*(?P<msg>.*)$")


def read_checks(logfile):
   """NEKLAB_CHECK lines -> {name: ("PASS"|"FAIL", message)}.

   Both forms are read, anywhere in the line (log prefixes and colours are fine):
      NEKLAB_CHECK <name> [tol = 1.0E-06] PASSED|FAILED [message]
      NEKLAB_CHECK PASS|FAIL <name> [message]
   A check reported several times counts as failed if it failed once.
   """
   checks = {}
   with open(logfile, errors="replace") as f:
      for line in f:
         line = _ANSI.sub("", line).rstrip()
         m = _CHECK_NEW.search(line) or _CHECK_OLD.search(line)
         if not m:
            continue
         status = "PASS" if m["status"].startswith("PASS") else "FAIL"
         name, message = m["name"].strip(), m["msg"].strip()
         if checks.get(name, ("PASS",))[0] != "FAIL":
            checks[name] = (status, message if status == "FAIL" else "")
   return checks


def log_contains(logfile, text):
   with open(logfile, errors="replace") as f:
      return any(text in _ANSI.sub("", line) for line in f)


def read_results(logfile):
   """Lines 'NEKLAB_RESULT <key> <value>' -> {key: float}. The last occurrence wins."""
   results = {}
   with open(logfile, errors="replace") as f:
      for line in f:
         line = _ANSI.sub("", line)
         if "NEKLAB_RESULT" not in line:
            continue
         parts = line[line.index("NEKLAB_RESULT"):].split()
         if len(parts) >= 3:
            try:
               results[parts[1]] = float(parts[2].replace("D", "E"))
            except ValueError:
               results[parts[1]] = parts[2]
   return results


def read_eigenvalues(path):
   """Converged eigenvalues from neklab's eigs_output.txt, in file order.

   Returns a list of dicts with keys re, im, modulus, residual.
   """
   eigs = []
   with open(path) as f:
      for line in f:
         parts = line.split()
         if len(parts) < 6 or parts[5] != "T":
            continue
         try:
            re_, im, mod, res = (float(x) for x in parts[1:5])
         except ValueError:
            continue
         eigs.append(dict(re=re_, im=im, modulus=mod, residual=res))
   return eigs
