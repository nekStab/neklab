"""Reading and editing case files, and reading results from a log.

Edits change only the values they are asked to change; all other lines,
including comments, stay exactly as they are.
"""
import re


# --- SIZE ----------------------------------------------------------------------

_SIZE_PARAM = r"^([^!\n]*\bparameter\s*\([^!\n]*?\b{name}\s*=\s*)([^,)\s!]+)"


def read_size(path):
   """parameter (name=value) declarations of a SIZE file -> {name: value as written}."""
   values = {}
   for line in open(path):
      code = line.split("!")[0]
      m = re.search(r"\bparameter\s*\((.*)\)", code, re.IGNORECASE)
      if m:
         for item in m.group(1).split(","):
            if "=" in item:
               name, value = item.split("=", 1)
               values[name.strip().lower()] = value.strip()
   return values


def set_size_parameters(path, params):
   """Change parameter (name=value) values in SIZE; comments are kept.

   Raises KeyError for an unknown name. Returns [(name, old, new)].
   """
   text = open(path).read()
   changes = []
   for name, value in params.items():
      pattern = re.compile(_SIZE_PARAM.format(name=re.escape(name)), re.IGNORECASE | re.MULTILINE)
      m = pattern.search(text)
      if not m:
         raise KeyError(f'"{name}" is not a parameter in {path}')
      changes.append((name, m.group(2), str(value)))
      text = pattern.sub(lambda m: m.group(1) + str(value), text, count=1)
   open(path, "w").write(text)
   return changes


# --- .par ----------------------------------------------------------------------

def read_par(path):
   """A .par file -> {SECTION: {key: value}}, comments removed, order kept."""
   sections, current = {}, None
   for line in open(path):
      code = line.split("#")[0].strip()
      if code.startswith("[") and code.endswith("]"):
         current = sections.setdefault(code[1:-1].upper(), {})
      elif "=" in code and current is not None:
         key, value = code.split("=", 1)
         current[key.strip()] = value.strip()
   return sections


def set_par_values(path, values):
   """Change `key = value` lines in a .par file: {"GENERAL": {"endTime": 0.05}}.

   Keys and sections are matched case-insensitively, inline comments are kept.
   A key that is missing is added at the end of its section (the section is
   created if needed). Returns the changes as a list of
   (section, key, old value or None if added, new value).
   """
   with open(path) as f:
      lines = f.read().splitlines()
   changes = []

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
               changes.append((section, key, m.group(2).strip(), str(value)))
               lines[i] = m.group(1) + str(value) + (m.group(3) or "")
               break
         else:
            insert_at = end
            while insert_at > start + 1 and not lines[insert_at - 1].strip():
               insert_at -= 1
            lines.insert(insert_at, f"{key} = {value}")
            changes.append((section, key, None, str(value)))

   with open(path, "w") as f:
      f.write("\n".join(lines) + "\n")
   return changes


# --- settings (initialised variables of the .usr) -------------------------------
#
# A setting is a declaration with an initial value, e.g. at the top of userchk
#    real(dp) :: tau = 0.1_dp
# `parameter`s are never settings. Names must be unique in the file: the first
# declaration with the name is the one that is read and changed. In a file with
# a module, only the part before the first `contains` is searched.

def _settings_part(text):
   """The text before the first `contains` (all of it if there is none)."""
   m = re.search(r"^\s*contains\b", text, re.IGNORECASE | re.MULTILINE)
   return text[:m.start()] if m else text


_SETTING = r"^([^!\n]*::\s*{name}\s*=\s*)(.*?)(\s*(?:!.*)?)$"


def read_settings(path):
   """Values of the settings in a .usr (or module) -> {name: value as written}."""
   values = {}
   for line in _settings_part(open(path).read()).splitlines():
      if line[:1] in "cC*" or "parameter" in line.lower():
         continue
      m = re.match(r"^[^!]*::\s*(\w+)\s*=\s*(.*?)\s*(?:!.*)?$", line)
      if m:
         values[m.group(1).lower()] = m.group(2)
   return values


def fortran_value(value):
   """Python value -> Fortran literal: floats get _dp, strings are quoted."""
   if isinstance(value, bool):
      return ".true." if value else ".false."
   if isinstance(value, float):
      return f"{value!r}_dp"
   if isinstance(value, str):
      return "'" + value.replace("'", "''") + "'"
   return str(value)


def set_settings(path, values):
   """Change the initial value of settings in a .usr (or module); comments are kept.

   Raises KeyError for a name that is not a setting. Returns [(name, old, new)].
   """
   text = open(path).read()
   part = _settings_part(text)
   rest = text[len(part):]
   changes = []
   for name, value in values.items():
      pattern = re.compile(_SETTING.format(name=re.escape(name)), re.IGNORECASE | re.MULTILINE)
      m = pattern.search(part)
      if not m:
         raise KeyError(f'"{name}" is not a setting in {path}; settings: {sorted(read_settings(path))}')
      new = fortran_value(value)
      changes.append((name, m.group(2), new))
      part = part[:m.start(2)] + new + part[m.end(2):]
   open(path, "w").write(part + rest)
   return changes


def same_value(a, b):
   """True if two values (as written in Fortran/.par or Python) are the same setting."""
   def norm(v):
      if isinstance(v, bool):
         return v
      v = str(v).strip().strip("'\"").lower()
      if v in (".true.", "t", "true", "yes"):
         return True
      if v in (".false.", "f", "false", "no"):
         return False
      try:
         return float(re.sub(r"_dp$", "", v).replace("d", "e"))
      except ValueError:
         return v
   return norm(a) == norm(b)


# --- the log -----------------------------------------------------------------------

_ANSI = re.compile(r"\x1b\[[0-9;]*m")
_CHECK = re.compile(                         # NEKLAB_CHECK <name> [tol = x] PASSED|FAILED [msg]
   r"NEKLAB_CHECK\s+(?P<name>.*?)\s+(?:tol\s*=\s*\S+\s+)?(?P<status>PASSED|FAILED)\b\s*(?P<msg>.*)$")


def read_checks(logfile):
   """NEKLAB_CHECK lines -> {name: ("PASS"|"FAIL", message)}.

   Read anywhere in the line (log prefixes and colours are fine):
      NEKLAB_CHECK <name> [tol = 1.0E-06] PASSED|FAILED [message]
   A check reported several times counts as failed if it failed once.
   """
   checks = {}                 # in the order they first appear in the log
   with open(logfile, errors="replace") as f:
      for line in f:
         line = _ANSI.sub("", line).rstrip()
         m = _CHECK.search(line)
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


def find_lines(logfile, markers):
   """Every line containing one of the markers, from the marker on, duplicates removed."""
   found = []
   with open(logfile, errors="replace") as f:
      for line in f:
         line = _ANSI.sub("", line).rstrip()
         for marker in markers:
            if marker in line:
               text = line[line.index(marker):]
               if text not in found:
                  found.append(text)
               break
   return found


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
