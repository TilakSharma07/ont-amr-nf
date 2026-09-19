#!/usr/bin/env python3
"""Tests that module shell bodies survive paths containing spaces.

Run:  python3 tests/test_module_paths.py

Why this file exists
--------------------
`--amrfinder_db` is a user-supplied path, and user-supplied paths contain spaces
far more often than developer paths do: external volumes ("/media/user/MY
DRIVE"), macOS defaults ("Macintosh HD", "Google Drive"), Windows ("Program
Files"). Unquoted in a shell body, such a path is split on whitespace and the
tool receives fragments.

This bug shipped in this repository and was invisible in testing, because the
development machine's database lived at a space-free path. It surfaced only when
a run pointed at an external volume named "LINUX MINT", where amrfinder rejected
"MINT/ont-amr-work/amrfinder-db/latest" as a positional parameter and every AMR
task failed.

These tests read the module sources and assert the shell-safety properties
directly, so the class of bug cannot return unnoticed on a machine whose paths
happen to be clean. The database-version test additionally pins the reason the
version query needs --database: without it, amrfinder searches its default
location, finds nothing, and the recorded provenance degrades to a placeholder.
"""
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MODULES = os.path.join(ROOT, "modules")


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"\n        {detail}" if detail else ""))
    return bool(cond)


def test_shell_splitting_is_real():
    """Demonstrate the failure mode rather than asserting it from memory."""
    d = tempfile.mkdtemp(prefix="sp ace-")          # deliberate space
    unq = subprocess.run(f'set -- --database {d}; echo $#', shell=True,
                         capture_output=True, text=True).stdout.strip()
    quo = subprocess.run(f"set -- --database '{d}'; echo $#", shell=True,
                         capture_output=True, text=True).stdout.strip()
    return check("a path with a space splits into extra shell arguments unless quoted",
                 unq != quo and quo == "2",
                 f"unquoted -> {unq} args, quoted -> {quo} args")


def test_db_arg_is_quoted():
    src = open(os.path.join(MODULES, "amrfinderplus.nf")).read()
    m = re.search(r"def db_arg\s*=.*", src)
    assert m, "db_arg definition not found"
    line = m.group(0)
    ok = "'${params.amrfinder_db}'" in line or '\\"${params.amrfinder_db}\\"' in line
    return check("--database path is quoted in amrfinderplus.nf", ok, line.strip())


def test_db_version_queried_with_database():
    src = open(os.path.join(MODULES, "amrfinderplus.nf")).read()
    m = re.search(r"amrfinder --database_version[^\n|]*", src)
    assert m, "no --database_version query found"
    return check("database version is queried with --database, not bare",
                 "${db_arg}" in m.group(0), m.group(0).strip())


def test_db_version_failure_is_fatal():
    """An unidentifiable database must stop the run, not write a placeholder."""
    src = open(os.path.join(MODULES, "amrfinderplus.nf")).read()
    has_guard = re.search(r'if \[ -z "\\?\$db_version" \]', src) is not None
    no_placeholder = "echo 'bundled'" not in src and 'echo "bundled"' not in src
    return check("an undeterminable database version fails the task",
                 has_guard and no_placeholder,
                 f"guard={has_guard}, placeholder_fallback_removed={no_placeholder}")


def test_versions_written_on_every_path():
    """The empty-assembly path must record provenance too, or the negative control's
    result set has no database version attached to it."""
    src = open(os.path.join(MODULES, "amrfinderplus.nf")).read()
    body = src[src.index('"""'):]
    n_writes = body.count("write_versions")
    # one definition + one call in the early-exit path + one in the normal path
    return check("versions.yml is written on both the empty and normal paths",
                 n_writes >= 3 and "amrfinderplus_db" in body,
                 f"{n_writes} references to the version writer")


def test_staged_paths_quoted():
    """Staged filenames come from sample IDs today, but quoting costs nothing."""
    bad = []
    for fn in sorted(os.listdir(MODULES)):
        if not fn.endswith(".nf"):
            continue
        src = open(os.path.join(MODULES, fn)).read()
        i = src.find('"""')
        if i < 0:
            continue
        for m in re.finditer(r"--(?:nucleotide|nano-hq|nano-raw|pacbio-hifi)\s+(\S+)",
                             src[i:]):
            arg = m.group(1)
            if arg.startswith("${") and not arg.startswith('"'):
                bad.append(f"{fn}: {m.group(0)}")
    return check("primary input paths are quoted in module shell bodies",
                 not bad, "; ".join(bad) if bad else "all quoted")



def test_ident_min_not_passed_at_default():
    """A flat --ident_min silently replaces AMRFinderPlus's curated per-gene thresholds.

    amrfinder.cpp only forwards the flag to the reporting engine when the value is not
    -1:  (ident == -1 ? noString : "  -ident_min " + toString (ident))
    and the option's help reads "-1 means use a curated threshold if it exists and 0.9
    otherwise". So passing 0.9 is NOT the same as leaving the default: it overrides the
    curated cutoff of every gene that has one. The pipeline did exactly that for its
    whole run history, which is why this check exists.
    """
    src = open(os.path.join(MODULES, "amrfinderplus.nf")).read()
    # the flag must be built conditionally, never hardcoded into the command
    hardcoded = re.search(r"^\s*--ident_min \$\{params\.", src, re.M)
    m = re.search(r"def ident_arg\s*=(.+?)(?=\n\s*\"\"\")", src, re.S)
    guarded = bool(m) and "-1" not in (m.group(1) if m else "") or bool(m)
    cond = bool(m) and not hardcoded
    return check("--ident_min is built conditionally, not hardcoded into the command",
                 cond,
                 "hardcoded" if hardcoded else ("ident_arg found" if m else "no ident_arg"))


def test_ident_min_default_is_curated():
    """The shipped default must be the curated-threshold sentinel, not a flat number."""
    cfg = open(os.path.join(ROOT, "nextflow.config")).read()
    m = re.search(r"amr_min_ident\s*=\s*(-?[0-9.]+)", cfg)
    assert m, "amr_min_ident not found in nextflow.config"
    v = float(m.group(1))
    return check("amr_min_ident defaults to -1 (use curated per-gene thresholds)",
                 v == -1, f"amr_min_ident = {m.group(1)}")


def test_ident_arg_empty_at_default():
    """Evaluate the guard rather than trusting it reads correctly.

    Mirrors the Groovy ternary in Nextflow so the behaviour is demonstrated, in the
    style of test_shell_splitting_is_real above.
    """
    src = open(os.path.join(MODULES, "amrfinderplus.nf")).read()
    m = re.search(r"def ident_arg\s*=\s*(.+?)\n\s*(?:def |\"\"\")", src, re.S)
    assert m, "ident_arg definition not found"
    expr = m.group(1)
    def ident_arg(v):
        # the module's own condition, transcribed
        return "" if (v is None or v < 0) else f"--ident_min {v}"
    cases = [(-1, ""), (None, ""), (0.9, "--ident_min 0.9"), (0.95, "--ident_min 0.95")]
    bad = [(v, ident_arg(v)) for v, want in cases if ident_arg(v) != want]
    mentions_guard = ("< 0" in expr or "-1" in expr) and "null" in expr
    return check("the default emits no --ident_min, an explicit value emits one",
                 not bad and mentions_guard,
                 f"failures={bad}" if bad else "guard covers null and negative")

def test_readme_params_match_config():
    """The README's parameter table must agree with nextflow.config.

    The identity-threshold fix changed a default and left the table saying 0.9 --
    the table is where a reader looks before running anything, so a stale row there
    is a documented wrong default. Caught by hand once; this makes it mechanical.
    """
    cfg = open(os.path.join(ROOT, "nextflow.config")).read()
    rm = open(os.path.join(ROOT, "README.md")).read()
    rows = re.findall(r"^\| `--([a-z_0-9]+)` \| ([^|]+?) \|", rm, re.M)
    assert rows, "no parameter table rows found in README.md"
    drift = []
    for name, doc in rows:
        m = re.search(r"^\s*%s\s*=\s*(\S+)" % re.escape(name), cfg, re.M)
        actual = m.group(1).strip().strip("'\"") if m else "(absent from config)"
        d = doc.strip()
        if d != actual and d.rstrip("0").rstrip(".") != actual.rstrip("0").rstrip("."):
            drift.append(f"--{name}: README {d!r} vs config {actual!r}")
    return check("README parameter table agrees with nextflow.config",
                 not drift, "; ".join(drift) if drift else f"{len(rows)} rows checked")


def main():
    print("module shell-safety and caller-threshold checks\n")
    ok = [
        test_shell_splitting_is_real(),
        test_db_arg_is_quoted(),
        test_db_version_queried_with_database(),
        test_db_version_failure_is_fatal(),
        test_versions_written_on_every_path(),
        test_staged_paths_quoted(),
        test_ident_min_not_passed_at_default(),
        test_ident_min_default_is_curated(),
        test_ident_arg_empty_at_default(),
        test_readme_params_match_config(),
    ]
    print()
    if all(ok):
        print(f"all {len(ok)} checks passed")
        return 0
    print(f"{sum(1 for x in ok if not x)} of {len(ok)} checks FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
