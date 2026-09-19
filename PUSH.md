# Pushing this repo to GitHub

The repo is committed and clean: 8 commits on `main`, nothing uncommitted. GitHub is
unreachable from the sandbox this was built in, so these commands run on your machine,
under your own credentials. Nothing was pushed for you.

## Already done for you

- **Author identity fixed.** The commits were originally authored as
  `tilak@users.noreply.github.com`, which belongs to nobody — all 8 would have shown on
  GitHub as unlinked, and none would have counted toward your contribution graph. They are
  now `Tilak Sharma <TilakSharma07@users.noreply.github.com>`. Verify:

  ```bash
  cd ~/Downloads/ont-amr-nf
  git log --format='%an <%ae>' | sort -u
  ```

  Expected, exactly one line: `Tilak Sharma <TilakSharma07@users.noreply.github.com>`

- **Branch renamed** `master` → `main`.

- **`example_results/` committed** (188 KB). This is what makes the repo reviewable
  without installing anything — see below.

## 1. Create the empty repo

On GitHub, **New repository** → name `ont-amr-nf` → **Public** → do *not* add a README,
.gitignore or licence (the repo already has them; an initialising commit would force you
to merge).

Suggested description:

> Nextflow ONT AMR calling pipeline with two control layers and mutation-tested gates

## 2. Push

```bash
cd ~/Downloads/ont-amr-nf
git remote add origin https://github.com/TilakSharma07/ont-amr-nf.git
git push -u origin main
```

If it asks for a password, use a personal access token, not your account password
(<https://github.com/settings/tokens> → classic → scope `repo`).

## 3. Verify it is reviewable without tools

This is the point of `example_results/`. On a fresh clone, with no conda env, no
Nextflow, and no AMRFinderPlus installed:

```bash
git clone https://github.com/TilakSharma07/ont-amr-nf && cd ont-amr-nf
python3 tests/run_all.py example_results/results_main
```

Expected:

```
  [ok  ] test_module_paths.py         all 6 checks passed
  [skip] test_samplesheet.py          SKIP  nextflow cannot launch: ...
  [ok  ] test_validation_gate.py      all 6 checks passed
  [ok  ] test_control_gate.py         all 10 checks passed
  [ok  ] test_aggregate.py            all 5 checks passed
  [ok  ] test_figures.py              all 9 checks passed

  1 suite(s) skipped: test_samplesheet.py
all 6 suites passed (1 skipped)
```

36 checks run with nothing installed. The 5 skipped ones drive the real workflow with
malformed samplesheets, so they need a working Nextflow; the skip states that reason
rather than passing silently.

You can also re-render all three figures from the committed tables:

```bash
python3 bin/make_figures.py example_results/results_main /tmp/fig
python3 bin/make_figures.py example_results/results_titration /tmp/fig
```

## Note on the external drive

The full run outputs — assemblies, BAMs, the read cache — were on
`/media/tilak/LINUX MINT`, which became unreadable near the end of the session. They are
not in this repo and are not needed for the tests or figures. The conda environment at
`/tmp/ont-env` (AMRFinderPlus + database, Flye, minimap2, samtools) was cleared by a
workspace sweep; `environment.yml` rebuilds it when you want to run the pipeline itself
rather than check it.