# Pushing this repo to GitHub

The repo is committed and clean, on `main`, nothing uncommitted. GitHub is unreachable
from the sandbox this was built in, so these commands run on your machine, under your
own credentials. Nothing was pushed for you.

```bash
cd ~/Downloads/ont-amr-nf && git rev-list --count HEAD && git status --short
```

This file used to state the commit count. It was wrong by six, I corrected it, and the
correcting commit made it wrong again — twice, because fixing that made it wrong a
third time. A number that changes every time you touch the file cannot be maintained by
hand, so it is a command now. The same goes for the check counts in section 3: they are
there to tell you roughly what to expect, not to be trusted over the runner's output.

## Already done for you

- **Author identity fixed.** The commits were originally authored as
  `tilak@users.noreply.github.com`, which belongs to nobody — every one would have shown
  on GitHub as unlinked, and none would have counted toward your contribution graph. They are
  now `Tilak Sharma <TilakSharma07@users.noreply.github.com>`. Verify:

  ```bash
  cd ~/Downloads/ont-amr-nf
  git log --format='%an <%ae>' | sort -u
  ```

  Expected, exactly one line: `Tilak Sharma <TilakSharma07@users.noreply.github.com>`

- **Branch renamed** `master` → `main`.

- **`example_results/` committed** (188 KB). This is what makes the repo reviewable
  without installing anything — see below.

## 1. Push

The repo exists on GitHub and `origin` in this clone already points at it. Nothing to
create, nothing to add:

```bash
cd ~/Downloads/ont-amr-nf
git push -u origin main
```

If it asks for a username and password, `git` has no credentials for github.com yet.
`gh` does, so hand them over once and retry:

```bash
gh auth setup-git
git push -u origin main
```

Your GitHub web password will not work even if you type it correctly — it has not been
accepted for git over HTTPS for years. `gh auth setup-git` is the fix; a personal
access token used as the password also works
(<https://github.com/settings/tokens> -> classic -> scope `repo`).

Then verify the push landed, rather than trusting that the command printed nothing
alarming:

```bash
git log --oneline -1 origin/main    # must match: git log --oneline -1
```

### Why `gh repo create` left you with an empty repo

`gh repo create ont-amr-nf --public --source=. --remote=origin --push` creates and
pushes in one step, **but only in a clone that has no `origin` yet.** This clone has
one, so `--remote=origin` fails with `Unable to add remote "origin"`, and because that
step failed `--push` never ran. The output reports a success and a failure together and
leaves an empty repository on GitHub. An earlier version of this file recommended that
command and claimed `--source=.` would reuse the existing remote. It does not.

If that has already happened you are exactly where section 1 starts: just push. Do not
run the `git remote add origin ...` line GitHub shows on the empty-repo page either --
it will say `remote origin already exists`. Check before adding:

```bash
git remote -v
```

## 2. Verify it is reviewable without tools

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
  [ok  ] test_aggregate.py            all 6 checks passed
  [ok  ] test_figures.py              all 10 checks passed

  1 suite(s) skipped: test_samplesheet.py
all 6 suites passed (1 skipped)
```

38 checks run with nothing installed. The 5 skipped ones drive the real workflow with
malformed samplesheets, so they need a working Nextflow; the skip states that reason
rather than passing silently.

You can also re-render all three figures from the committed tables. Two commands are
needed because figure 3 is the depth titration, which the main run does not produce:

```bash
python3 bin/make_figures.py example_results/results_main /tmp/fig
python3 bin/make_figures.py example_results/results_titration /tmp/fig
```

The second command prints `skipping fig1` and `skipping fig2`. That is correct and it
matters: both write fixed filenames, so before those guards existed the second command
silently replaced the real six-sample figures with one-isolate, no-control versions of
themselves — and both renders reported success. `tests/test_figures.py` now fails if
that regresses.

## Note on the external drive

The full run outputs — assemblies, BAMs, the read cache — are on
`/media/tilak/LINUX MINT` (153 MB under `ont-amr-work/results_main`). The drive stopped
responding for part of the session but reads fine now. They are not in this repo and
are not needed for the tests or figures: the 188 KB of tables in `example_results/` is
everything the checks and figures consume.

The conda environment (AMRFinderPlus + database, Flye, minimap2, samtools) is archived
on the same drive as `ont-amr-work/ont-env.tar.gz`, 871 MB. It was originally unpacked
at `/tmp/ont-env` and disappeared — not a workspace sweep, as this file previously
claimed, but because `/tmp` is a tmpfs that is cleared between runs. Unpack it
somewhere persistent, or let `environment.yml` rebuild it, when you want to run the
pipeline itself rather than check it.