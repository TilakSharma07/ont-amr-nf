# Pushing this repo to GitHub

The repo is committed and clean; it just needs a remote. These commands run on your
machine, under your own credentials — nothing here was pushed for you.

## 1. Check the commit author email first

GitHub only links commits to your profile when the author email is one your account
owns. Check what is recorded:

```bash
cd ont-amr-nf
git log --format='%an <%ae>' | sort -u
```

Expected: `Tilak <tilak@users.noreply.github.com>`.

**If your GitHub username is not `tilak`**, that email belongs to nobody and the commits
will show as unlinked. Fix it before pushing — find your real noreply address at
<https://github.com/settings/emails> (it looks like `12345678+yourname@users.noreply.github.com`),
then:

```bash
git config user.email 'YOUR_REAL_NOREPLY@users.noreply.github.com'

FILTER_BRANCH_SQUELCH_WARNING=1 git filter-branch -f --env-filter '
  export GIT_AUTHOR_EMAIL="YOUR_REAL_NOREPLY@users.noreply.github.com"
  export GIT_COMMITTER_EMAIL="YOUR_REAL_NOREPLY@users.noreply.github.com"
' -- --all

git log --format='%an <%ae>' | sort -u   # verify
```

A backup of the history before any rewrite is at
`ont-amr-work/ont-amr-nf-backup.bundle` on the USB drive. To restore from it:
`git clone ont-amr-nf-backup.bundle recovered-repo`.

## 2. Create the repo and push

With the `gh` CLI:

```bash
gh repo create ont-amr-nf --public --source=. --remote=origin --push
```

Or by hand — create an empty repo named `ont-amr-nf` on GitHub (no README, no
.gitignore, no licence; this repo has its own), then:

```bash
git remote add origin git@github.com:YOUR_USERNAME/ont-amr-nf.git
git branch -M main
git push -u origin main
```

## 3. Repo settings worth two minutes

- **Description**: `Validated Nextflow pipeline for ONT long-read AMR calling, with a shuffled-assembly caller control and mutation-tested validation gates.`
- **Topics**: `nextflow` `nanopore` `amr` `bioinformatics` `long-read-sequencing` `genomics` `reproducibility`
- Leave Issues on. Turn Wikis and Projects off — an empty wiki tab reads as abandoned.

## 4. What a reviewer sees first

The README opens with the problem, then the design decisions, then results on real
data with every verdict traceable to `results/validation/`. The point of the repo is
that **the pipeline's own controls are checked by code, not asserted in prose** — so
if a reviewer reads one thing beyond the README, point them at
`tests/test_control_gate.py`.

## Not in this repo, deliberately

- **FASTQ files.** Accessions in `assets/samplesheet.csv` are the canonical reference
  and the pipeline re-fetches byte-identical input from SRA. Committing ~800 MB of
  reads to a portfolio repo is a mistake reviewers notice.
- **Basecalling.** Documented in the README as a step this pipeline does not perform;
  the public data is already basecalled. Claiming otherwise would be the one thing in
  here that isn't reproducible.
