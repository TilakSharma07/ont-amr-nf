# ont-amr-nf

**Oxford Nanopore (R10.4.1) bacterial reads → validated antimicrobial-resistance gene calls.**

A Nextflow DSL2 pipeline that takes long-read sequencing data from clinical bacterial
isolates and produces resistance-gene calls **together with a verdict on whether those
calls are interpretable**. Run on real carbapenemase- and ESBL-producing clinical
isolates from public NCBI SRA data, with a synthetic negative control and a depth
titration that establishes where the method stops working.

---

## The problem this addresses

A resistance-gene caller will happily return a clean-looking table from a bad assembly.
Two failure modes matter clinically and neither raises an error:

1. **False negative from insufficient data.** Absence of `blaKPC` in a 6× assembly is not
   absence of `blaKPC`. Reported as "no carbapenemase detected", it is a wrong answer
   that looks like a right one.
2. **False positive from composition.** A caller matching on sequence composition rather
   than gene identity will report determinants on input that contains no genes at all.

So this pipeline does not just call genes. Every sample passes a **validation gate** that
can mark its result *uninterpretable*, and the run includes controls designed to make both
failure modes visible if they occur.

---

## Design decisions worth defending

**Depth is measured, not assumed.** Reads are remapped onto their own assembly with
minimap2 and depth comes from `samtools depth`. Estimating coverage as
`input_bases / expected_genome_size` is wrong twice over: reads removed by filtering are
still counted, and content that failed to assemble is counted as if it assembled.

**The negative control shares one code path with the real isolates.** `MAKE_DECOY`
shuffles the bases within each read of a real sample, preserving the read-length
distribution and per-read base composition and destroying only sequence order. It then
flows through the *identical* filter → assemble → call → gate path. A control processed
by different code proves nothing about the pipeline that processed the samples.

**There are two negative controls, because one is not enough.** The read-level decoy
does not assemble — correctly, that is what shuffled reads should do. But that means the
gene caller only ever receives an empty FASTA from it and returns an empty call set
trivially. That control tests the **assembler**; it says nothing about the **caller**.
So `SHUFFLE_ASSEMBLY` takes a real, good assembly and shuffles bases within each contig:
contig count, contig lengths and GC content are preserved exactly, gene content is
destroyed. The result goes through the same `AMRFINDERPLUS` process with the same
parameters. A call there is a call driven by composition rather than gene identity — a
false positive by construction — and the run summary says so in those terms.

**The decoy is run without `--organism`, deliberately.** No organism claim can be made
about shuffled sequence, and the goal is the most permissive call mode possible on it. If
anything survives there, it should be seen rather than filtered away.

**A failed assembly still emits a file.** `FLYE_ASSEMBLE` is not declared `optional`.
If it were, the decoy — which is *expected* to assemble poorly or not at all — would drop
out of the channel and its "negative control is empty" check would never execute. The
pipeline would then look clean precisely by omitting its own control. Instead a
zero-sequence FASTA plus a status file carries the failure forward to the gate.

**"No calls" is disambiguated.** The gate reads the assembly status alongside the call
set, because *no calls because nothing assembled* and *no calls in a good assembly* are
different findings that must not collapse into the same row.

**Role-specific expectations.** The isolates are known clinical carbapenemase/ESBL
producers, so zero calls on them indicates pipeline failure, not a susceptible isolate —
the gate asserts `positive_expectation_met`. For the decoy only `negative_control_is_empty`
is binding: a decoy that assembles badly has not failed, that is the point of it.

**Titration subsamples to a target depth, not a read count.** Read-length distributions
differ between isolates, so a fixed read count yields different coverage per sample and
the resulting curve would not be comparable.

---

## Quick start

```bash
# stub run — validates the whole DAG in seconds, no data or tools needed
nextflow run . -stub-run

# real run, tools resolved per process by conda
nextflow run . -profile conda

# real run against reads already on disk (stays offline)
nextflow run . -profile conda --reads_dir /path/to/fastq
```

Reads are fetched from SRA by accession if not found locally; every fetch is md5-recorded
in `results/reads/*.md5`.

## Outputs

| File | Contents |
|---|---|
| `results/amr_calls.tsv` | every resistance determinant, per sample, with identity/coverage |
| `results/validation_summary.tsv` | PASS/FAIL per sample and per check, with the failing detail |
| `results/assembly_metrics.tsv` | contigs, N50, total length, GC, measured depth, breadth |
| `results/depth_titration.tsv` | determinants recovered vs. depth (when `--run_titration`) |
| `results/controls/` | control composition reports (source vs. shuffled) |
| `results/run_summary.md` | human-readable overview, including both control outcomes |
| `results/pipeline_info/` | execution report, trace, DAG, resolved tool versions |

## Key parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--min_read_len` | 1000 | read length floor (bp) |
| `--min_read_q` | 10 | read quality floor |
| `--min_depth_x` | 20 | depth below which AMR *absence* is not reported |
| `--min_n50` | 50000 | assembly N50 floor for an interpretable result |
| `--amr_min_ident` | 0.9 | AMRFinderPlus identity floor |
| `--amr_min_cov` | 0.5 | AMRFinderPlus reference-coverage floor |
| `--make_decoy` | true | read-level control: shuffled reads through the full path |
| `--caller_control` | true | caller-level control: shuffled assembly, re-called |
| `--run_titration` | false | run the depth-titration experiment |

Every threshold that can change a result is a parameter — none are buried in a script.
`nextflow run . --help` lists them all.

## Pipeline steps

```
samplesheet
   ↓
FETCH_READS ──────── md5-recorded, local copy reused if present
   ↓
   ├─→ MAKE_DECOY ──────────── read-level control  (tests the assembler)
   ├─→ SUBSAMPLE_DEPTH ─────── depth titration points (optional)
   ↓
CALL_AMR  (one identical path for isolates, decoy and titration points)
   │
   ├─ NANOQ_FILTER ─────── length/quality filter
   ├─ FLYE_ASSEMBLE ────── --nano-hq (correct mode for R10.4.1 + SUP)
   ├─ ASSEMBLY_STATS ───── contiguity + depth measured by remapping
   ├─ AMRFINDERPLUS ────── organism-aware for isolates, permissive for controls
   └─ VALIDATION_GATE ──── PASS/FAIL + amr_result_interpretable
   │
   └─→ SHUFFLE_ASSEMBLY ──→ AMRFINDERPLUS   caller-level control
                                            (same process, same parameters)
   ↓
AGGREGATE_RESULTS ──── tables + run summary
```

## Data

Four real R10.4.1 clinical isolates, selected so that the sample set tests something:
a same-outbreak pair (call concordance), an independent lineage, and a second species
(organism-aware calling). Chemistry was confirmed from each submission's own kit and
basecaller fields rather than from a keyword match. Accessions, BioProjects, originating
studies and expected phenotypes are in
[`docs/data_provenance.md`](docs/data_provenance.md).

## Scope and limits

- **Basecalling is not performed here.** The public input is already basecalled with
  dorado SUP. Raw-signal basecalling needs a GPU and the POD5 signal data; the pipeline
  documents the chemistry and basecaller of its inputs rather than re-deriving them.
- **Resistance genotype is not resistance phenotype.** A detected determinant is not an
  MIC. Expected phenotypes in the provenance doc come from the originating studies and
  are used here as positive-control priors, not as validated susceptibility results.
- **Not a diagnostic device.** This is a reproducible research pipeline. Clinical use
  requires validation, accreditation and controls far beyond its scope.

## Citations

Tool and reference-database citations: [`CITATIONS.md`](CITATIONS.md).

## License

MIT — see [`LICENSE`](LICENSE).
