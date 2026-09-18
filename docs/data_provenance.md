# Data provenance

All sequencing data in this project is **public, real clinical isolate data** retrieved from the NCBI Sequence Read Archive on 2026-09-18. Nothing is simulated except the explicitly-labelled negative control.

## Why these isolates

Oxford Nanopore R10.4.1 whole-genome runs of carbapenemase- or ESBL-producing *Enterobacterales* from clinical surveillance studies. The selection is deliberately structured so the pipeline can be *validated*, not merely executed:

| Sample | Run | Role in the validation design |
|---|---|---|
| `KP_ES_7636` | SRR36388034 | Carbapenemase-producing K. pneumoniae ST5994, Spain outbreak; expect carbapenemase gene |
| `KP_ES_7983` | SRR36388030 | Second isolate of the same ST5994 outbreak clone; AMR profile should agree with KP_ES_7636 |
| `KP_BG_81` | SRR38062537 | CRE surveillance isolate, Bulgaria; independent lineage, expect carbapenem resistance determinants |
| `EC_PE_M09449` | SRR38125693 | Carbapenemase co-producing E. coli, Peru; tests organism-aware calling in a second species |
| `NEG_DECOY` | synthetic | MUST yield zero AMR calls; guards against reference-bias / spurious hits |

## Per-isolate metadata

| Sample | Run | BioProject | BioSample | Species | Strain | Source | Country | Collected | Reads | Mean len | Est. depth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `KP_ES_7636` | [SRR36388034](https://www.ncbi.nlm.nih.gov/sra/SRR36388034) | PRJNA1378463 | SAMN53838466 | *Klebsiella pneumoniae* | #7636 | urine | Spain: Cantabria | 2023-12-13 | 56,169 | 3,331 bp | 35.3x |
| `KP_ES_7983` | [SRR36388030](https://www.ncbi.nlm.nih.gov/sra/SRR36388030) | PRJNA1378463 | SAMN53838470 | *Klebsiella pneumoniae* | #7983 | urine | Spain: Cantabria | 2023-05-06 | 87,529 | 2,038 bp | 33.7x |
| `KP_BG_81` | [SRR38062537](https://www.ncbi.nlm.nih.gov/sra/SRR38062537) | PRJNA1452497 | SAMN57234492 | *Klebsiella pneumoniae* | 81 | urine | Bulgaria: Sofia | 2025-01-14 | 16,975 | 12,811 bp | 41.0x |
| `EC_PE_M09449` | [SRR38125693](https://www.ncbi.nlm.nih.gov/sra/SRR38125693) | PRJNA1443524 | SAMN56767967 | *Escherichia coli* | M09-449 | urine | Peru:Lima | 2024 | 38,768 | 4,449 bp | 34.5x |

## Source studies

- **PRJNA1378463** — Emergence of CP-K. Pneumoniae ST5994 in Cantabria, Spain
- **PRJNA1452497** — CRE surveillance in University Hospital Lozenetz - Sofia, Bulgaria 2022-2025
- **PRJNA1443524** — Genomic characterization of clinical Enterobacterales isolates co-producing carbapenemases in Peru

## Chemistry verification

R10.4.1 was **confirmed from submission metadata**, not inferred from the search query. The SRA XML for these runs records `SQK-RBK114` / `SQK-RBK114-96` ligation-free barcoding kits and, where the submitter recorded it, Dorado SUP basecalling. Kit/chemistry strings were checked per run before selection.

## Expected resistance phenotype (the positive-control logic)

Every source study describes its isolates as carbapenemase-producing, carbapenem-resistant, or ESBL-producing. That is an independent, *a priori* expectation from the study design — so a pipeline run that fails to recover beta-lactamase determinants in these genomes is wrong, regardless of whether it completes without error. This is what makes them usable as positive controls.

## Negative controls

There are two, because they test different things.

`NEG_DECOY` (read-level) is generated inside the pipeline by shuffling real read sequences while preserving read length and base composition. It contains no true biological ORFs, so any AMR call against it is a false positive. It flows through the identical filter → assemble → call → gate path as the real isolates.

`CALLER_CONTROL` (caller-level) exists because the read-level decoy does not assemble — which is the correct outcome for shuffled reads, but it means the gene caller receives an empty FASTA and returns an empty call set trivially. That tests the assembler, not the caller. So a real, good assembly has its bases shuffled *within each contig*: contig count, contig lengths and GC content are preserved exactly, gene content is destroyed. It then goes through the same AMRFinderPlus process with the same parameters. A call there is driven by composition rather than gene identity. The measured composition of source vs. shuffled sequence is published to `results/controls/caller_control_composition.tsv` so the claim is checkable rather than asserted.

## Reproducing the retrieval

```bash
nextflow run . -profile conda --samplesheet assets/samplesheet.csv
```

Accessions are pinned in `assets/samplesheet.csv`; the pipeline fetches FASTQ from the NCBI Traces endpoint and records an md5 per run in `results/reads/<sample_id>.md5`, so a later run can be shown to have used byte-identical input.

The `Est. depth` column above is an *estimate from input bases over expected genome size* — it is the figure used to choose the isolates, not a result. The pipeline never uses it: realised depth is measured by remapping reads onto their own assembly and is reported in `results/assembly_metrics.tsv`.
