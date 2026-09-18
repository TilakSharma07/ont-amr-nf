# Citations

## Workflow manager

- **Nextflow** — Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C.
  *Nextflow enables reproducible computational workflows.* Nat Biotechnol 35, 316–319 (2017).
  doi:10.1038/nbt.3820

## Pipeline tools

- **nanoq** — Steinig E, Coin L. *Nanoq: ultra-fast quality control for nanopore reads.*
  J Open Source Softw 7(69), 2991 (2022). doi:10.21105/joss.02991
  Used for read length/quality filtering and read-set statistics.

- **Flye** — Kolmogorov M, Yuan J, Lin Y, Pevzner PA. *Assembly of long, error-prone reads
  using repeat graphs.* Nat Biotechnol 37, 540–546 (2019). doi:10.1038/s41587-019-0072-8
  Run in `--nano-hq` mode, which is the mode intended for R10.4.1 chemistry with
  high-accuracy basecalling.

- **minimap2** — Li H. *Minimap2: pairwise alignment for nucleotide sequences.*
  Bioinformatics 34(18), 3094–3100 (2018). doi:10.1093/bioinformatics/bty191
  Used to remap reads onto their own assembly so depth is measured rather than assumed.

- **SAMtools** — Danecek P, Bonfield JK, Liddle J, et al. *Twelve years of SAMtools and
  BCFtools.* GigaScience 10(2), giab008 (2021). doi:10.1093/gigascience/giab008

- **SeqKit** — Shen W, Le S, Li Y, Hu F. *SeqKit: a cross-platform and ultrafast toolkit
  for FASTA/Q file manipulation.* PLoS One 11(10), e0163962 (2016).
  doi:10.1371/journal.pone.0163962
  Used for assembly statistics and for seeded depth subsampling in the titration experiment.

- **AMRFinderPlus** — Feldgarden M, Brover V, Gonzalez-Escalona N, et al. *AMRFinderPlus
  and the Reference Gene Catalog facilitate examination of the genomic links among
  antimicrobial resistance, stress response, and virulence.* Sci Rep 11, 12728 (2021).
  doi:10.1038/s41598-021-91456-0
  Run with `--organism` for the test isolates so species-specific rules and intrinsic-gene
  suppression apply.

## Reference data

- **NCBI Reference Gene Catalog / AMRFinderPlus database** — the resistance-gene reference
  set distributed with AMRFinderPlus. The exact database version used in a run is recorded
  in `results/pipeline_info/software_versions.yml`.

- **NCBI Sequence Read Archive** — Leinonen R, Sugawara H, Shumway M.
  *The Sequence Read Archive.* Nucleic Acids Res 39(suppl_1), D19–D21 (2011).
  doi:10.1093/nar/gkq1019

## Source datasets

Per-isolate accessions, BioProjects and the studies that generated them are listed in
[`docs/data_provenance.md`](docs/data_provenance.md). Please cite the originating studies
if you reuse these isolates.
