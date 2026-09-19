# Run summary

- samples evaluated: **6**
- passed validation: **6/6**
- resistance determinants called: **93**
- other elements called (stress, virulence): **95**

## Per-sample

| Sample | Role | Verdict | Depth | N50 | AMR genes |
|---|---|---|---|---|---|
| `CALLER_CONTROL` | caller_control | **PASS** | n/a | n/a | 0 |
| `KP_ES_7983` | clonal_replicate | **PASS** | 26.1x | 5,310,109 | 20 |
| `EC_PE_M09449` | cross_species | **PASS** | 31.2x | 4,740,491 | 26 |
| `NEG_DECOY` | negative_control | **PASS** | 0.0x | 0 | 0 |
| `KP_BG_81` | test | **PASS** | 36.5x | 5,437,787 | 23 |
| `KP_ES_7636` | test | **PASS** | 29.9x | 5,425,447 | 20 |

## Control outcomes

**Read-level decoy** (tests the assembler) — `NEG_DECOY` returned **0** elements of any type (as required); flye status `no_assembly`, which is the expected outcome for shuffled reads.

**Caller-level control** (tests the gene caller) — a real assembly with its bases shuffled within each contig, identical in contig count, length and GC, returned **0** calls. The caller is matching on gene identity, not composition.
