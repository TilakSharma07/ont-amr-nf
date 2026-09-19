# Run summary

- samples evaluated: **5**
- passed validation: **2/5**
- total AMR determinant calls: **205**

## Per-sample

| Sample | Role | Verdict | Depth | N50 | AMR genes |
|---|---|---|---|---|---|
| `KP_ES_7636` | test | **PASS** | 29.9x | 5,425,444 | 20 |
| `KP_ES_7636_d10_r1` | titration | **FAIL** | 8.4x | 203,873 | 18 |
| `KP_ES_7636_d20_r1` | titration | **FAIL** | 16.9x | 5,425,465 | 20 |
| `KP_ES_7636_d40_r1` | titration | **PASS** | 29.9x | 5,425,446 | 20 |
| `KP_ES_7636_d5_r1` | titration | **FAIL** | 4.7x | 53,206 | 17 |

## Depth titration

Recovery fraction is measured against the full-depth call set of the same parent isolate.

| Point | Target | Realised | Genes | Recovery |
|---|---|---|---|---|
| `KP_ES_7636_d10_r1` | 10x | 8.4x | 18/20 | 75% |
| `KP_ES_7636_d20_r1` | 20x | 16.9x | 20/20 | 100% |
| `KP_ES_7636_d40_r1` | 40x | 29.9x | 20/20 | 100% |
| `KP_ES_7636_d5_r1` | 5x | 4.7x | 17/20 | 70% |

## Control outcomes

No read-level decoy was run (`--make_decoy false`).

No caller-level control was run (`--caller_control false`).
