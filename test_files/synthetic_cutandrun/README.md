# Synthetic TP73 CUT&RUN Containment Example

This fixture is small enough to inspect directly. All coordinates are BED
0-based half-open. Each synthetic motif is four bases long:

| Motif | Interval | Score | Expected strict depth |
| --- | --- | ---: | ---: |
| TP73_A | `[10,14)` | 10 | 1 |
| TP73_B | `[20,24)` | 8 | 1 |
| TP73_C | `[30,34)` | 6 | 2 |
| TP73_D | `[40,44)` | 4 | 0 |
| TP73_E | `[50,54)` | 2 | 0 |
| TP73_F | `[60,64)` | 0 | 1 |
| TP73_G | `[70,74)` | -2 | 1 |
| TP73_H | `[80,84)` | -4 | 0 |

A merged positive-coverage component supports a motif only under the strict
rule:

```text
component_start < motif_start AND component_end > motif_end
```

Overlapping and directly adjacent input intervals are merged before applying
the rule. `equal_left_D=[40,46)` covers all of `TP73_D=[40,44)`, but starts
exactly at the motif boundary. `equal_right_E=[48,54)` also covers all of
`TP73_E=[50,54)`, but ends exactly at its boundary. Their ordinary coverage
depth is therefore 1 while their effective depth is 0. `join_G_left=[68,72)`
and `join_G_right=[72,76)` are adjacent and jointly form `[68,76)`, which
immerses TP73_G although neither source interval does so alone. The two
overlapping intervals covering TP73_C produce effective maximum depth 2.

`tp73_coverage.bedGraph` encodes the same support pattern as positive
run-length coverage. Its fourth column is depth rather than a fragment name.
The expected ordinary/effective maxima for A through H are respectively
`3/3`, `2/2`, `4/4`, `5/0`, `6/0`, `7/7`, `9.5/9.5`, and `0/0`. The fractional
value deliberately verifies that bedGraph signal is neither counted as one row
nor truncated to an integer.

The five supported motifs have scores 10, 8, 6, 0, and -2. The three
unsupported motifs have scores 4, 2, and -4. This gives a non-perfect ROC AUC
of 0.7333 and average precision of 0.8762. Maximum Youden J selects threshold
6, while maximum F1 selects -2.

Run the retained local example with:

```sh
make synthetic_cutandrun_example
```

The generated audit tables are written under
`dry_runs/synthetic_cutandrun_coverage_union/`.
The test uses an ephemeral directory:

```sh
make check_cutandrun_containment
```
