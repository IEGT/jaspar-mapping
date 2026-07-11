# Five-Base Dense-Score Example

This fixture is deliberately small enough to inspect without software:

```text
reference:  A C G T N
start:      0 1 2 3 4
```

The one-column `A_RICH` motif has counts `A=3` and `C=G=T=1`. Its column sum
is 6 and the test uses pseudocount 0.

For `log2_relative_risk`:

```text
A:       log2((3/6) / 0.25) =  1
C/G/T:   log2((1/6) / 0.25) = -0.5849625007
```

For `log_odds`:

```text
A:       log2((3/(6-3)) / (0.25/0.75)) =  1.5849625007
C/G/T:   log2((1/(6-1)) / (0.25/0.75)) = -0.7369655942
```

The plus orientation therefore scores `A` highly at start 0. The minus
orientation scores the reverse complement, so reference `T` scores highly at
start 3. The `N` at start 4 is skipped and represented as `NULL`, not zero.
All expected values are visible in [`expected_scores.tsv`](expected_scores.tsv).

Run the ephemeral test:

```sh
make check_synthetic_dense
```

Or retain a tiny package for interactive inspection:

```sh
make synthetic_dense_example
cd dry_runs/synthetic_dense
duckdb -readonly synthetic_dense.duckdb \
  "SELECT * FROM synthetic_comparison ORDER BY start;"
```

The comparison table places expected and actual float32 values side by side
and ends each row with `matches=true`.
