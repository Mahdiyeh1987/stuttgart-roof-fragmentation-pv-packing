# Reproducibility guide

## Reproducibility levels

This repository distinguishes three types of material so that archived results are not confused with newly recomputed outputs.

### 1. Directly executable from the historical raw inputs

The empirical Stuttgart geometry/packing chain can be rerun from `data_raw/` using the repository-relative scripts. The recommended command is:

```bash
python scripts/reproduce_empirical.py
```

This command:

1. verifies analytical input hashes;
2. parses the supplied CityGML roof geometry;
3. repairs/excludes negligible geometry;
4. removes nearly nested coplanar surfaces;
5. consolidates adjacent coplanar CityGML partitions into physical roof facets;
6. transforms each physical facet to local roof-plane coordinates;
7. applies the external roof-envelope setback;
8. packs complete module rectangles using the production search;
9. regenerates the base building/facet tables into `reproduced/`;
10. regenerates the 14-facet denser-search benchmark; and
11. compares the regenerated aggregate shortfall with the canonical archived result.

Canonical publication tables under `results/` are never overwritten by this command.

### 2. Deterministically regenerated from archived result tables

`python scripts/make_figures.py` regenerates analytical Figures 1 and 3–7 locally from canonical repository tables. The generated files are written to `figures/`, which is intentionally excluded from version control. Figure 2 is not redistributed because it includes municipal contextual aerial imagery and is not part of the quantitative computation.

`python scripts/verify_reported_results.py` checks the manuscript crosswalk, including the 31/340 object counts, 14.26% aggregate shortfall, bootstrap interval, primary Spearman coefficient, FDR-adjusted p-value, controlled-experiment summary values, benchmark agreement, tilt subsets, ALKIS–CityGML screens, and leave-one-building-out range.

### 3. Archived controlled-experiment publication outputs

The controlled subdivision experiment is fully specified in the manuscript and supplementary material (outer geometry, setback, module, gap, seed, number of realizations, split selection rule, rejected-piece threshold, and descriptors). The recovered analysis bundle, however, contains the final controlled summary/correlation tables rather than the exact historical realization-level generator source.

Accordingly, this repository does **not** invent a replacement random generator and label it as the historical one. It archives and verifies the exact publication tables:

- `results/controlled_production_summary.csv`
- `results/controlled_production_correlations.csv`

and regenerates Figure 3 from them. This preserves the reported result faithfully while making the boundary of code-level regeneration explicit.

## Recommended commands

Fast repository audit:

```bash
python scripts/check_inputs.py
python scripts/verify_reported_results.py
python -m compileall scripts
```

Empirical base-case reproduction:

```bash
python scripts/run_scenario.py base
```

Packing-search benchmark:

```bash
python scripts/run_benchmark.py
```

Full convenience run:

```bash
python scripts/reproduce_empirical.py
```

Figure regeneration:

```bash
python scripts/make_figures.py  # optional local generation; output is git-ignored
```

## Random seeds

The submitted study uses seed `20260913` for the controlled fragmentation experiment and building-level bootstrap calculations. The canonical bootstrap table is included in `results/bootstrap_percentile.csv`.

## Numerical search used in the submitted benchmark

The production packing search evaluates edge-derived candidate orientations using a 9 × 9 offset grid and 10 mm clearance. The verification benchmark adds ±2.5° perturbations around production orientations and uses an 11 × 11 offset grid. The deterministic benchmark contains 14 physical facets spanning approximately 1.96–2,896.88 m².

The two searches agree on 12 of 14 facets. The denser search finds one extra module in each of the remaining two cases, for a pooled count difference of approximately 0.08%. This is a numerical sensitivity check, not a proof of global rectangle-packing optimality.

## Platform notes

The public entry-point scripts use `if __name__ == "__main__"` guards so Python multiprocessing works on Windows as well as Unix-like systems. The implementation modules are intentionally non-executable as standalone scripts; use the documented entry points instead.
