# Manuscript-to-repository results crosswalk

This table maps the principal submitted results to their canonical repository files.

| Manuscript quantity | Value | Repository source |
|---|---:|---|
| Buildings | 31 | `results/base_buildings.csv` |
| Physical roof facets | 340 | `results/base_facets.csv` |
| Setback-adjusted area | 63,587.5 m² | `results/base_buildings.csv` |
| Continuous-area capacity | 14.328 MWp | `results/base_buildings.csv` |
| Explicit packed capacity | 12.284 MWp | `results/base_buildings.csv` |
| Aggregate packing shortfall | 14.26% | `results/base_buildings.csv`, `results/sensitivity_summary.csv` |
| Aggregate bootstrap interval | 12.67–16.60% | `results/bootstrap_percentile.csv` |
| Median building shortfall | 14.50% | `results/base_buildings.csv` |
| Median bootstrap interval | 13.16–17.68% | `results/bootstrap_percentile.csv` |
| Facet-density Spearman ρ | 0.631 | `results/correlations_FDR.csv` |
| Facet-density FDR p | 0.0010 | `results/correlations_FDR.csv` |
| Facet-density bootstrap interval | 0.337–0.825 | `results/bootstrap_percentile.csv` |
| HC3 facet-density coefficient | 4.22 percentage points | `results/regression_HC3.csv` |
| HC3 95% CI | 1.48–6.96 | `results/regression_HC3.csv` |
| HC3 p | 0.0026 | `results/regression_HC3.csv` |
| HC3 density + log-area model R² / adjusted R² | 0.641 / 0.615 | `results/regression_HC3_full.csv` |
| Raw facet-count HC3 coefficient | 1.205 percentage points | `results/regression_HC3_full.csv` |
| Raw facet-count HC3 p | 0.3763 | `results/regression_HC3_full.csv` |
| Raw-count + log-area model R² / adjusted R² | 0.140 / 0.079 | `results/regression_HC3_full.csv` |
| Partial Spearman rho (size-adjusted) | 0.622 | `results/partial_spearman_size_adjusted.csv` |
| Partial Spearman two-sided p | 0.00025 | `results/partial_spearman_size_adjusted.csv` |
| Coplanar tolerance physical facets | 337–342 | `results/coplanar_tolerance_sensitivity.csv` |
| Coplanar tolerance aggregate shortfall | 14.23–14.28% | `results/coplanar_tolerance_sensitivity.csv` |
| Coplanar tolerance facet-density rho | 0.631–0.657 | `results/coplanar_tolerance_sensitivity.csv` |
| 14-facet benchmark agreement | 12/14 | `results/packing_benchmark.csv` |
| Benchmark pooled count difference | 0.08% | `results/packing_benchmark.csv` |
| No-setback shortfall | 14.40% | `results/sensitivity_summary.csv` |
| 0.50 m-setback shortfall | 14.44% | `results/sensitivity_summary.csv` |
| 575 W module shortfall | 15.46% | `results/sensitivity_summary.csv` |
| 0/10/20/40 mm clearance | 13.01/14.26/15.28/17.86% | `results/sensitivity_summary.csv` |
| Near-flat (<5°) shortfall | 14.23% | derived from `results/base_facets.csv`; recorded in `results/sensitivity_summary.csv` |
| Pitched (≥5°) shortfall | 14.53% | derived from `results/base_facets.csv`; recorded in `results/sensitivity_summary.csv` |
| Remove largest / three largest | 14.95 / 16.29% | `results/leave_largest_out.csv` |
| IoU ≥0.80 / ≥0.90 | 13.69 / 13.68% | `results/alignment_sensitivity.csv` |
| Leave-one-building-out ρ range | 0.593–0.677 | `results/leave_one_out.csv` |
| Controlled loss: 2 facets | 12.5% | `results/controlled_production_summary.csv` |
| Controlled loss: 3 facets | 18.8% | `results/controlled_production_summary.csv` |
| Controlled loss: 6 facets | 32.3% | `results/controlled_production_summary.csv` |
| Controlled loss: 10 facets | 43.8% | `results/controlled_production_summary.csv` |
| Controlled loss: 12 facets | 47.9% | `results/controlled_production_summary.csv` |
| Controlled internal-boundary ρ | 0.968 | `results/controlled_production_correlations.csv` |

Run `python scripts/verify_reported_results.py` to check these crosswalk values programmatically.
