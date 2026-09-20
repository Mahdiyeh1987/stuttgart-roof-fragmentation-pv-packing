# Results directory

This directory contains the **canonical processed outputs used by the submitted manuscript**. Public reproduction commands write to `reproduced/` and do not overwrite these files.

Important tables:

- `base_buildings.csv` — building-level base-case metrics and capacities.
- `base_facets.csv` — physical-facet-level base-case area, tilt, lineage count, and module count.
- `physical_facets.csv` — consolidated physical roof-facet geometry/lineage table.
- `surface_to_physical_facet.csv` — source CityGML surface → physical facet mapping.
- `geometry_basic_exclusions.csv`, `nested_overlap_audit.csv`, `coplanar_adjacency_audit.csv` — geometry audit trail.
- `correlations_FDR.csv` — seven pre-specified geometric-descriptor Spearman tests with BH-FDR correction.
- `regression_HC3.csv` — original size-adjusted HC3-robust density + log-area regression.
- `regression_HC3_full.csv` — final two-model HC3 table: density + log area and raw facet count + log area.
- `partial_spearman_size_adjusted.csv` — residual-rank partial Spearman check controlling for ranked log gross 3D roof area.
- `coplanar_tolerance_sensitivity.csv` — seven-pair normal/vertical coplanar-consolidation tolerance grid with facet count, aggregate shortfall, and facet-density Spearman rho.
- `bootstrap_percentile.csv` — 10,000-resample building-level percentile bootstrap results.
- `packing_benchmark.csv` — submitted deterministic 14-facet search check.
- `controlled_production_summary.csv`, `controlled_production_correlations.csv` — canonical controlled-experiment publication outputs.
- `sensitivity_summary.csv` — installation/setback/module/tilt sensitivity summary. The near-flat and pitched rows are computed from the **same base clipping**, matching the submitted manuscript.
- `alignment_sensitivity.csv`, `alkis_citygml_alignment_bestsubset.csv` — ALKIS–CityGML planimetric consistency results.
- `leave_largest_out.csv`, `leave_one_out.csv` — influence checks.
- `input_file_manifest_sha256.csv` — exact analytical input integrity manifest.

See `docs/RESULTS_CROSSWALK.md` for manuscript values and their source files.
