# Changelog

All notable changes to the public replication package are documented here.

## 1.1.0 - 2026-09-20

- Added the final manuscript robustness checks without adding any publication image files.
- Added HC3 regression outputs for both facet-density + log-area and raw-facet-count + log-area specifications.
- Added the residual-rank partial Spearman result controlling for ranked log gross 3D roof area (rho = 0.622; two-sided p = 0.00025).
- Added the seven-pair coplanar-consolidation tolerance grid, including physical-facet counts, aggregate packing shortfall, and facet-density Spearman rho.
- Updated the manuscript-to-repository result crosswalk, verification script, reproducibility documentation, and supplementary material to match the final submission files.
- Corrected the bootstrap terminology to “95% bootstrap percentile interval”.
- Kept publication figures excluded from the repository and Zenodo package.

## 1.0.0 - 2026-09-16

- Prepared the public GitHub/Zenodo replication release for the submitted *Electricity* manuscript.
- Preserved exact historical Stuttgart analytical inputs and SHA-256 integrity manifest.
- Included the final 31-building / 340-physical-facet empirical results.
- Included the submitted 14-facet numerical packing-search benchmark.
- Included canonical controlled-fragmentation summary and correlation tables.
- Removed stale `base_fast_*` and independently recalculated tilt-subset tables that did not correspond to the submitted same-clipping sensitivity definition.
- Recomputed near-flat and pitched subset summaries directly from the base physical-facet table to match the manuscript (14.23% and 14.53%).
- Replaced absolute-path figure code with repository-relative figure generation.
- Added Windows-safe public multiprocessing entry points.
- Added manuscript result verification, data provenance, licensing, citation, Zenodo metadata, and CI integrity checks.
