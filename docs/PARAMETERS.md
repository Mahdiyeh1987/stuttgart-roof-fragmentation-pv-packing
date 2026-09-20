# Analysis parameters

## Base empirical scenario

| Parameter | Base value |
|---|---|
| CRS | EPSG:25832 (ETRS89 / UTM zone 32N) |
| External roof-envelope setback | 0.30 m in plan view |
| Base module | LONGi LR5-54HTB-440M |
| Module dimensions | 1.722 × 1.134 m |
| Module DC nameplate power | 440 W |
| Inter-module clearance | 10 mm |
| Packing orientation candidates | principal/dominant polygon direction, longest-edge directions, and orthogonal alternatives |
| Production translation search | 9 × 9 offsets over one module pitch |
| Module/facet rule | complete rectangle must lie within one usable physical facet; modules cannot cross physical-facet boundaries |

## Geometry audit and physical-facet consolidation

| Parameter | Value |
|---|---:|
| Negligible plan-area threshold | ≤0.05 m² |
| Nested-coplanar overlap threshold | ≥99% of smaller polygon |
| Plane-normal tolerance | ≤1° |
| Vertical separation tolerance | ≤0.10 m |
| Minimum shared boundary for adjacent merge | >0.05 m |
| Nearby tolerance check | 0.5–2° and 0.05–0.20 m |

The submitted analytical lineage is 450 raw roof polygons → 431 retained source polygons after audit/nested-coplanar exclusion → 340 physical facets after adjacent-coplanar consolidation.

## Controlled fragmentation experiment

| Parameter | Value |
|---|---:|
| Outer roof | 12 × 10 m |
| External setback | 0.30 m |
| Usable unfragmented area | 107.16 m² |
| Unfragmented module count | 48 |
| Facet counts | 2–12 (plus 1-facet reference) |
| Realizations per fragmented facet count | 50 |
| Seed | 20260913 |
| Facet selected for splitting | probability proportional to facet area |
| Minimum accepted piece | 0.8 m² |
| Internal setback | none |
| Module / gap / production packing rule | same base module, 10 mm gap, and 9 × 9 production search |

## Statistical analysis

| Component | Specification |
|---|---|
| Primary association | Spearman rank correlation: physical-facet density per 100 m² gross 3D roof area vs. building packing shortfall |
| Exploratory family | 7 geometric descriptors |
| Multiple-testing correction | Benjamini–Hochberg FDR across the 7 descriptors |
| Bootstrap | 10,000 building-level resamples with replacement |
| Bootstrap seed | 20260913 |
| Size-adjusted check A | HC3-robust OLS with standardized facet density and standardized log gross 3D roof area |
| Size-adjusted check B | HC3-robust OLS with standardized raw physical-facet count and standardized log gross 3D roof area |
| Partial-rank check | Residual-rank partial Spearman between facet density and packing shortfall, controlling for ranked log gross 3D roof area |
| Influence check | leave-one-building-out primary Spearman coefficient |

## Sensitivity scenarios

- external setback: 0.00, 0.30, 0.50 m;
- module: 440 W base and 575 W larger module;
- clearance: 0, 10, 20, 40 mm (0 mm only a theoretical lower bound);
- roof tilt: near-flat <5° and pitched ≥5°, using the same base clipping;
- influence: remove largest one / largest three buildings;
- spatial consistency: ALKIS–CityGML best-match IoU ≥0.80 and ≥0.90;
- coplanar consolidation: normal/vertical tolerance pairs (0.5°,0.05 m), (0.5°,0.10 m), (1.0°,0.05 m), (1.0°,0.10 m), (1.0°,0.20 m), (2.0°,0.10 m), and (2.0°,0.20 m).

## Numerical packing-search check

The submitted deterministic check selects 12 area/compactness-cell medoids and adds the smallest and largest eligible facets if not already represented, giving 14 benchmark facets. The denser reference search adds ±2.5° orientation perturbations and uses an 11 × 11 offset grid.
