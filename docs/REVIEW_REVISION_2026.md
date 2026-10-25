# Electricity 2026 Major Revision – verification and sensitivity addendum

This repository update accompanies the revised manuscript **Roof-Facet Fragmentation and Rooftop PV Packing: A Controlled LoD2 Study with a Stuttgart Pilot** (Tabatabaei and Cremers). It supplements, but does **not** overwrite, the original v1.1.0 canonical `results/` tables or change historical analytical raw data.

## Additional manuscript-revision analyses

1. **Geometric orientation restriction** (not tilt/mounting design): `scripts/run_orientation_ablation.py`, with archived outputs in `revision_results/orientation_axis_ablation_*`. Production unrestricted search retains 14.264835% shortfall; the dominant-axis / orthogonal-axis / best-of-two restrictions produced approximately 18.35% / 15.79% / 15.27% shortfall.
2. **Two-factor grid**: `scripts/run_revision_joint_sensitivity.py` for all 3 setbacks (0, 0.30, 0.50 m) × 2 module clearances (10, 40 mm). Results are stored in `revision_results/setback_clearance_joint_sensitivity.csv`.
3. **Independent analytic oracle**: `python validation/independent_exact_benchmark.py`. Test geometries have exact module-count witnesses/area upper bounds. Under strict containment 12/17 passed; with `1e-8 m` numerical boundary tolerance 17/17 passed. Not a proof of optimality on arbitrary shapes.
4. **Full empirical precision audit**: `python validation/empirical_precision_audit.py`. Recomputes 340 physical facets using the numerical tolerance; no module count changed, and the archived 14.264835% shortfall was preserved.
5. **New controlled Monte Carlo replication**: `python validation/controlled_independent.py --repeats 50`. Generates 50 independent splits per facet count 2–12 (550 total). This is a **new generator** under the reported rules, not a recovery of the historical random realizations. The new 12-facet median is 47.92%; internal-boundary density Spearman rho is 0.955, versus 0.968 originally.
6. **Precision robustness of new controlled experiment**: `python validation/controlled_independent_precision.py`. The same 550 independently generated splits give identical module counts using the numerical tolerance.

### Important limitations

- The historical random-realization source and historical random polygons were not recovered; the original `results/controlled_production_*` publications are retained as historical archive.
- The `1e-8 m` boundary tolerance applies **only** in the validation scripts through a runtime override. It does not silently change the original production algorithm or physical setback/clearance.
- These tests concern geometric packing, not as-built installation, irradiance, shadowing, structure, electrical wiring, roof-access rules or global combinatorial optimality.
- Runtime observations depend on the computing environment; they are not a hardware-neutral performance benchmark.

### Execution

From the root of the existing v1.1.0 repository, with its original `scripts/`, `data_raw/`, and `results/` preserved:

```bash
python scripts/check_inputs.py
python scripts/verify_reported_results.py
python validation/independent_exact_benchmark.py
python validation/empirical_precision_audit.py
python validation/controlled_independent.py --repeats 50
python validation/controlled_independent_precision.py
python scripts/run_orientation_ablation.py
python scripts/run_revision_joint_sensitivity.py 0.30 0.01
```

For all combinations of the joint grid, repeat the last command with `(setback, gap)` of `(0, .01)`, `(0, .04)`, `(.30, .01)`, `(.30, .04)`, `(.50, .01)`, and `(.50, .04)`. The sensitivity scripts write to a local `revised_checks/` directory; the immutable tested outputs are archived under `revision_results/`.

### Release handling

Do not publish an invented new Zenodo DOI. After review, tag a new GitHub version and create the corresponding Zenodo version; use its actual DOI in the final manuscript. Keep the submitted article, public supplementary file, and release documentation consistent. This addendum is a candidate update, **not evidence that any GitHub/Zenodo release has already been published**.
