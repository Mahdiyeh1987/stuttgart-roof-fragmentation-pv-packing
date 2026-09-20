"""Verify archived repository outputs against the submitted manuscript values.

This is intentionally a fast integrity/crosswalk check. It does not recompute the
full geometry workflow; use ``run_scenario.py`` and ``run_benchmark.py`` for that.
"""
from __future__ import annotations

from pathlib import Path
import math
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"

CHECKS = []


def check(name, actual, expected, atol=1e-9):
    ok = math.isclose(float(actual), float(expected), rel_tol=0.0, abs_tol=atol)
    CHECKS.append((name, actual, expected, ok))
    if not ok:
        raise AssertionError(f"{name}: actual={actual!r}, expected={expected!r}, atol={atol}")


def main() -> int:
    buildings = pd.read_csv(RES / "base_buildings.csv")
    facets = pd.read_csv(RES / "base_facets.csv")
    boot = pd.read_csv(RES / "bootstrap_percentile.csv").set_index("metric")
    corr = pd.read_csv(RES / "correlations_FDR.csv").set_index("predictor")
    control = pd.read_csv(RES / "controlled_production_summary.csv").set_index("facet_count")
    ccorr = pd.read_csv(RES / "controlled_production_correlations.csv").set_index("descriptor")
    bench = pd.read_csv(RES / "packing_benchmark.csv")
    align = pd.read_csv(RES / "alkis_citygml_alignment_bestsubset.csv")
    sens = pd.read_csv(RES / "sensitivity_summary.csv").set_index("scenario")
    leave = pd.read_csv(RES / "leave_one_out.csv")
    reg_full = pd.read_csv(RES / "regression_HC3_full.csv")
    partial = pd.read_csv(RES / "partial_spearman_size_adjusted.csv")
    coplanar = pd.read_csv(RES / "coplanar_tolerance_sensitivity.csv")

    check("buildings", len(buildings), 31, 0)
    check("physical_facets", len(facets), 340, 0)
    check("usable_area_m2", buildings.usable_area_m2.sum(), 63587.465980, 1e-5)
    check("continuous_kwp", buildings.continuous_kwp.sum(), 14327.749935, 1e-6)
    check("packed_kwp", buildings.packed_kwp.sum(), 12283.92, 1e-8)
    aggregate = 100 * (1 - buildings.packed_kwp.sum() / buildings.continuous_kwp.sum())
    check("aggregate_shortfall_pct", aggregate, 14.26483533482279, 1e-10)
    check("median_building_shortfall_pct", buildings.packing_shortfall_pct.median(), 14.496689, 1e-6)

    check("bootstrap_aggregate_low", boot.loc["aggregate", "percentile_2.5"], 12.665321, 1e-6)
    check("bootstrap_aggregate_high", boot.loc["aggregate", "percentile_97.5"], 16.595697, 1e-6)
    check("rho_facet_density", corr.loc["facet_density_per_100m2_gross3d", "rho"], 0.631452, 1e-6)
    check("fdr_facet_density", corr.loc["facet_density_per_100m2_gross3d", "p_fdr_bh"], 0.000975, 1e-6)

    for n, expected in {2: 0.125, 3: 0.1875, 6: 0.3229166666666667, 10: 0.4375, 12: 0.4791666666666667}.items():
        check(f"controlled_median_loss_f{n}", control.loc[n, "median_loss"], expected, 1e-12)
    check(
        "controlled_rho_internal_boundary",
        ccorr.loc["internal_boundary_density_m_per_m2", "rho"],
        0.968349,
        1e-6,
    )
    check("controlled_rho_facet_count", ccorr.loc["facet_count", "rho"], 0.938962, 1e-6)

    check("benchmark_n", len(bench), 14, 0)
    check("benchmark_agreement", (bench.gap_modules == 0).sum(), 12, 0)
    check("benchmark_max_gap", bench.gap_modules.max(), 1, 0)
    pooled = 100 * (bench.stress_count.sum() - bench.production_count.sum()) / bench.stress_count.sum()
    check("benchmark_pooled_difference_pct", pooled, 0.08431703204047218, 1e-10)
    check("benchmark_min_area_m2", bench.usable_area_m2.min(), 1.96008597, 1e-7)
    check("benchmark_max_area_m2", bench.usable_area_m2.max(), 2896.881405, 1e-6)

    check("alignment_iou_ge_080", (align.best_iou >= 0.80).sum(), 27, 0)
    check("alignment_iou_ge_090", (align.best_iou >= 0.90).sum(), 26, 0)
    check("alignment_median_iou", align.best_iou.median(), 0.998, 0.001)

    # Tilt subsets are derived from the same base-case clipping, exactly as in the manuscript.
    check("near_flat_facet_count", (facets.tilt_deg < 5).sum(), 196, 0)
    check("pitched_facet_count", (facets.tilt_deg >= 5).sum(), 144, 0)
    check("near_flat_usable_area_m2", facets.loc[facets.tilt_deg < 5, "usable_area_m2"].sum(), 55812.213727, 1e-6)
    check("near_flat_shortfall_pct", sens.loc["flat", "aggregate_shortfall_pct"], 14.22734555, 1e-6)
    check("pitched_shortfall_pct", sens.loc["pitched", "aggregate_shortfall_pct"], 14.53394400, 1e-6)

    check("loo_rho_min", leave.rho_facet_density.min(), 0.593325917686318, 1e-12)
    check("loo_rho_max", leave.rho_facet_density.max(), 0.6765294771968854, 1e-12)

    # Final robustness checks added in v1.1.0.
    row = reg_full[(reg_full.model == "A_density_plus_log_area") & (reg_full.term == "facet_density_z")].iloc[0]
    check("hc3_density_coef_pp", row.coefficient_pp, 4.219599, 1e-6)
    check("hc3_density_p", row.p_value, 0.002562095, 1e-9)
    row = reg_full[(reg_full.model == "B_count_plus_log_area") & (reg_full.term == "facet_count_z")].iloc[0]
    check("hc3_raw_count_coef_pp", row.coefficient_pp, 1.204537, 1e-6)
    check("hc3_raw_count_p", row.p_value, 0.3762769, 1e-7)
    check("partial_spearman_rho", partial.loc[0, "rho"], 0.6216572157004039, 1e-12)
    check("partial_spearman_p", partial.loc[0, "p_two_sided"], 0.00024552347053714504, 1e-15)
    check("coplanar_facets_min", coplanar.physical_facets.min(), 337, 0)
    check("coplanar_facets_max", coplanar.physical_facets.max(), 342, 0)
    check("coplanar_shortfall_min", coplanar.aggregate_packing_shortfall_pct.min(), 14.231285321638376, 1e-12)
    check("coplanar_shortfall_max", coplanar.aggregate_packing_shortfall_pct.max(), 14.277119189634547, 1e-12)
    check("coplanar_rho_min", coplanar.spearman_rho_facet_density.min(), 0.6314516129032258, 1e-12)
    check("coplanar_rho_max", coplanar.spearman_rho_facet_density.max(), 0.6568548387096774, 1e-12)

    for name, actual, expected, ok in CHECKS:
        print(f"[{'OK' if ok else 'FAIL'}] {name}: {actual} (expected {expected})")
    print(f"\nVerified {len(CHECKS)} manuscript crosswalk checks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
