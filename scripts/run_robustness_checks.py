"""Reproduce the robustness checks added for the final submitted manuscript.

Outputs (written to ``reproduced/`` by default):
- regression_HC3_full.csv
- partial_spearman_size_adjusted.csv
- coplanar_tolerance_sensitivity.csv

The canonical archived copies are stored under ``results/``.  The regression and
partial-rank checks are computed from ``results/base_buildings.csv``.  The
coplanar-tolerance grid is recomputed from the archived historical CityGML input;
for efficiency, only buildings whose physical-facet grouping changes relative to
the base 1.0 degree / 0.10 m rule are repacked.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, t as student_t
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))
import analysis_core as ac  # noqa: E402

BASE_BUILDINGS = ROOT / "results" / "base_buildings.csv"

TOLERANCE_GRID = [
    (0.5, 0.05),
    (0.5, 0.10),
    (1.0, 0.05),
    (1.0, 0.10),
    (1.0, 0.20),
    (2.0, 0.10),
    (2.0, 0.20),
]


def _z(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=0)


def regression_checks(base: pd.DataFrame) -> pd.DataFrame:
    d = base[[
        "packing_shortfall_pct",
        "facet_density_per_100m2_gross3d",
        "n_facets",
        "gross_roof_area_3d_m2",
    ]].dropna().copy()
    d["log_area"] = np.log(d["gross_roof_area_3d_m2"])
    d["facet_density_z"] = _z(d["facet_density_per_100m2_gross3d"])
    d["facet_count_z"] = _z(d["n_facets"])
    d["log_area_z"] = _z(d["log_area"])

    specs = [
        ("A_density_plus_log_area", ["facet_density_z", "log_area_z"]),
        ("B_count_plus_log_area", ["facet_count_z", "log_area_z"]),
    ]
    rows = []
    for label, predictors in specs:
        model = sm.OLS(
            d["packing_shortfall_pct"], sm.add_constant(d[predictors])
        ).fit(cov_type="HC3")
        ci = model.conf_int()
        for term in model.params.index:
            rows.append({
                "model": label,
                "term": term,
                "coefficient_pp": float(model.params[term]),
                "hc3_se": float(model.bse[term]),
                "ci_95_low": float(ci.loc[term, 0]),
                "ci_95_high": float(ci.loc[term, 1]),
                "p_value": float(model.pvalues[term]),
                "r_squared": float(model.rsquared),
                "adjusted_r_squared": float(model.rsquared_adj),
                "n": int(model.nobs),
            })
    return pd.DataFrame(rows)


def partial_spearman_check(base: pd.DataFrame) -> pd.DataFrame:
    d = base[[
        "packing_shortfall_pct",
        "facet_density_per_100m2_gross3d",
        "gross_roof_area_3d_m2",
    ]].dropna().copy()
    log_area = np.log(d["gross_roof_area_3d_m2"].to_numpy(float))
    rx = rankdata(d["facet_density_per_100m2_gross3d"].to_numpy(float))
    ry = rankdata(d["packing_shortfall_pct"].to_numpy(float))
    rz = rankdata(log_area)

    design = sm.add_constant(rz)
    ex = sm.OLS(rx, design).fit().resid
    ey = sm.OLS(ry, design).fit().resid
    rho = float(np.corrcoef(ex, ey)[0, 1])

    # One control variable -> df = n - k - 2 = n - 3.
    n = len(d)
    controls = 1
    df = n - controls - 2
    t_stat = rho * math.sqrt(df / (1.0 - rho * rho))
    p_two_sided = float(2.0 * student_t.sf(abs(t_stat), df))

    return pd.DataFrame([{
        "method": "partial_spearman_residualized_ranks",
        "x": "facet_density_per_100m2_gross3d",
        "y": "packing_shortfall_pct",
        "control": "log_gross_roof_area_3d_m2",
        "rho": rho,
        "t_statistic": float(t_stat),
        "df": int(df),
        "p_two_sided": p_two_sided,
        "n": int(n),
    }])


def _group_signature(facets: list[dict]) -> dict[str, tuple[str, ...]]:
    by_building: dict[str, list[str]] = {}
    for r in facets:
        members = tuple(sorted(str(r["source_surface_uids"]).split("|")))
        token = "|".join(members)
        by_building.setdefault(str(r["building_id"]), []).append(token)
    return {bid: tuple(sorted(tokens)) for bid, tokens in by_building.items()}


def coplanar_tolerance_checks(base: pd.DataFrame, workers: int) -> pd.DataFrame:
    _raw, cleaned, _excluded = ac.initial_rows()
    retained, _removed, _audit = ac.remove_nested(
        cleaned, angle_tol=1.0, z_tol=0.10, contain=0.99
    )
    base_facets, _pairs, _membership = ac.merge_coplanar_adjacent(
        retained, angle_tol=1.0, z_tol=0.10
    )
    base_sig = _group_signature(base_facets)
    base_indexed = base.set_index("building_id", drop=False)

    rows = []
    for angle_tol, z_tol in TOLERANCE_GRID:
        facets, _pairs, _membership = ac.merge_coplanar_adjacent(
            retained, angle_tol=angle_tol, z_tol=z_tol
        )
        sig = _group_signature(facets)
        affected = sorted(
            bid for bid in set(base_sig) | set(sig)
            if base_sig.get(bid) != sig.get(bid)
        )

        combined = base.copy()
        if affected:
            variant_facets = [r for r in facets if str(r["building_id"]) in affected]
            recalculated, _facet_rows = ac.scenario(
                variant_facets, setback=0.30, mkey=ac.BASE, gap=0.01,
                workers=max(1, workers)
            )
            combined = combined[~combined["building_id"].isin(affected)].copy()
            combined = pd.concat([combined, recalculated], ignore_index=True)

        packed = float(combined["packed_kwp"].sum())
        continuous = float(combined["continuous_kwp"].sum())
        shortfall = 100.0 * (1.0 - packed / continuous)
        rho, p_raw = spearmanr(
            combined["facet_density_per_100m2_gross3d"],
            combined["packing_shortfall_pct"],
        )
        rows.append({
            "normal_tolerance_deg": angle_tol,
            "vertical_tolerance_m": z_tol,
            "physical_facets": len(facets),
            "affected_buildings": len(affected),
            "aggregate_packing_shortfall_pct": shortfall,
            "spearman_rho_facet_density": float(rho),
            "p_raw": float(p_raw),
        })
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "reproduced",
        help="Directory for regenerated outputs (default: reproduced/).",
    )
    parser.add_argument(
        "--skip-coplanar", action="store_true",
        help="Run only the fast regression and partial-Spearman checks.",
    )
    parser.add_argument(
        "--workers", type=int, default=4,
        help="Worker processes for repacking affected buildings in the tolerance grid.",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    base = pd.read_csv(BASE_BUILDINGS)
    reg = regression_checks(base)
    partial = partial_spearman_check(base)
    reg.to_csv(args.output_dir / "regression_HC3_full.csv", index=False)
    partial.to_csv(args.output_dir / "partial_spearman_size_adjusted.csv", index=False)

    print(reg.to_string(index=False))
    print("\n", partial.to_string(index=False), sep="")

    if not args.skip_coplanar:
        tol = coplanar_tolerance_checks(base, workers=args.workers)
        tol.to_csv(args.output_dir / "coplanar_tolerance_sensitivity.csv", index=False)
        print("\n", tol.to_string(index=False), sep="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
