"""Regenerate the 14-facet packing-search verification benchmark.

The benchmark uses the production candidate orientations plus +/-2.5 degree
perturbations and an 11 x 11 offset grid, matching the submitted manuscript.
Results are written to ``reproduced/packing_benchmark.csv`` by default.
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from shapely import wkt
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import analysis_core as v


def _worker(payload):
    uid, polygon_wkt, angle = payload
    polygon = wkt.loads(polygon_wkt)
    count = v.count_angle_gap(polygon, angle, v.ha.MODULES[v.BASE], 0.01, 11)
    return uid, count


def build_benchmark(workers: int) -> pd.DataFrame:
    _, retained, _ = v.initial_rows()
    retained, _, _ = v.remove_nested(retained)
    facets, _, _ = v.merge_coplanar_adjacent(retained)

    module = v.ha.MODULES[v.BASE]
    by_building = defaultdict(list)
    for row in facets:
        by_building[row["building_id"]].append(row)

    records = []
    geometry = {}
    for building_id, rows in by_building.items():
        envelope = unary_union([wkt.loads(r["plan_wkt"]) for r in rows]).buffer(0)
        envelope = wkt.loads(envelope.wkt)
        for row in rows:
            usable = v.local_for_facet(row, envelope, 0.30)
            if usable.is_empty or usable.area < module["L"] * module["W"]:
                continue
            compactness = 4 * np.pi * usable.area / usable.length**2 if usable.length else np.nan
            uid = row["physical_facet_uid"]
            records.append(
                {
                    "building_id": building_id,
                    "surface_uid": uid,
                    "area": usable.area,
                    "compactness": compactness,
                }
            )
            geometry[uid] = usable

    candidates = pd.DataFrame(records)
    candidates["area_stratum"] = pd.qcut(
        candidates.area.rank(method="first"), 4, labels=False
    )
    candidates["comp_stratum"] = pd.qcut(
        candidates.compactness.rank(method="first"), 3, labels=False
    )

    selected = []
    for area_class in range(4):
        for compactness_class in range(3):
            subset = candidates[
                (candidates.area_stratum == area_class)
                & (candidates.comp_stratum == compactness_class)
            ].copy()
            za = (subset.area - subset.area.median()) / (subset.area.std(ddof=0) + 1e-12)
            zc = (subset.compactness - subset.compactness.median()) / (
                subset.compactness.std(ddof=0) + 1e-12
            )
            selected.append(subset.loc[(za * za + zc * zc).idxmin()])

    for ext_idx in [candidates.area.idxmin(), candidates.area.idxmax()]:
        row = candidates.loc[ext_idx]
        if not any(x.surface_uid == row.surface_uid for x in selected):
            selected.append(row)

    canonical = pd.read_csv(ROOT / "results" / "base_facets.csv").set_index(
        "physical_facet_uid"
    )

    tasks = []
    meta = {}
    for row in selected:
        uid = row.surface_uid
        polygon = geometry[uid]
        base_angles = [float(x) for x in v.edge_angles(polygon)]
        angles = sorted(
            set(((angle + delta) % 180) for angle in base_angles for delta in (-2.5, 0, 2.5))
        )
        meta[uid] = row
        tasks.extend((uid, polygon.wkt, angle) for angle in angles)

    best = defaultdict(int)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_worker, task) for task in tasks]
        for future in as_completed(futures):
            uid, count = future.result()
            best[uid] = max(best[uid], count)

    rows = []
    for index, row in enumerate(selected, 1):
        uid = row.surface_uid
        production = int(canonical.loc[uid, "modules"])
        reference = max(production, best[uid])
        rows.append(
            {
                "selection_index": index,
                "building_id": row.building_id,
                "surface_uid": uid,
                "usable_area_m2": float(row.area),
                "compactness": float(row.compactness),
                "area_stratum": int(row.area_stratum),
                "compactness_stratum": int(row.comp_stratum),
                "production_count": production,
                "stress_count": reference,
                "gap_modules": reference - production,
                "gap_pct": 100 * (reference - production) / reference if reference else 0,
                "reference_search": "production edge angles plus ±2.5° perturbations; 11x11 offset grid",
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "reproduced" / "packing_benchmark.csv",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, min(5, os.cpu_count() or 1)),
    )
    args = parser.parse_args()

    output = build_benchmark(args.workers)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)

    pooled = 100 * (output.stress_count.sum() - output.production_count.sum()) / output.stress_count.sum()
    print(output.to_string(index=False))
    print(f"area_range_m2={output.usable_area_m2.min():.8f}..{output.usable_area_m2.max():.8f}")
    print(f"agreement={int((output.gap_modules == 0).sum())}/{len(output)}")
    print(f"max_gap_modules={int(output.gap_modules.max())}")
    print(f"pooled_difference_pct={pooled:.12f}")
    print(f"output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
