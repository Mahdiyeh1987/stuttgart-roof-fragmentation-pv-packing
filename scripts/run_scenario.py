"""Run one empirical packing scenario from the historical Stuttgart inputs.

Outputs are written to ``reproduced/`` by default so the archived canonical
results in ``results/`` are never overwritten accidentally.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import analysis_core as v

SCENARIOS = {
    "base": (0.30, v.BASE, 0.01, None),
    "setback0": (0.00, v.BASE, 0.01, None),
    "setback05": (0.50, v.BASE, 0.01, None),
    "large575": (0.30, "LONGi_LR5_72HGD_575M", 0.01, None),
    "gap0": (0.30, v.BASE, 0.00, None),
    "gap02": (0.30, v.BASE, 0.02, None),
    "gap04": (0.30, v.BASE, 0.04, None),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenario", nargs="?", default="base", choices=sorted(SCENARIOS))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "reproduced",
        help="Directory for regenerated CSV files (default: %(default)s)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, min(5, os.cpu_count() or 1)),
        help="Parallel workers for facet packing (default: up to 5)",
    )
    args = parser.parse_args()

    raw, retained, exclusions = v.initial_rows()
    retained2, _, _ = v.remove_nested(retained)
    facets, _, _ = v.merge_coplanar_adjacent(retained2)

    setback, module, gap, tilt_filter = SCENARIOS[args.scenario]
    buildings, facet_table = v.scenario_fast(
        facets, setback, module, gap, tilt_filter, workers=args.workers
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    buildings.to_csv(args.output_dir / f"{args.scenario}_buildings.csv", index=False)
    facet_table.to_csv(args.output_dir / f"{args.scenario}_facets.csv", index=False)

    shortfall = 100.0 * (1.0 - buildings.packed_kwp.sum() / buildings.continuous_kwp.sum())
    print(f"scenario={args.scenario}")
    print(f"raw_roof_polygons={len(raw)}")
    print(f"basic_retained_polygons={len(retained)}")
    print(f"physical_facets={len(facets)}")
    print(f"buildings={len(buildings)}")
    print(f"aggregate_shortfall_pct={shortfall:.12f}")
    print(f"outputs={args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
