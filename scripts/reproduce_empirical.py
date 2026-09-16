"""Convenience runner for the empirical base case and numerical benchmark.

Usage:
    python scripts/reproduce_empirical.py

The command first checks raw-input SHA-256 hashes, regenerates the base empirical
packing tables and the 14-facet numerical benchmark into ``reproduced/``, and then
compares key regenerated quantities with the archived canonical outputs.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


def run(cmd):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run([str(x) for x in cmd], check=True, cwd=ROOT)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    run([sys.executable, SCRIPTS / "check_inputs.py"])
    run([sys.executable, SCRIPTS / "run_scenario.py", "base", "--workers", args.workers])
    run([sys.executable, SCRIPTS / "run_benchmark.py", "--workers", args.workers])

    canonical = pd.read_csv(ROOT / "results" / "base_buildings.csv")
    regenerated = pd.read_csv(ROOT / "reproduced" / "base_buildings.csv")
    c_short = 100 * (1 - canonical.packed_kwp.sum() / canonical.continuous_kwp.sum())
    r_short = 100 * (1 - regenerated.packed_kwp.sum() / regenerated.continuous_kwp.sum())
    delta = abs(c_short - r_short)
    print(f"canonical_shortfall_pct={c_short:.12f}")
    print(f"regenerated_shortfall_pct={r_short:.12f}")
    print(f"absolute_difference_pp={delta:.12g}")
    if delta > 1e-8:
        raise SystemExit("Regenerated base case does not match the archived canonical result.")
    print("Empirical reproduction check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
