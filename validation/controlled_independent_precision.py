"""Sensitivity counterpart: new controlled generator + 10 nm numerical tolerance packing.
Historical published stochastic realization identities remain unavailable.
"""
from pathlib import Path
import controlled_independent as ci
import analysis_core as core
from empirical_precision_audit import robust_count_angle_gap
core.old.count_angle_gap=robust_count_angle_gap
if __name__=='__main__':
 ci.run(repeats=50,seed=20260913,outdir=Path(__file__).parent/'outputs'/'precision_controlled')
