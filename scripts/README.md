# Scripts

## Public entry points

- `check_inputs.py` — validates the exact historical analytical inputs by byte size and SHA-256.
- `verify_reported_results.py` — fast manuscript/result crosswalk assertions.
- `run_scenario.py` — reruns one empirical scenario into `reproduced/`.
- `run_benchmark.py` — reruns the submitted 14-facet numerical search benchmark.
- `reproduce_empirical.py` — convenience runner for input check + base empirical rerun + benchmark.
- `make_figures.py` — optionally regenerates analytical Figures 1 and 3–7 locally from canonical tables; generated outputs are intentionally not versioned.

## Internal implementation modules

- `analysis_core.py`
- `final_revision_v5.py`
- `hardening_analysis.py`
- `final_packing_analysis.py`

The internal modules preserve the analysis implementation used during development and are imported by the public entry points. Their historical standalone `main()` routes are deliberately disabled to prevent accidental generation of superseded intermediate outputs. Use the public entry points above.
