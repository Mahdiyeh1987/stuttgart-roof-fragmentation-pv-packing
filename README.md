# Roof-facet fragmentation and rooftop PV packing in Stuttgart

[![Python 3.13](https://img.shields.io/badge/Python-3.13-blue.svg)](https://www.python.org/)
[![Code license: MIT](https://img.shields.io/badge/code%20license-MIT-green.svg)](LICENSE)

Replication data, analysis code, and processed results for:

> **Quantifying Roof-Facet Fragmentation Effects on Rooftop PV Packing Using LoD2 Geometry: A Controlled Experiment and Stuttgart Pilot**  
> Mahdiyeh Tabatabaei and Jan Cremers  
> Submitted to *Electricity* (2026).

## What this repository reproduces

The study separates two related geometric effects:

1. a **controlled fragmentation experiment**, in which the outer roof geometry remains fixed while non-crossable internal subdivisions are added; and
2. an **empirical Stuttgart LoD2 pilot**, in which complete rectangular PV modules are packed independently on audited physical roof facets.

The submitted base case contains **31 buildings and 340 physical roof facets**. With a 0.30 m external planimetric setback, a LONGi LR5-54HTB-440M module, and a 10 mm inter-module clearance, the archived result is:

| Quantity | Submitted value |
|---|---:|
| Setback-adjusted roof-plane area | 63,587.5 m² |
| Continuous-area capacity reference | 14.328 MWp |
| Explicitly packed capacity | 12.284 MWp |
| Aggregate packing shortfall | 14.26% |
| 95% bootstrap percentile interval | 12.67–16.60% |
| Spearman ρ: physical-facet density vs. shortfall | 0.631 |
| FDR-adjusted p | 0.0010 |

The empirical test window is a bounded historical Stuttgart example and is **not a city-representative probability sample**. Irradiance, shading, rooftop obstacles, structural adequacy, electrical design, energy yield, and economics are outside the scope of the reported packing correction.

## Repository structure

```text
.
├── data_raw/                 Historical Stuttgart input files and provenance notes
├── results/                  Canonical processed tables used in the manuscript
├── scripts/                  Reproduction, verification, and figure-generation code
├── supplementary/            Supplementary Materials submitted with the manuscript
├── docs/                     Reproducibility, provenance, and result crosswalk documentation
├── .github/workflows/        Lightweight integrity checks for GitHub Actions
├── CITATION.cff              Citation metadata
├── .zenodo.json              Zenodo deposit metadata
├── environment.yml           Reproducible environment specification
└── requirements.txt          Exact Python package versions
```

## Quick verification

Create an isolated Python 3.13 environment and install the exact dependencies:

```bash
python -m venv .venv
```

Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

macOS/Linux:

```bash
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Then run the fast integrity checks:

```bash
python scripts/check_inputs.py
python scripts/verify_reported_results.py
python scripts/run_robustness_checks.py --skip-coplanar  # fast final robustness checks
python scripts/make_figures.py  # optional; generates figures locally into an ignored figures/ folder
```

`check_inputs.py` verifies the exact historical analytical inputs using their byte sizes and SHA-256 hashes. `verify_reported_results.py` checks the archived result tables against the numerical values reported in the submitted manuscript, including the final size-adjusted robustness checks. `run_robustness_checks.py` regenerates the two HC3 models, the residual-rank partial Spearman check, and the seven-pair coplanar-tolerance sensitivity grid. `make_figures.py` can regenerate analytical Figures 1 and 3–7 locally from repository-relative inputs. Generated figure files are intentionally excluded from this public repository. Figure 2 is not redistributed because it contains contextual municipal aerial imagery and is not a computational input.

## Re-run the empirical analysis from raw inputs

To regenerate the empirical base case and numerical packing benchmark without overwriting the canonical archived tables:

```bash
python scripts/reproduce_empirical.py
```

Regenerated files are written to the ignored `reproduced/` directory. Individual scenarios can also be run, for example:

```bash
python scripts/run_scenario.py base
python scripts/run_scenario.py gap04
python scripts/run_benchmark.py
```

The empirical pipeline performs CityGML parsing, geometry repair, nested-coplanar exclusion, adjacent-coplanar consolidation, roof-plane transformation, setback clipping, complete-module packing, and aggregation. See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for the exact reproducibility scope and entry points.

## Final robustness checks

Version 1.1.0 adds the robustness outputs reported in the final manuscript: (i) HC3-robust OLS using standardized facet density plus log gross 3D roof area; (ii) the complementary raw physical-facet-count plus log-area model; (iii) a residual-rank partial Spearman check controlling for ranked log gross 3D roof area; and (iv) the seven-pair coplanar-consolidation tolerance grid. Canonical outputs are archived under `results/`, and `python scripts/run_robustness_checks.py` regenerates them.

No publication image files are distributed in this repository. Figure-generation code is retained only for reproducibility and writes to a Git-ignored local directory.

## Controlled fragmentation experiment

The reported controlled experiment uses a 12 × 10 m outer roof, a 0.30 m external setback, 50 pseudo-random realizations per facet count, facet counts 2–12, seed `20260913`, a minimum accepted split piece of 0.8 m², and the same 10 mm-clearance production packing rule used in the empirical calculation. Its canonical publication outputs are archived in:

- `results/controlled_production_summary.csv`
- `results/controlled_production_correlations.csv`

The recovered source bundle contains the final controlled outputs and full methodological specification but not the exact historical realization-level generator source. The repository therefore **verifies and regenerates the controlled publication figure from the archived canonical tables**, while the empirical CityGML-to-capacity chain and numerical benchmark are directly executable from raw inputs. This boundary is documented explicitly rather than silently reconstructing an unverified generator.

## Data provenance

The empirical analysis uses an official Stuttgart historical test-data package exported on **27 November 2018**, in ETRS89 / UTM zone 32N (`EPSG:25832`). The package includes the exact filenames, byte sizes, and SHA-256 checksums used in the study.

Key upstream sources and current documentation:

- Landeshauptstadt Stuttgart, Stadtmessungsamt — Open Geodata and test data: https://www.stuttgart.de/leben/bauen/geoportal/open-data-und-testdaten
- Stuttgart ALKIS OpenData: https://opendata.stuttgart.de/dataset/alkis
- LGL Baden-Württemberg Open GeoData: https://www.lgl-bw.de/Produkte/Open-Data/
- OGC CityGML 2.0: https://www.ogc.org/standards/citygml/

See [`docs/DATA_PROVENANCE.md`](docs/DATA_PROVENANCE.md) and [`THIRD_PARTY_DATA.md`](THIRD_PARTY_DATA.md).

## Software environment

The submitted calculations were carried out with:

- Python 3.13.5
- NumPy 2.3.5
- pandas 2.2.3
- Shapely 2.1.2
- SciPy 1.17.0
- statsmodels 0.14.6
- Matplotlib 3.10.8
- lxml 6.1.1

Exact pins are in `requirements.txt` and `environment.yml`.

## Reuse and licensing

The authors' **source code and original repository documentation** are released under the MIT License; see [`LICENSE`](LICENSE).

Third-party geodata, municipal imagery, product information, and other source materials retain their original upstream terms and are **not relicensed under MIT**. In particular, LGL open geodata requires the attribution:

> Datenquelle: LGL, www.lgl-bw.de, dl-de/by-2-0

See [`THIRD_PARTY_DATA.md`](THIRD_PARTY_DATA.md) for source-data provenance and redistribution notes. Publication figure files are intentionally not included in this repository.

## Citation

Use the repository citation metadata in [`CITATION.cff`](CITATION.cff). For the final manuscript, cite the Zenodo DOI corresponding to the exact archived release used at submission. Version 1.1.0 supersedes 1.0.0 for the final robustness-expanded package.

## Authors

- **Mahdiyeh Tabatabaei** — Department of Architecture, Alma Mater Studiorum – University of Bologna, Bologna, Italy
- **Jan Cremers** — Faculty of Architecture and Design, HFT Stuttgart – University of Applied Sciences, Stuttgart, Germany

## Contact

Correspondence: `mahdiyeh.tabatabaei2@unibo.it`
