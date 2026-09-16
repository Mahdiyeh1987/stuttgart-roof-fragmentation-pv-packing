# Data provenance and lineage

## Historical analytical inputs

The empirical calculation uses a Stuttgart test-data export dated **27 November 2018**. The exact analytical files are retained in `data_raw/` and identified by SHA-256 in `results/input_file_manifest_sha256.csv`.

| File | Role |
|---|---|
| `CityGML_Alles_2018_11_27_07_15_56.zip` | Historical CityGML / LoD2 roof geometry |
| `2018-11-27_07-16-14_Liegenschaftskarte_ohne_Eigentuemer.xml` | Historical ALKIS/NAS geometry without owner information |
| `Auswertung.log` | Supplied export summary |
| `Koordinatenreferenzsystem_EPSG_25832_ETRS89_UTM.txt` | Supplied CRS note |
| `Information_zum_NAS_Schema.txt` | Supporting NAS-schema information; not an analytical hash target |

All metric geometry calculations use **ETRS89 / UTM zone 32N (EPSG:25832)**.

## Geometry lineage

The final empirical analytical path is:

```text
raw CityGML roof polygon
→ validity / sliver audit
→ nested-coplanar exclusion
→ adjacent-coplanar consolidation
→ physical roof facet
→ external-envelope setback clipping
→ local roof-plane coordinates
→ complete-module packing
→ building aggregation
→ statistical analysis
```

Key lineage counts in the submitted analysis:

- 450 raw CityGML roof polygons
- 431 retained source polygons after the geometry audit and nested-coplanar exclusion
- 340 physical roof facets after adjacent-coplanar consolidation
- 31 CityGML buildings with explicit roof geometry
- 58 `AX_Gebaeude` footprint features in the supplied NAS extract

The 94 `Gebäude` records in the export summary, 58 NAS features, and 31 CityGML buildings refer to different object definitions and are not sequential filtering-stage counts.

## Current upstream documentation

The City of Stuttgart currently documents open geodata and test datasets, including a CityGML 3D-city-model test dataset, on its Geoportal page:

https://www.stuttgart.de/leben/bauen/geoportal/open-data-und-testdaten

The current Stuttgart ALKIS OpenData entry is:

https://opendata.stuttgart.de/dataset/alkis

The current LGL Baden-Württemberg Open GeoData information is:

https://www.lgl-bw.de/Produkte/Open-Data/

These current pages are documentation/source routes. The numerical analysis itself retains the exact historical 2018 files listed above.

## Integrity verification

Run:

```bash
python scripts/check_inputs.py
```

The script compares each analytical input against both the expected byte size and SHA-256 checksum in `results/input_file_manifest_sha256.csv`.
