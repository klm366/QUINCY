"""
PFT13 Mangrove Trait-Climate Extraction
=======================================

Author:
Dr. Sri Nandini-Weiss

Date:
19 May 2026

Purpose
-------
Extract CHELSA BIO01 and BIO12 climate predictors at TRY observation
coordinates for future SLA, Vcmax25, and Jmax25 regressions.

Scientific Logic
----------------
TRY trait observations are point measurements. To fit climate-trait
relationships, each observation must be paired with local climate.

CHELSA predictors:
- BIO01: mean annual near-surface air temperature, stored as Kelvin x 10
- BIO12: annual precipitation, stored as mm/year

The converted predictors are:
- bio01_temperature_C = BIO01 / 10.0 - 273.15
- bio12_precipitation_mm = BIO12

Assumptions
-----------
- TRY table contains latitude and longitude columns.
- Coordinates are in EPSG:4326.
- CHELSA rasters are in EPSG:4326.
- TRY Request 49898 will later provide the raw trait table.

Limitations
-----------
- This script is prepared before the TRY table is available.
- Column names may need adjustment after inspecting the actual TRY file.
- No regression is fitted here.
- No JSBACH input file is produced here.

Expected Output
---------------
processed/climate_extractions/try_traits_with_chelsa.csv
"""

from pathlib import Path
import sys

import pandas as pd
import rasterio


BASE = Path("/work/bm1548/u241304/preprocessing/pft13_traits")
RAW_TRY = BASE / "raw_try" / "request_49898"
RAW_CLIMATE = BASE / "raw_climate"
OUT_DIR = BASE / "processed" / "climate_extractions"
OUT_FILE = OUT_DIR / "try_traits_with_chelsa.csv"


def require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def find_one(pattern: str, folder: Path, label: str) -> Path:
    matches = sorted(folder.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one {label} matching {pattern} in {folder}, found {len(matches)}"
        )
    return matches[0]


def sample_raster(path: Path, lon, lat):
    with rasterio.open(path) as src:
        if str(src.crs).upper() != "EPSG:4326":
            raise ValueError(f"{path.name} CRS is {src.crs}, expected EPSG:4326")

        values = []
        for value in src.sample(list(zip(lon, lat))):
            values.append(float(value[0]))
    return values


def main() -> int:
    print("PFT13 CHELSA extraction workflow")
    print("--------------------------------")

    bio01 = find_one("*bio01*.tif", RAW_CLIMATE, "CHELSA BIO01 file")
    bio12 = find_one("*bio12*.tif", RAW_CLIMATE, "CHELSA BIO12 file")

    print(f"BIO01 file: {bio01}")
    print(f"BIO12 file: {bio12}")

    try_files = sorted(RAW_TRY.glob("*.csv"))
    if len(try_files) == 0:
        print("\nNo TRY CSV file found yet.")
        print(f"Expected future location: {RAW_TRY}")
        print("Stop here until TRY Request 49898 data are available.")
        return 0

    if len(try_files) > 1:
        raise RuntimeError(f"Expected one TRY CSV file, found {len(try_files)}: {try_files}")

    try_file = try_files[0]
    print(f"TRY file: {try_file}")

    df = pd.read_csv(try_file, low_memory=False)

    print("\nTRY table preview")
    print("-----------------")
    print("Rows:", len(df))
    print("Columns:", list(df.columns))

    lat_candidates = ["Latitude", "latitude", "Lat", "lat"]
    lon_candidates = ["Longitude", "longitude", "Lon", "lon"]

    lat_col = next((c for c in lat_candidates if c in df.columns), None)
    lon_col = next((c for c in lon_candidates if c in df.columns), None)

    if lat_col is None or lon_col is None:
        raise KeyError("Could not identify latitude/longitude columns in TRY table.")

    df = df.dropna(subset=[lat_col, lon_col]).copy()

    df["bio01_raw"] = sample_raster(bio01, df[lon_col], df[lat_col])
    df["bio12_precipitation_mm"] = sample_raster(bio12, df[lon_col], df[lat_col])
    df["bio01_temperature_C"] = df["bio01_raw"] / 10.0 - 273.15

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_FILE, index=False)

    print("\nOutput written:")
    print(OUT_FILE)
    print("Rows:", len(df))

    return 0


if __name__ == "__main__":
    sys.exit(main())
