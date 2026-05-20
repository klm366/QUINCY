"""
PFT13 Mangrove Trait Regression Framework
-----------------------------------------

Author:
Dr. Sri Nandini-Weiss

Date:
19 May 2026

Purpose
-------
Prepare CHELSA climate rasters for extraction at TRY observation
coordinates for future mangrove trait-climate regressions.

Scientific Logic
----------------
Mangrove physiology traits vary with environmental gradients,
particularly temperature and precipitation.

This workflow follows the conceptual logic of:
- Verheijen et al. trait-enabled JSBACH approaches
- trait-climate regression frameworks
- ORCHIDEE-MAN style parameterization concepts

TRY observations provide:
- SLA
- Vcmax25
- Jmax25

CHELSA provides:
- BIO01 temperature
- BIO12 precipitation

These climate rasters will later be sampled at each TRY coordinate.

Current Stage
-------------
Infrastructure preparation only.

No TRY extraction is performed yet because:
- TRY Request 49898 is still pending approval.

Future Workflow
---------------
1. Read TRY tables
2. Extract CHELSA climate at TRY coordinates
3. Merge climate + traits
4. Fit regressions
5. Create gridded trait maps
6. Integrate into JSBACH

Limitations
-----------
- Currently static climatologies only
- No dynamic annual climate yet
- No salinity or water-table physics yet
- No coupled ICON simulations yet
"""

from pathlib import Path
import rasterio

BASE = Path("/work/bm1548/u241304/preprocessing/pft13_traits")

bio01 = next((BASE / "raw_climate").glob("*bio01*.tif"))
bio12 = next((BASE / "raw_climate").glob("*bio12*.tif"))

print("\nPFT13 climate sampling preparation")
print("----------------------------------")

for path in [bio01, bio12]:

    with rasterio.open(path) as src:

        print("\nFILE:", path.name)
        print("CRS:", src.crs)
        print("Resolution:", src.res)
        print("Bounds:", src.bounds)
        print("Width:", src.width)
        print("Height:", src.height)

print("\nClimate rasters ready for future TRY coordinate extraction.")
