Future Mangrove Trait Maps
==========================

Author:
Dr. Sri Nandini-Weiss

Date:
19 May 2026

Purpose
-------
Plan future gridded mangrove physiological
trait products for JSBACH.

Planned Outputs
---------------
- bc_mangrove_sla.nc
- bc_mangrove_vcmax25.nc
- bc_mangrove_jmax25.nc

Scientific Logic
----------------
Trait-climate regressions will later be applied
globally across the mangrove mask.

Pipeline
--------
TRY observations
    ↓
CHELSA climate extraction
    ↓
Trait-climate regressions
    ↓
Global prediction surfaces
    ↓
Mangrove mask filtering
    ↓
NetCDF export
    ↓
Future JSBACH ingestion

Assumptions
-----------
- Existing fract_pft13 is the correct mangrove mask.
- Trait maps should be created outside the JSBACH source tree.
- Raw and large NetCDF products should not be committed to Git.

Limitations
-----------
- No JSBACH trait-map reader exists yet.
- No dynamic annual trait update exists yet.
- No salinity or water-table physics is included yet.

Long-Term Goal
--------------
Replace static lctlib fallback traits with
spatially varying and eventually annually
dynamic mangrove physiology.
