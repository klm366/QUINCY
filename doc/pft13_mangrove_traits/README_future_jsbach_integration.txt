Future JSBACH Integration Plan
==============================

Author:
Dr. Sri Nandini-Weiss

Date:
19 May 2026

Purpose
-------
Document the planned transition from static PFT13 mangrove
parameters toward spatially varying and eventually annually
dynamic physiological traits in JSBACH.

Scientific Logic
----------------
The current PFT13 implementation is structurally complete but
physiologically static.

Current trait values are still read from:
data/lctlib_nlct22_mangrove.def

The next scientific goal is to replace fixed fallback trait values
with climate-dependent trait fields derived from TRY and CHELSA.

Target Traits
-------------
- SLA
- Vcmax25
- Jmax25

Current JSBACH Pathways
-----------------------
CarboxRate:
lctlib -> mo_assimi_interface.f90 -> mo_assimi_process.f90 -> VCMAX

ETransport:
lctlib -> mo_assimi_interface.f90 -> mo_assimi_process.f90 -> JMAX

specificLeafArea_C:
lctlib -> carbon, radiation, and phenology pathways

Implementation Stages
---------------------
Stage 1:
Static lctlib fallback cleanup for PFT13.

Stage 2:
Static gridded NetCDF trait maps.

Stage 3:
Runtime JSBACH trait-map ingestion.

Stage 4:
Annual climate-dependent trait updates.

Stage 5:
Salinity and water-table process coupling.

Limitations
-----------
- No NetCDF trait reader implemented yet.
- No dynamic annual update implemented yet.
- No salinity or water-table stress implemented yet.
- Current offline JSBACH simulations use static fallback traits.

Expected Future Outputs
-----------------------
- bc_mangrove_sla.nc
- bc_mangrove_vcmax25.nc
- bc_mangrove_jmax25.nc
- restart-safe trait diagnostics
- future salinity and water-table stress diagnostics
