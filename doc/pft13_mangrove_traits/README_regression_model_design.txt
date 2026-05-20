PFT13 Regression Model Design
=============================

Author:
Dr. Sri Nandini-Weiss

Date:
19 May 2026

Purpose
-------
Define the statistical framework for mangrove trait-climate regressions.

Scientific Logic
----------------
Mangrove physiological traits vary with temperature, precipitation, and coastal hydroclimate gradients.
The first-order regression framework links observed TRY trait data to CHELSA climatological predictors.

Scientific Basis
----------------
- Verheijen et al. trait-enabled JSBACH framework
- Trait-enabled DGVM concepts
- ORCHIDEE-MAN mangrove parameterization logic

Target Traits
-------------
- SLA
- Vcmax25
- Jmax25

Climate Predictors
------------------
CHELSA BIO01:
Mean annual temperature, stored as Kelvin x 10.

CHELSA BIO12:
Annual precipitation, stored as mm/year.

Initial Statistical Form
------------------------
log(SLA)     = b0 + b1 * BIO01_temperature_C + b2 * BIO12_precipitation_mm
log(Vcmax25) = b0 + b1 * BIO01_temperature_C + b2 * BIO12_precipitation_mm
log(Jmax25)  = b0 + b1 * BIO01_temperature_C + b2 * BIO12_precipitation_mm

Assumptions
-----------
- TRY observations contain usable latitude and longitude.
- CHELSA BIO01 and BIO12 are appropriate first-order climate predictors.
- Static climatology is acceptable for the first offline JSBACH implementation.
- Mangrove filtering will be applied after the broad global trait table is downloaded.

Limitations
-----------
- No TRY table has been received yet.
- No regression has been fitted yet.
- No salinity predictor is included yet.
- No water-table predictor is included yet.
- No dynamic annual update is implemented yet.

Expected Outputs
----------------
- merged TRY-CHELSA regression table
- regression coefficients for SLA, Vcmax25, Jmax25
- diagnostic plots
- future gridded trait maps
