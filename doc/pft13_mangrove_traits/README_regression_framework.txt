PFT13 Mangrove Trait Regression Framework
=========================================

Author:
Dr. Sri Nandini-Weiss

Date:
19 May 2026

Purpose
-------
Develop climate-dependent mangrove trait regressions
for future integration into JSBACH.

Scientific Motivation
---------------------
Mangrove physiology varies strongly across:
- coastal temperature gradients
- precipitation gradients
- tropical hydroclimate conditions

Static lctlib values cannot represent this variability.

Therefore:
SLA, Vcmax25, and Jmax25 will later become
climate-dependent physiological traits.

Scientific Basis
----------------
This framework follows the conceptual logic of:

- Verheijen et al. (2013, 2015)
- trait-enabled DGVM approaches
- ORCHIDEE-MAN style mangrove parameterization

Current Workflow
----------------
1. Stable structural PFT13 implementation completed
2. Mangrove mask already implemented using GMW
3. CHELSA climate predictors downloaded
4. TRY Request 49898 submitted
5. Regression preprocessing framework created

Climate Predictors
------------------
CHELSA BIO01:
Mean annual temperature

CHELSA BIO12:
Annual precipitation

Trait Targets
--------------
- Specific leaf area (SLA)
- Vcmax25
- Jmax25

Future Workflow
---------------
1. Receive TRY tables
2. Clean trait observations
3. Extract climate at TRY coordinates
4. Merge climate + traits
5. Fit regressions
6. Generate gridded trait maps
7. Integrate into JSBACH
8. Later evolve toward dynamic annual traits

Limitations
-----------
Current framework:
- uses static climatologies
- does not yet include salinity
- does not yet include water-table physics
- currently offline JSBACH only

Future Extensions
-----------------
- Dynamic annual traits
- Salinity stress
- Water-table coupling
- Coupled ICON-A simulations
- Future ICON-ESM integration
