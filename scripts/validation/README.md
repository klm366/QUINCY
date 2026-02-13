# Evaluation of simulation results

While the buildbot scripts conducts technical tests to ensure the quality of the developed model code, scientific model development needs to be evaluated in comparison to real world data. Alternatively to developing analysis on your own, there are a number of benchmarking and simulation analysis tools available for this purpose. Some of those needs some preparation to be setup and used with ICON-LAND. Thus, we collect guidelines for such tools in this section of the wiki.

## ILAMB

The [International LAnd Model Benchmarking](https://ilamb.org/) tool is specifically tailored towards evaluation of land surface states and fluxes. It specializes on comparing simulations and provides scores and plots on how the models perform with respect to each other. An example for an ILAMB analysis can be found [here](https://swift.dkrz.de/v1/dkrz_d7120498-8237-4b35-8a40-439be4e569cb/ICON-DEV/SoilTemperatureFix/index.html).

### Overview and basic use

<details><summary>Setup of ILAMB</summary>
ILAMB can be setup directly on levante. Detailed instructions can be found at their [github page](https://github.com/rubisco-sfa/ILAMB/blob/master/doc/install.rst), but will be summarized in the following paragraphs. The most simply way is to use conda or mamba. Guidelines to setup an conda/mamba environment on levante are [provided by DKRZ](https://docs.dkrz.de/doc/levante/code-development/python.html#). Following the general setup, you need to create a new environment for ILAMB, download the tool and prepare its directory structure:

```bash
mamba create --name ilamb
conda activate ilamb
mamba install ilamb matplotlib=3.8.4
```

Next, you need to create a directory for your ILAMB analysis and the directory structure within

```bash
mkdir <ilamb-directory>
cd <ilamb-directory>
mkdir MODELS
ln -s /pool/data/JSBACH/validation/ILAMB/DATA DATA
export ILAMB_ROOT=$PWD
```

The ILAMB tool uses the `$ILAMB_ROOT` variable as reference to the location of the configure scipt, as well as the `MODELS` and `DATA` directories. Therefore, you either need to `export ILAMB_ROOT=$PWD` (assuming you are in the directory) prior to any analysis in a new terminal, or - even better - add an `export` to your `.bashrc`

ILAMB provides a large set of [observational dataset](https://github.com/rubisco-sfa/ILAMB-Data). While the data can be downloaded or updated running `ilamb-fetch` in your `ILAMB_ROOT` directory, you should instead link to the data already downloaded in the DKRZ data pool at `/pool/data/JSBACH/validation/ILAMB/DATA`.

</details>

ILAMB requires a fixed directory structure for observation data as well as simulation data. While the observation data structure is created by ILAMB itself, the user needs to take care about simulation data postprocessing. ILAMB expects only one variable per files structured as in this example:
```
MODELS
  |------ Sim1
  |        |------ hfls
  |        |        |------ sim1_hfls.nc
  |        |
  |        |------ hfss
  |                 |------ sim1_hfss.nc
  |------ Sim2
  |        |------ hfls
  |        |        |------ sim2_hfls.nc
  |        |
  |        |------ hfss
  |                 |------ sim2_hfss.nc
```
which needs to be extended to all simulations and variables that should be part of the evaluation. The script `postproc_jsbach-for-ilamb.sh` can extract the variables from each simulation and sort them into this directory structure. It needs to get at least an experiment name and the path to the `*lnd_basic*.nc` output stream files. Please use `./postproc_jsbach-for-ilamb.sh --help` to get an overview about its optional parameters. Those need to be adapted if you either analyse simulations which don't use the 0049-R2B4 grid or want to only use a subset of variables or timesteps.

As next step, a configure script needs to be create to tell ILAMB which analysis to perform. It's basic structure follows this syntax
```yaml
[h1: Hydrology Cycle]
bgcolor = "#E6F9FF"

[h2: Evapotranspiration]
variable="et"
alternate_vars = "evspsbl"
[MOD16A2]
source="DATA/evspsbl/MOD16A2/et.nc"
[GLEAMv3.3a]
source="DATA/evspsbl/GLEAMv3.3a/et.nc"

[h2: Snow Water Equivalent]
variable="swe"
[CanSISE]
source="DATA/swe/CanSISE/swe.nc"
```
Different topics are set using the `h1` level, while individual variables use the `h2` level. The entries `variable` and `alternate_vars` provide names which are used in the simulation and observation datasets, the respective dataset names for comparison are given as `[DATASET_NAME]` and `source=` needs to contain the path to the dataset. Further variable options are:
- `weight`
- `cmap`
- `mass_weighting`
The file `ilamb-for-jsbach.cfg` contains all currently supported variables.

After providing all this information, ILAMB is run as
```bash
ilamb-run --config <my_config.cfg> --model_root $ILAMB_ROOT/MODELS/ --regions global
```
The analysis results are stored in `$ILAMB_ROOT/_build`. The most convenient way is to open it via web-browser (which mean you need to mount the directory locally, if you used it on levante). Setup an server with `python -m http.server` and go to the page provided by python (usually `http://0.0.0.0:8000/`)

A ICON-LAND post-processing script suitable for ILAMB as well as an ICON-LAND config file are currently under preparation and will be added to the scripts folder of the ICON-LAND repository.

### Simulation data post-processing

:warning: under construction

**Radiation and Energy cycle**

| Variable | ILAMB name | ILAMB unit | ICON-LAND post-processing |
| ------ | ------ | ------ | ------ |
| 2m temperature     | tas    | K     | `seb_t_box` |
| Ground heat flux | hfdsl | W m-2 | `sse_grnd_hflx_box` * -1 |
| Latent heat flux   | hfls   | W m-2 | `seb_latent_hflx_box` * -1 |
| Sensible heat flux | hfss   | W m-2 | `seb_sensible_hflx_box` * -1 |
| Soil Layer Temperature | tsl | K | `sse_t_soil_sl_box` |
| Surface Downward Longwave Radiation | rlds | W m-2 | `a2l_lw_srf_down_box` |
| Surface Downward Shortwave Radiation | rsds | W m-2 | `a2l_swvis_srf_down_box + a2l_swnir_srf_down_box` |
| Surface Net LW Radiation | rlns | W m-2 | `rad_lw_srf_net_box` |
| Surface Net Radiation | rns | W m-2 | `rad_rad_srf_net_box` |
| Surface Net SW Radiation | rsns | W m-2 | `rad_sw_srf_net_box` |

**Hydrological cycle**

| Variable | ILAMB name | ILAMB unit | ICON-LAND post-processing |
| ------ | ------ | ------ | ------ |
| 2 meter relative humidity | rhums | % | `disturb_q_air_clim_yday_box` |
| Evapotranspiration | evspsbl | kg m-2 s-1 | `hydro_evapotrans_box` * -1 |
| Precipitation | pr | kg m-2 s-1 | `a2l_rain_box` + `a2l_snow_box` |
| Snow water equivalent | swe | m weq | `hydro_weq_snow_box` |
| Surface soil moisture | mrsos | kg m-2 | `hydro_wtr_soil_sl_box` + `hydro_ice_soil_sl_box`, 1st lay + 0.1378 * 2nd lay (min soil depth == 0.5m --> no soil depth factor needed) |
| Terrestrial Water storage anomaly | twsa | kg m-2 | (vertsum(`hydro_wtr_soil_sl_box`) + vertsum(`hydro_ice_soil_sl_box`) + `hydro_weq_snow_box`) * 1000 |
| Total runoff | mrro | kg m-2 s-1 | `hydro_runoff_box + hydro_drainage_box` |

**Ecosystem and carbon cycle**

| Variable | ILAMB name | ILAMB unit | ICON-LAND post-processing |
| ------ | ------ | ------ | ------ |
| Total vegetation biomass | biomass | [Mg(C) ha-1] | `carbon_c_sum_veg_ta_box` * 12.01 * 0.01 |
| Burnt Area Fraction | burntArea | % | `disturb_burned_fract_box` * 100 |
| Surface Carbon Dioxide Concentration | co2 | ppm | `a2l_co2_air_mol_box` * 1.0e+6 |
| Total Soil Carbon | cSoil | kg m-2 | (`carbon_c_sum_litter_bg_ta_box` + `carbon_c_sum_humus_ta_box`) * 12.01 / 1000 |
| Gross primary productivity | gpp | [g(C) m-2 d-1] | `carbon_gpp_ydaymean_ta_box` * 12.01 * 86400 |
| Wildfire disturbance carbon flux | dist | g(C) m-2 month-1 | `disturb_burned_fract_box` * `fuel_fuel_box` * 12.01 |
| Leaf area index | lai | [m2 m-2] | `pheno_lai_box` |
| Global Net Ecosystem Carbon Balance | nbp | g(C) -m-2 d-1 | ((`carbon_npp_pot_ydaymean_ta_box` + `carbon_soil_respiration_ta_box`) * 12.01 * 86400 - `carbon_co2_l2a_herb_ta_box` * 1000 * 8640) / dist |
| Net ecosystem respiration | reco | g(C) -m-2 d-1 | ( -`carbon_soil_respiration_ta_box` + `carbon_gpp_ydaymean_ta_box` - `carbon_npp_pot_ydaymean_ta_box`) * 12.01 * 86400 + `carbon_co2_l2a_herb_ta_box` * 1000 * 86400|

### Additional post-processing within ILAMB

*Add variables that are computed within ILAMB*

### Troubleshooting

ILAMB is a very sensible tool and might fail easily if the model data does not fit its expectation perfectly. Here are a number of possible solutions to encountered problems:

- `VarNotInModel`: If this happens during the `Running model-confrontation pairs` step, the most likely reason is that the time period of your simulation data has no overlap with the time period of the ilamb observations
- `IndexError`: ?
- `ValueError`: ?
- `AssertionError`: ?
- Whatever else error: first, try to delete the `_build` directory and run again. ILAMB stores data from previous runs which often causes conflicts
