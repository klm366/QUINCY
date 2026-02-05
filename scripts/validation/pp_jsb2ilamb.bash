#!/bin/bash
#SBATCH --account=
#SBATCH --job-name=pp_jsb2ilamb
#SBATCH --partition=shared
#SBATCH --time=20:00:00
#SBATCH --output=./pp_jsb2ilamb.%j.log
#SBATCH --error=./pp_jsb2ilamb.%j.log
set -e
# -----------------------------------------------------------------------------------------------------
# This program converts ICON-LAND (JSBACH4) output from one experiment such that it can be evaluated
# using iLamb (see https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2018MS001354)
# -----------------------------------------------------------------------------------------------------

VARLIST='biomass,burntArea,co2,cSoil,dist,evspsbl,gpp,hfdsl,hfls,hfss,lai,mrro,mrsos,nbp,pr,reco,rhums,rlds,rlns,rns,rsds,rsns,swe,tas,tsl,twsa'
TRGGRID='r360x180'
MOL2G_C=12.01
GM2_2_MGHA=0.01

cdo="cdo -s -r -f nc4c"

if [ "$1" == "-h" -o "$1" == "--help" -o ${#@} -eq 0 ] ; then
  echo ""
  echo "Usage: `basename $0` [-h] -p [YYYY-YYYY] outputfiles.nc"
  echo "with:"
  echo "  -e, --expid   experiment name to be used in iLamb"
  echo "  outputfiles   list of all outputfiles (of one stream) to be processed"
  echo "optionally provide:"
  echo "  -m, --mask    bc_land_fract.nc file for land sea mask"
  echo "  -s, --soil    bc_land_soil.nc file for soil depth"
  echo "  -p, --period  start and end year of validation period"
  echo "  -v, --varlist comma separated variable list, e.g. tas,swe,mrsos"
  echo ""
  exit 0
fi

PARACHECK=true
DATA="missing"
EXPID="missing"
PERD="none"
MASK=$(ls /pool/data/ICON/grids/public/mpim/0049/land/r0103/bc_land_frac.nc)
SOIL=$(ls /pool/data/ICON/grids/public/mpim/0049/land/r0103/bc_land_soil.nc)

while $PARACHECK; do
  case $1 in
    -e | --expid   ) EXPID=$2      && shift 2;;
    -m | --mask    ) MASK=$(ls $2) && shift 2;;
    -s | --soil    ) SOIL=$(ls $2) && shift 2;;
    -p | --period  ) PERD=$2       && shift 2;;
    -v | --varlist ) VARLIST=$2    && shift 2;;
    *              ) PARACHECK=false;;
  esac
done

OUTPUTFILES=$@

[[ "$EXPID" == "missing"  ]] && echo "ERROR: No experiment ID provided" && exit 1
[[ "$OUTPUTFILES" == ""   ]] && echo "ERROR: No outputfiles provided"   && exit 1

if [ "$PERD" != "none" ]; then
    PERIOD="-selyear,${PERD%-*}/${PERD#*-}"
else
    PERIOD=""
fi

# -----------------------------------------------------------------------------------------------------
# Create auxiliary data fields
LSM="-ifthen -gtc,0 -selname,notsea $MASK"
SL2="fac_sl2 = 0.035 / max(0.001, min(soil_depth - 0.065, 0.254));" # Currently not needed, as minimum soil depth == 0.5m

# -----------------------------------------------------------------------------------------------------
# Create target directory
printf "\n*** Starting post-processing ${EXPID} for validation with ILAMB\n"
DTRG="$(pwd)/MODELS/${EXPID}"
test -d $DTRG || mkdir -p $DTRG
EXPID=$(echo $EXPID | tr [A-Z] [a-z])

# -----------------------------------------------------------------------------------------------------
# Create command to extract or compute the variables
TASKLIST="$(mktemp -u --tmpdir=.)_tasklist.cdo"
for VAR in $(echo ${VARLIST//,/ }); do
    case $VAR in
        # *** Radiation and Energy Cycle
        hfdsl   ) echo "hfdsl   = sse_grnd_hflx_box * -1;" >> $TASKLIST ;;
        hfls    ) echo "hfls    = seb_latent_hflx_box * -1;"   >> $TASKLIST ;;
        hfss    ) echo "hfss    = seb_sensible_hflx_box * -1;" >> $TASKLIST ;;
        rlds    ) echo "rlds    = a2l_lw_srf_down_box;" >> $TASKLIST ;;
        rlns    ) echo "rlns    = rad_lw_srf_net_box;" >> $TASKLIST ;;
        rns     ) echo "rns     = rad_rad_srf_net_box;" >> $TASKLIST ;;
        rsds    ) echo "rsds    = a2l_swvis_srf_down_box + a2l_swnir_srf_down_box;" >> $TASKLIST ;;
        rsns    ) echo "rsns    = rad_sw_srf_net_box;" >> $TASKLIST ;;
        tas     ) echo "tas     = seb_t_box;"                  >> $TASKLIST ;;
        tsl     ) echo "tsl     = sse_t_soil_sl_box;" >> $TASKLIST ;;
        # *** Hydrology Cycle
        evspsbl ) echo "evspsbl = hydro_evapotrans_box * -1;"             >> $TASKLIST ;;
        mrro    ) echo "mrro    = hydro_runoff_box + hydro_drainage_box;" >> $TASKLIST ;;
        mrsos   ) echo "mrsos   = ((sellevel(hydro_wtr_soil_sl_box,0.0325) + sellevel(hydro_ice_soil_sl_box,0.0325)) + 0.1378 * (sellevel(hydro_wtr_soil_sl_box,0.192) + sellevel(hydro_ice_soil_sl_box,0.192))) * 1000;"       >> $TASKLIST ;;
        pr      ) echo "pr      = a2l_rain_box + a2l_snow_box;"           >> $TASKLIST ;;
        rhums   ) echo "rhums   = disturb_q_air_clim_yday_box;" >> $TASKLIST ;;
        swe     ) echo "swe     = hydro_weq_snow_box;"                    >> $TASKLIST ;;
        twsa    ) echo "twsa    = (vertsum(hydro_wtr_soil_sl_box) + vertsum(hydro_ice_soil_sl_box) + hydro_weq_snow_box) * 1000;" >> $TASKLIST ;;
        # *** Ecosystem and Carbon Cycle
        biomass   ) echo "biomass = carbon_c_sum_veg_ta_box * ${MOL2G_C} * ${GM2_2_MGHA};" >> $TASKLIST ;;
        burntArea ) echo "burntArea = disturb_burned_fract_box * 100;" >> $TASKLIST ;;
        co2       ) echo "co2     = a2l_co2_air_box * 1000 * 1000;"                        >> $TASKLIST ;;
        cSoil     ) echo "cSoil   = (carbon_c_sum_litter_bg_ta_box + carbon_c_sum_humus_ta_box) * ${MOL2G_C} / 1000;" >> $TASKLIST ;;
        dist      ) echo "dist    = disturb_burned_fract_box * fuel_fuel_box * ${MOL2G_C};" >> $TASKLIST ;;
        gpp       ) echo "gpp     = carbon_gpp_ydaymean_ta_box * ${MOL2G_C} * 86400;"      >> $TASKLIST ;;
        lai       ) echo "lai     = pheno_lai_box;"                                        >> $TASKLIST ;;
        nbp       ) echo "nbp     = (carbon_npp_pot_ydaymean_ta_box + carbon_soil_respiration_ta_box) * ${MOL2G_C} * 86400 - carbon_co2_l2a_herb_ta_box * 1000 * 86400;" >> $TASKLIST ;;
        reco      ) echo "reco    = ( - carbon_soil_respiration_ta_box + carbon_gpp_ydaymean_ta_box - carbon_npp_pot_ydaymean_ta_box) * ${MOL2G_C} * 86400 + carbon_co2_l2a_herb_ta_box * 1000 * 86400;" >> $TASKLIST ;;
        # *** Capture missing variables
        *       ) echo "Postprocessing for variable ${VAR} not yet implemented in $0" && exit 1 ;;
    esac
done

# -----------------------------------------------------------------------------------------------------
# Write parameter table for variable metadata
PARTAB="$(mktemp -u --tmpdir=.)_partab.txt"
if [ ! -f $PARTAB ]; then
cat > $PARTAB << EOF
&parameter
  name = biomass
  long_name = "Total vegetation biomass"
  units = "Mg ha-1"
/
&parameter
  name = burntArea
  long_name = "Burnt Area Fraction"
  units = "%"
/
&parameter
  name = co2
  long_name = "Surface Carbon Dioxide Concentration"
  units = "ppm"
/
&parameter
  name = cSoil
  long_name = "Total Soil Carbon"
  units = "kg m-2"
/
&parameter
  name = dist
  long_name = "Wildfire disturbance carbon flux"
  units = "g m-2 month-1"
/
&parameter
  name = evspsbl
  long_name = "Evapotranspiration"
  units = "kg m-2 s-1"
/
&parameter
  name = gpp
  long_name = "Gross primary productivity"
  units = "g m-2 d-1"
/
&parameter
  name = hfdsl
  long_name = "Downward Soil Heat Flux"
  units = "W m-2"
/
&parameter
  name = hfls
  long_name = "Latent heat flux towards atmosphere"
  units = "W m-2"
/
&parameter
  name = hfss
  long_name = "Sensible heat flux towards atmosphere"
  units = "W m-2"
/
&parameter
  name = lai
  long_name = "Leaf area index"
  units = "1"
/
&parameter
  name = mrro
  long_name = "Total runoff"
  units = "kg m-2 s-1"
/
&parameter
  name = mrsos
  long_name = "Surface soil moisture"
  units = "kg m-2"
/
&parameter
  name = nbp
  long_name = "Global Net Ecosystem Carbon Balance"
  units = "g m-2 d-1"
/
&parameter
  name = pr
  long_name = "Total precipitation"
  units = "kg m-2 s-1"
/
&parameter
  name = reco
  long_name = "Net ecosystem respiration"
  units = "g m-2 d-1"
/
&parameter
  name = rhums
  long_name = "2 meter relative humidity"
  units = "%"
/
&parameter
  name = rlds
  long_name = "Surface Downward Longwave Radiation"
  units = "W m-2"
/
&parameter
  name = rlns
  long_name = "Surface Net LW Radiation"
  units = "W m-2"
/
&parameter
  name = rns
  long_name = "Surface Net Radiation"
  units = "W m-2"
/
&parameter
  name = rsds
  long_name = "Surface Downward Shortwave Radiation"
  units = "W m-2"
/
&parameter
  name = rsns
  long_name = "Surface Net SW Radiation"
  units = "W m-2"
/
&parameter
  name = swe
  long_name = "Snow water equivalent"
  units = "m"
/
&parameter
  name = tas
  long_name = "Surface temperature"
  units = "K"
/
&parameter
  name = tsl
  long_name = "Soil Layer Temperature"
  units = "K"
/
&parameter
  name = twsa
  long_name = "Terrestrial water storage anomaly"
  units = "kg m-2"
/
EOF
fi

# -----------------------------------------------------------------------------------------------------
# Loop over output files and compute variables
TEMPID="$(mktemp -u --tmpdir=.)"
PPFILES=""
for FILE in ${OUTPUTFILES}; do
    OUT="${TEMPID}_$(basename ${FILE})"
    printf "Processing $FILE\r"
    $cdo exprf,$TASKLIST $FILE $OUT
    PPFILES="${PPFILES} $OUT"
done



# -----------------------------------------------------------------------------------------------------
# Merge, remap and add metadata
$cdo setpartabn,$PARTAB -remapycon,${TRGGRID} ${LSM} ${PERIOD} -mergetime ${PPFILES} ${TEMPID}_data.nc && rm ${PPFILES}
startdate=($(cdo showdate ${TEMPID}_data.nc))
newdate="$(echo ${startdate[0]} | cut -c 1-7)-01"
cdo -s --reduce_dim -copy ${TEMPID}_data.nc ${TEMPID}_squeezed.nc && rm ${TEMPID}_data.nc && \
       mv ${TEMPID}_squeezed.nc ${TEMPID}_data.nc

# -----------------------------------------------------------------------------------------------------
# Special post-processing steps for individual variables
if [[ $VARLIST =~ "nbp" ]]; then
    $cdo merge -sub -selname,nbp ${TEMPID}_data.nc -divdpm -selname,dist ${TEMPID}_data.nc -delname,nbp ${TEMPID}_data.nc ${TEMPID}_special.nc && \
      mv ${TEMPID}_special.nc ${TEMPID}_data.nc
fi

# -----------------------------------------------------------------------------------------------------
# Rewrite time axis to make iLamb happy
cdo -s -splitvar -setreftime,1850-01-01,00:00:00 -settunits,days -settaxis,${newdate},00:00:00,1month ${TEMPID}_data.nc ${EXPID}_ && rm ${TEMPID}_data.nc
for VAR in $(echo ${VARLIST//,/ }); do
    mkdir -p ${DTRG}/${VAR}
    mv ${EXPID}_${VAR}.nc ${DTRG}/${VAR}/${EXPID}_${VAR}_ilamb.nc
done
printf "\nFinished post-processing for EXPID ${EXPID}\n"


# -----------------------------------------------------------------------------------------------------
# Cleanup
rm $PARTAB $TASKLIST
