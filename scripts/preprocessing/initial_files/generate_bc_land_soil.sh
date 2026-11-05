#!/bin/ksh

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2026, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------
set -e
#-----------------------------------------------------------------------------
#   Generate
#     ICON-Land soil boundary conditions file
#-----------------------------------------------------------------------------
#
# The script is part of a script series to generate jsbach initial files
#   started by master script "create_jsbach_ini_files.sh".
#
#   bc_land_soil:
#   - Soil textures from extpar
#   - Soil properties from FAO
#   - Soil organic carbon fraction based on WISE data
#   - Deprecated variable maxmoist calculated.

# Variables exported from master script create_icon-land_ini_files.sh
## extpar
extpar_file=${extpar_file}
icon_grid=${icon_grid}
## grid
atmGridID=${atmGridID}
## directories
path_bc=${bc_file_dir}
work_dir=${work_dir}
scripts_dir=${scripts_dir}

# Variables exported from script create_icon-land_ic_bc_files.sh
pool_prepare=${pool_prepare}

# Directory for input data specific to IQJ (QUINCY)
quincyInputDataDirPath="/work/mj0143/icon_land_data_pool/preprocessing/bc_soil_quincy"

# set defaults for major commandline tools
cdo=${cdo}
rm=$rm
cp=$cp
# set QUINCY specific cdo command options (silent, create no history and use double precision)
cdoQ="cdo -s --no_history -b 64"

# Local variables
prog=$(basename $0)
file_bc_soil=bc_land_soil
clean_up=true  # Remove intermediate files


cd ${work_dir}

# ----------------------------------------------------------------------------
#  Soil data file
# -----------------------------
echo "${prog}: Generating ${file_bc_soil}:"

# The order of variables should not be changed
var_list="FR_SAND FR_SILT FR_CLAY SUB_FR_SAND SUB_FR_SILT SUB_FR_CLAY \
 root_depth maxmoist soil_depth fract_org_sl bclapp bclapp_mineral soil_field_cap \
 soil_field_cap_mineral heat_capacity heat_capacity_mineral heat_conductivity \
 heat_conductivity_mineral hyd_cond_sat hyd_cond_sat_mineral moisture_pot \
 moisture_pot_mineral pore_size_index pore_size_index_mineral soil_porosity \
 soil_porosity_mineral wilting_point wilting_point_mineral residual_water \
 residual_water_mineral fao BLDFIE PHIHOX AWCtS \
 Labile_Inorganic_P Seconday_Mineral__P Occluded_P Apatite_P TAXOUSDA_10km TAXNWRB qmax_org"

# Preparations that need to be done outside the parallel loop
${cdo} splitvar ${pool_prepare}/05/soil_parameter_05_apr2024.nc   05_

# Process root depth, as it is needed for soil depth calculation
if [[ ${var_list} == *"soil_depth"* || ${var_list} == *"root_depth"* ]]; then
  # ----------------------
  var=root_depth
  # ----------------------
  echo "   ${var} ..."
  minimum_root_depth=0.15  # Note: root depth of GLCC data used in extpar ranges from 0.3 to 2.0 m.
  # Get root depth from extpar file and extrapolate land values to ocean area
  ${cdo} -setmisstodis -setmissval,0 -selvar,ROOTDP ${extpar_file} extpar_ROOTDP.nc
  # Set minimum root depth - also use this value for complete ocean cells. Root depth of glaciers is zero.
  ${cdo} -setattribute,root_depth@long_name="Root depth" -setattribute,root_depth@units="m" \
      -setname,root_depth -mul notglac.nc -maxc,${minimum_root_depth} \
      -mul non-glac-land-mask.nc extpar_ROOTDP.nc  ${atmGridID}_root_depth.nc
  # Copy for later usage in generate_ic_land_soil
  $cp ${atmGridID}_root_depth.nc ${atmGridID}_root_depth.keep
fi

# averaging of input values for input variables: bulk density, ph, awc
python ${scripts_dir}/generate_bc_land_soil_hlp.py ${quincyInputDataDirPath} ${work_dir}
cp ${scripts_dir}/generate_bc_land_soil_hlp.py ${path_bc}/scripts/
echo ">> calculations (weighted average, medians in grid cells etc) were successfully performed by generate_bc_land_soil_hlp.py"

# main loop over var_list (runs variable processing in parallel, see '&' after closing '}')
for var in ${var_list}; do
  {
  echo "   ${var} ..."
  case ${var} in
  # ----------------------
    FR_SAND | FR_SILT | FR_CLAY )
  # ----------------------
      # Soil texture fractions in the upper 30 cm of the soil column
      # Convert from percent to fraction and remove 'institution' attribute (i.e. DWD)
      ${cdo} mulc,0.01 -selvar,${var} ${extpar_file} ${atmGridID}_${var}.tmp
      ncatted -a institution,global,d,c,sng       ${atmGridID}_${var}.tmp

      # Generate mask for valid data, used below to mask out grid cells
      # where any of the texture variables are missing
      ${cdo} gec,0. ${atmGridID}_${var}.tmp ${atmGridID}_${var}_MASK.nc
      ;;
  # ----------------------
    SUB_FR_SAND | SUB_FR_SILT | SUB_FR_CLAY )
  # ----------------------
      # Soil texture fractions in deeper soil (below 30 cm)
      # Convert from percent to fraction and remove institution attribute (i.e. DWD)
      ${cdo} mulc,0.01 -selvar,${var} ${extpar_file} ${atmGridID}_${var}.tmp
      ncatted -a institution,global,d,c,sng       ${atmGridID}_${var}.tmp

      # Generate mask for valid data, used below to mask out grid cells
      # where any of the texture variables are missing
      ${cdo} gec,0. ${atmGridID}_${var}.tmp ${atmGridID}_${var}_MASK.nc
      ;;
  # ----------------------
    root_depth )
  # ----------------------
      # Root depth already calculated above ...
      ;;
  # ----------------------
    maxmoist )
  # ----------------------
      # Maxmoist calculated below ...
      ;;
  # ----------------------
    soil_depth )
  # ----------------------
      # In case soil depth is smaller than root depth it needs to be adapted: As the root depth
      # data set in extpar is more recent than the soil data from echam we rather adapt soil
      # depth to root depth - and not vice versa as we did in previous versions.
      ${cdo} selvar,soildepth ${pool_prepare}/05/soil_parameters_05.nc 05_soil_depth.nc
      ${cdo} remap,${icon_grid},wgt05.nc -setmisstodis 05_soil_depth.nc \
          ${atmGridID}_soil_depth.tmp
      ${cdo} -setattribute,soil_depth@long_name="Soil depth until bedrock" -setattribute,soil_depth@units="m" \
          -setvar,soil_depth -setmisstoc,0. -ifthen non-glac-land-mask.nc \
          -max ${atmGridID}_soil_depth.tmp ${atmGridID}_root_depth.nc ${atmGridID}_soil_depth.nc
      ncatted -a code,${var},d,c, ${atmGridID}_${var}.nc
      # Copy for later usage in generate_ic_land_soil
      $cp ${atmGridID}_soil_depth.nc ${atmGridID}_soil_depth.keep
      ;;
  # ----------------------
    fract_org_sl )
  # ----------------------
      wise_data=${pool_prepare}/05/wise_som-fract_05deg_jsbsoillayers.nc
      ${cdo} -selvar,${var} -sellonlatbox,0,360,-90,90 $wise_data 05_${var}.tmp
      # Extrapolation
      ${cdo} -setmisstodis -ifthen 05_non-glac-land.nc 05_${var}.tmp 05_${var}.nc
      ${cdo} remap,${icon_grid},wgt05.nc 05_$var.nc  ${atmGridID}_${var}.tmp
      ${cdo} -setmisstoc,0 -ifthen non-glac-land-mask.nc ${atmGridID}_${var}.tmp \
          ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    bclapp | bclapp_mineral | soil_field_cap | soil_field_cap_mineral | heat_capacity | heat_capacity_mineral \
    | heat_conductivity | heat_conductivity_mineral | hyd_cond_sat | hyd_cond_sat_mineral | moisture_pot \
    | moisture_pot_mineral | pore_size_index | pore_size_index_mineral | soil_porosity | soil_porosity_mineral \
    | wilting_point | wilting_point_mineral | residual_water | residual_water_mineral )
  # ----------------------
      ${cdo} -setvar,${var} -mul -remap,${icon_grid},wgt05.nc -invertlat -setmisstodis 05_${var}.nc \
                                 notglac.nc    ${atmGridID}_${var}.tmp
      case ${var} in
        bclapp* )
          ${cdo} -setmisstoc,4.5    -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        soil_field_cap* )
          ${cdo} -setmisstoc,0.229  -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          # Copy for later usage in ic_soil script
          $cp ${atmGridID}_${var}.nc ${atmGridID}_${var}.keep
          ;;
        heat_capacity* )
          ${cdo} -setmisstoc,2.e+6  -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        heat_conductivity* )
          ${cdo} -setmisstoc,7.     -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        hyd_cond_sat* )
          ${cdo} -setmisstoc,5.e-6  -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        moisture_pot* )
          ${cdo} -setmisstoc,-0.15  -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        pore_size_index* )
          ${cdo} -setmisstoc,0.2    -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        soil_porosity* )
          ${cdo} -setmisstoc,0.45   -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        wilting_point* )
          ${cdo} -setmisstoc,0.15   -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        residual_water* )
          ${cdo} -setmisstoc,0.05   -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        * )
          ${cdo} -setmissval,-9.e33 -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
      esac
      ;;
  # ----------------------
    fao )
  # ----------------------
      # We use old ECHAM5 data here! TODO: replace
      ${cdo} -remaplaf,${icon_grid} -selvar,FAO ${pool_prepare}/05/ECHAM6/05_jan_surf.nc \
          ${atmGridID}_fao.tmp
      ${cdo} setname,fao -setmisstoc,0 -ifthen non-glac-land-mask.nc ${atmGridID}_fao.tmp \
          ${atmGridID}_fao.nc
      ;;
  # ----------------------
    BLDFIE )
  # ----------------------
      # bulk density: SPREAD(sse_init_vars%bulk_density) -> bulk_dens_sl
      # *.tmp file provided by generate_bc_land_soil_hlp.py
      # default value
      bdDefaultVal=1500
      ${cdoQ} -setmisstoc,${bdDefaultVal} BLDFIE_3600_1800.tmp ${atmGridID}_${var}.tmp2
      ${cdoQ} -remapcon,${icon_grid} ${atmGridID}_${var}.tmp2 ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    PHIHOX )
  # ----------------------
      # soil ph: SPREAD(sb_init_vars%soil_ph) -> soil_ph_sl
      # *.tmp file provided by generate_bc_land_soil_hlp.py
      # default value
      phDefaultVal=65
      ${cdoQ} -setmisstoc,${phDefaultVal} PHIHOX.soilgrid.3600.1800.tmp ${atmGridID}_${var}.tmp2
      ${cdoQ} -remapcon,${icon_grid} ${atmGridID}_${var}.tmp2 ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    AWCtS )
  # ----------------------
      # saturated water content (volumetric fraction): not used in ICON-Land at the moment
      # *.tmp file provided by generate_bc_land_soil_hlp.py
      # default value
      awcDefaultVal=300
      ${cdoQ} -setmisstoc,${awcDefaultVal} AWCtS.soilgrid.3600.1800.tmp ${atmGridID}_${var}.tmp2
      ${cdoQ} -remapcon,"${icon_grid}" ${atmGridID}_${var}.tmp2 ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    Labile_Inorganic_P )
  # ----------------------
      # soil P labile inorganic: sb_init_vars%soil_p_labile
      # default value
      pLabileDefaultVal=50
      ${cdoQ} -setmisstoc,${pLabileDefaultVal} ${quincyInputDataDirPath}/Labile_Inorganic_P.nc ${atmGridID}_${var}.tmp
      ${cdoQ} -remapcon,${icon_grid} ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    Seconday_Mineral__P )
  # ----------------------
      # soil P slow (secondary mineral P): sb_init_vars%soil_p_slow
      # default value
      pSlowDefaultVal=70
      ${cdoQ} -setmisstoc,${pSlowDefaultVal} ${quincyInputDataDirPath}/Seconday_Mineral__P.nc ${atmGridID}_${var}.tmp
      ${cdoQ} -remapcon,${icon_grid} ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    Occluded_P )
  # ----------------------
      # soil P occluded: sb_init_vars%soil_p_occluded
      # default value
      pOcclDefaultVal=150
      ${cdoQ} -setmisstoc,${pOcclDefaultVal} ${quincyInputDataDirPath}/Occluded_P.nc ${atmGridID}_${var}.tmp
      ${cdoQ} -remapcon,${icon_grid} ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    Apatite_P )
  # ----------------------
      # soil P primary (apatite P): sb_init_vars%soil_p_primary
      # default value
      pPrimDefaultVal=100
      ${cdoQ} -setmisstoc,${pPrimDefaultVal} ${quincyInputDataDirPath}/Apatite_P.nc ${atmGridID}_${var}.tmp
      ${cdoQ} -remapcon,${icon_grid} ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    TAXOUSDA_10km )
  # ----------------------
      # USDA taxonomy soil class: sb_init_vars%usda_taxonomy_class (REAL) -> sb_init_vars%usda_taxonomy_class_int (INTEGER)
      # categorial data, not used at the moment
      # default value
      usdaTaxDefaultVal=1
      ${cdoQ} -setmisstoc,${usdaTaxDefaultVal} ${quincyInputDataDirPath}/TAXOUSDA_10km.soilgrid.3600.1800.nc ${atmGridID}_${var}.tmp
      ${cdoQ} -remaplaf,${icon_grid} ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    TAXNWRB )
  # ----------------------
      # NWRB taxonomy soil class: sb_init_vars%nwrb_taxonomy_class (REAL) -> sb_init_vars%nwrb_taxonomy_class_int (INTEGER)
      # categorial data, used to calculate 'qmax_org'
      # default value
      nwrbTaxDefaultVal=1
      ${cdoQ} -setmisstoc,${nwrbTaxDefaultVal} ${quincyInputDataDirPath}/TAXNWRB.soilgrid.3600.1800.nc ${atmGridID}_${var}.tmp
      ${cdoQ} -remaplaf,${icon_grid} ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    qmax_org )
  # ----------------------
      # qmax organic fine particle calculated below
      ;;
  # ----------------------
    * )
  # ----------------------
      echo "ERROR: Variable ${var} not known to be a bc_land_soil file variable."
      exit 1
      ;;
  esac
  } &  # comment out '&' for serial processing
done
wait


# Further process soil texture data (needs to be done outside the above parallel loop)
# --------
if [[ ${var_list} == *"FR_SAND"* || ${var_list} == *"FR_SILT"* || ${var_list} == *"FR_CLAY"* ]]; then
  if [[ -f ${atmGridID}_SUB_FR_SAND_MASK.nc \
     && -f ${atmGridID}_SUB_FR_SILT_MASK.nc \
     && -f ${atmGridID}_SUB_FR_CLAY_MASK.nc  ]]; then
    # Generate mask with grid cells where any of the soil texture data is missing
    ${cdo} -O ensmin \
        ${atmGridID}_SUB_FR_SAND_MASK.nc ${atmGridID}_SUB_FR_SILT_MASK.nc ${atmGridID}_SUB_FR_CLAY_MASK.nc \
        ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc

    # Apply this mask to the texture data
    ${cdo} mul ${atmGridID}_SUB_FR_SAND.tmp ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc ${atmGridID}_SUB_FR_SAND_0.nc
    ${cdo} mul ${atmGridID}_SUB_FR_SILT.tmp ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc ${atmGridID}_SUB_FR_SILT_0.nc
    ${cdo} mul ${atmGridID}_SUB_FR_CLAY.tmp ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc ${atmGridID}_SUB_FR_CLAY_0.nc

    # Fill data gaps with 50 percent sand, 25 percent silt and 25 percent clay fraction
    ${cdo} add ${atmGridID}_SUB_FR_SAND_0.nc -mulc,-0.5  -addc,-1. ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_SUB_FR_SAND.nc
    ${cdo} add ${atmGridID}_SUB_FR_SILT_0.nc -mulc,-0.25 -addc,-1. ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_SUB_FR_SILT.nc
    ${cdo} add ${atmGridID}_SUB_FR_CLAY_0.nc -mulc,-0.25 -addc,-1. ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_SUB_FR_CLAY.nc
  else
    echo "ERROR: Deep soil texture variables incomplete"
    exit 1
  fi
  if [[ -f ${atmGridID}_FR_SAND_MASK.nc \
     && -f ${atmGridID}_FR_SILT_MASK.nc \
     && -f ${atmGridID}_FR_CLAY_MASK.nc  ]]; then
    # Generate mask with grid cells where any of the soil textures are missing.
    ${cdo} -O ensmin \
        ${atmGridID}_FR_SAND_MASK.nc ${atmGridID}_FR_SILT_MASK.nc ${atmGridID}_FR_CLAY_MASK.nc \
        ${atmGridID}_SAND_SILT_CLAY_MASK.nc

    # Apply this mask to the texture data
    ${cdo} mul ${atmGridID}_FR_SAND.tmp ${atmGridID}_SAND_SILT_CLAY_MASK.nc ${atmGridID}_FR_SAND_0.nc
    ${cdo} mul ${atmGridID}_FR_SILT.tmp ${atmGridID}_SAND_SILT_CLAY_MASK.nc ${atmGridID}_FR_SILT_0.nc
    ${cdo} mul ${atmGridID}_FR_CLAY.tmp ${atmGridID}_SAND_SILT_CLAY_MASK.nc ${atmGridID}_FR_CLAY_0.nc

    # Fill data gaps with 50 percent sand, 25 percent silt and 25 percent clay fraction
    ${cdo} add ${atmGridID}_FR_SAND_0.nc -mulc,-0.5  -addc,-1. ${atmGridID}_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_FR_SAND.nc
    ${cdo} add ${atmGridID}_FR_SILT_0.nc -mulc,-0.25 -addc,-1. ${atmGridID}_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_FR_SILT.nc
    ${cdo} add ${atmGridID}_FR_CLAY_0.nc -mulc,-0.25 -addc,-1. ${atmGridID}_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_FR_CLAY.nc
  else
    echo "ERROR: Upper soil texture variables incomplete"
    exit 1
  fi
fi

# Maxmoist: water content of the root zone at field capacity
# --------
# Note: maxmoist is calculated here - outside the variable loop - because it depends on
#       soil_field_cap and root_depth generated above.
# Note: maxmoist is no longer read from initial files in recent ICON-Land hydrology versions
# as it depends on soil texture and organic layer fractions. It is calculated at runtime
# from field capacity and root depth. For backward compatibility we re-calculate it here to
# be consistent with root depth.

${cdo} -setattribute,maxmoist@long_name="Maximum amount of soil moisture" -setattribute,maxmoist@units="m" \
    -setname,maxmoist -mul ${atmGridID}_root_depth.nc ${atmGridID}_soil_field_cap.nc \
    ${atmGridID}_maxmoist.nc

# qmax organic fine particle
# --------
# sb_init_vars%qmax_org_fine_particle
# create based on NWRB taxonomy soil class (TAXNWRB)
qmaxVarName=qmax_org   # needs to be member of 'var_list', see below variable merging into output file
qmaxInputFile=${quincyInputDataDirPath}/TAXNWRB.soilgrid.3600.1800.nc
qmaxDefaultVal=6.5537
# set values based on category as documented in https://gitlab.dkrz.de/quincy-community/quincy-model-developers/qs-development/-/blob/develop/contrib/forcing_scripts/qmax_org_values_per_nwrb_category_20180515.csv
${cdoQ} -setvals,1,3.66601005727987,2,3.66601005727987,3,3.66601005727987,4,3.66601005727987,5,3.66601005727987,6,3.66601005727987,7,6.5537123684561,8,6.5537123684561,9,6.5537123684561,10,6.5537123684561,11,6.5537123684561,12,12.7381110963101,13,12.7381110963101,14,12.7381110963101,15,6.24417210603436,16,6.24417210603436,17,6.24417210603436,18,6.24417210603436,19,6.24417210603436,20,6.24417210603436,21,6.5537123684561,22,6.5537123684561,23,6.5537123684561,24,6.5537123684561,25,6.5537123684561,26,3.66601005727987,27,6.5537123684561,28,6.5537123684561,29,6.5537123684561,30,6.5537123684561,31,6.5537123684561,32,6.5537123684561,33,6.5537123684561,34,6.5537123684561,35,6.5537123684561,36,6.5537123684561,37,6.5537123684561,38,6.5537123684561,39,6.24417210603436,40,6.24417210603436,41,6.24417210603436,42,3.66601005727987,43,3.66601005727987,44,3.66601005727987,45,3.66601005727987,46,3.66601005727987,47,3.66601005727987,48,6.5537123684561,49,6.5537123684561,50,6.5537123684561,51,6.5537123684561,52,6.5537123684561,53,6.24417210603436,54,6.24417210603436,55,6.24417210603436,56,6.24417210603436,57,6.24417210603436,58,6.24417210603436,59,6.5537123684561,60,6.5537123684561,61,6.24417210603436,62,6.24417210603436,63,6.24417210603436,64,6.24417210603436,65,6.24417210603436,66,6.5537123684561,67,6.5537123684561,68,6.5537123684561,69,6.5537123684561,70,6.5537123684561,71,6.5537123684561,72,6.5537123684561,73,3.66601005727987,74,3.66601005727987,75,3.66601005727987,76,6.5537123684561,77,6.5537123684561,78,6.5537123684561,79,6.5537123684561,80,6.5537123684561,81,6.5537123684561,82,6.5537123684561,83,6.5537123684561,84,6.5537123684561,85,3.66601005727987,86,3.66601005727987,87,6.5537123684561,88,6.5537123684561,89,6.5537123684561,90,6.5537123684561,91,6.5537123684561,92,6.5537123684561,93,6.5537123684561,94,6.5537123684561,95,3.66601005727987,96,3.66601005727987,97,6.5537123684561,98,6.5537123684561,99,6.5537123684561,100,6.5537123684561,101,6.5537123684561,102,6.5537123684561,103,6.5537123684561,104,6.5537123684561,105,6.5537123684561,106,6.5537123684561,107,6.5537123684561,108,6.5537123684561,109,6.5537123684561,110,6.5537123684561,111,6.5537123684561,112,6.24417210603436,113,6.5537123684561,114,6.5537123684561,115,6.5537123684561,116,6.5537123684561,117,6.5537123684561,118,6.5537123684561 "${qmaxInputFile}" "Q_max_org_values.tmp"
${cdoQ} -setmisstoc,${qmaxDefaultVal} "Q_max_org_values.tmp" "Q_max_org_values.tmp2"
${cdoQ} -chname,TAXNWRB,${qmaxVarName} -remaplaf,"${icon_grid}" "Q_max_org_values.tmp2" ${atmGridID}_${qmaxVarName}.nc


# Merge newly calculated variables (*.nc files) to the bc soil file
# --------
# remove temporary file_bc_soil file (from previous runs)
[[ -f ${file_bc_soil}.tmp ]] && ${rm} ${file_bc_soil}.tmp
# loop over variable list
for var in ${var_list}; do
  if [[ ! -f ${file_bc_soil}.tmp ]]; then
    # first variable
    ${cp} ${atmGridID}_${var}.nc  ${file_bc_soil}.tmp
  else
    # following variable
    mv ${file_bc_soil}.tmp ${file_bc_soil}.tmp2
    ${cdo} -O merge ${file_bc_soil}.tmp2  ${atmGridID}_${var}.nc ${file_bc_soil}.tmp
  fi
done

${cdo} --no_history setattribute,history="${history_att}" ${file_bc_soil}.tmp \
       ${path_bc}/${file_bc_soil}.nc

# Clean up
# --------
# remove tmp files from current directory (${work_dir})
if [[ ${clean_up} == true ]]; then
  ${rm} *.tmp*
fi
echo "${prog}:     ${file_bc_soil}.nc         done"

exit 0
