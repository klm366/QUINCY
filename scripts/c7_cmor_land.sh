#!/bin/bash
#-------------------------------------------------------------------------------------------------
# CMOR-like postprocessing of jsbach output
#
#    The output files are not cmor compliant but have cmor variable names and units.
#
#       --- WARNING: not tested, incomplete and preliminary ---
#
#                                                              Veronika Gayler  MPIM, August 2025
#-------------------------------------------------------------------------------------------------
#
#  12 PFT setup
# --------------
#  pft01:  Tropical evergreen trees
#  pft02:  Tropical deciduous trees
#  pft03:  Temperate broadl. evergreen + evergreen_conifer
#  pft04:  Temperate broadl. deciduous + deciduous conifer
#  pft05:  Raingreen shrub
#  pft06:  Deciduous shrub
#  pft07:  C3-grass
#  pft08:  C4-grass
#  pft09:  C3-pasture
#  pft10:  C4-pasture
#  pft11:  C3-crop
#  pft12:  C4-crop
#
set -e

cdo="cdo -s"
#-------------------------------------------------------------------------------------------------
# CMOR variable list: MIP "CAP7" - Realm "land"
#--------------------

varlist_cap7="baresoilFrac burntFractionAll cLand cLitter cLitterCwd cProduct cSoil cSoilPools \
    cVeg cVegGrass cVegShrub cVegTree cropFrac evspsblsoi fFireAll gpp gppGrass gppShrub gppTree \
    grassFrac lai landCoverFrac mrfso mrro mrrob mrros mrso mrsofc mrsol mrsos nbp \
    netAtmosLandCO2Flux npp nppGrass nppShrub nppTree orog pastureFrac ra raGrass raShrub raTree \
    residualFrac rh rhGrass rhShrub rhTree rootd sftgif sftlaf shrubFrac slthick snc snd snw \
    tran treeFrac tsl vegFrac"
# Additional CMOR variables not requested in CAP7
varlist_cmor="${varlist_cap7} c3PftFrac c4PftFrac cnc cropFracC3 cropFracC4 depthsl grassFracC3 \
    grassFracC4 mrsfl mrsll pastureFracC3 pastureFracC4"
varlist_all="${varlist_cmor}"
#varlist_all="rh"

test_list="test01 test02 test03 test04 test05 test06 test07 test08 test09 test10 \
           test11 test12 test13 test14 test15 test16 test17 test18 test19 test20 \
           test21 test22 test23 test24        test26 test27"
#test_list="test26"
#-------------------------------------------------------------------------------------------------

data=/scratch/local1/m220053/CMIP7
datadir=${data}/testdata                      # Directory containing model output

bc_soil=${datadir}/bc_land_soil.nc
bc_sso=${datadir}/bc_land_sso.nc

# As soon as basic output is complete, the second file should no longer be needed
# TODO: Adapt script to run with fractions only included in separate file
lnd_2d=${datadir}/cmip7-production-cdriven-PiCtl-0007_lnd_P1M_mean_21510101T000000Z.nc
lnd_frac=${lnd_2d}

do_tests=yes

cmordir=${data}/cmor                    # Directory for the CMOR-like output
testdir=${data}/tests                   # Directory for consitency test arrays (with do_tests=yes)

workdir=/scratch/local1/m220053/CMIP7/work    # Temporary; will be deleted
[[ -d ${workdir} ]] || mkdir -p ${workdir}
cd ${workdir}; rm -f *

#-------------------------------------------------------------------------------------------------
# Preparations

# Land fraction rel. to box tile, i.e. 1-lake fraction
$cdo selvar,fract_land_box ${lnd_frac} fract_land_box.nc

# Land fraction rel. to grid cell
$cdo selvar,fract_box ${lnd_frac} fract_box.nc

# Grid cells with a land fraction
$cdo gtc,0 fract_land_box.nc with_land.nc

# Unit conversion factors
#------------------------
# mol(C) to kg(C)
C_mol2kg=0.012011
# kg(CO2) to kg(C) : molar Mass CO2: 44.0095;  molar Mass C: 12.01
CO2_to_C=0.272912

rhoh2o=1000.  # [kg/m3]
rhoice=916.7  # [kg/m3]

#-------------------------------------------------------------------------------------------------
# CMOR processing

for var in ${varlist_all}; do
echo "Processing variable ${var}"
case $var in
  baresoilFrac )
    # Percentage of entire grid cell that is covered by bare soil.
    $cdo -setattribute,${var}@long_name="Bare Soil Percentage Area Coverage" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*(1.-pheno_fract_fpc_max_box-fract_glac_box)*fract_land" \
               ${lnd_frac} ${var}.nc
    ;;
  burntFractionAll )
    # Percentage of grid cell burned due to all fires including natural and anthropogenic fires and
    # those associated with anthropogenic land use change.
    $cdo -setattribute,${var}@long_name="Percentage of Entire Grid Cell That Is Covered by Burnt Vegetation (All Classes)" \
        -setattribute,${var}@units="%" \
        -muldpm -expr,"${var}=100.*(disturb_burned_fract_box)*fract_box" \
        ${lnd_2d} ${var}.nc
   ;;
  c3PftFrac )
    # Percentage of entire grid cell that is covered by C3 PFTs (including grass, crops, and trees).
    $cdo -setattribute,${var}@long_name="Percentage Cover by C3 Plant Functional Type" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*(fract_pft01+fract_pft02+fract_pft03+fract_pft04+fract_pft05+fract_pft06 \
                           +fract_pft07+fract_pft09+fract_pft11)*pheno_fract_fpc_max_box" \
        ${lnd_2d} ${var}.nc
   ;;
  c4PftFrac )
    # Percentage of entire grid cell that is covered by C4 PFTs (including grass, crops, and trees).
    $cdo -setattribute,${var}@long_name="Percentage Cover by C4 Plant Functional Type" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*(fract_pft08+fract_pft10+fract_pft12)*pheno_fract_fpc_max_box" \
        ${lnd_2d} ${var}.nc
   ;;
  cLand )
    # Report missing data over ocean grid cells. For fractional land report value averaged over the
    # land fraction.
    # Note: c_sum_natural_ta also contains the (anthropogenic) crop harvest pool.
    $cdo -setattribute,${var}@long_name="Total Carbon in All Terrestrial Carbon Pools" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${C_mol2kg}*(carbon_c_sum_natural_ta_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  cLitter )
    # "Litter" is dead plant material in or above the soil. [...] The sum of the quantities with standard
    # names surface_litter_mass_content_of_carbon and subsurface_litter_mass_content_of_carbon has the
    # standard name litter_mass_content_of_carbon.
    # Note: Below ground yasso pools are counted as soil carbon. Thus we cannot distinguish between
    #       cLitterSurf and cLitterSubSurf.
    $cdo -setattribute,${var}@long_name="Carbon Mass in Litter Pool" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${C_mol2kg}*(carbon_c_sum_litter_ag_ta_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  cLitterCwd )
    # "Wood debris" means dead organic matter composed of coarse wood. It is distinct from fine litter.
    $cdo -setattribute,${var}@long_name="Carbon Mass in Coarse Woody Debris" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${C_mol2kg}*(carbon_c_acid_ag2_ta_box+carbon_c_water_ag2_ta_box
                                       +carbon_c_ethanol_ag2_ta_box+carbon_c_nonsoluble_ag2_ta_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  cnc )
    # Canopy covered fraction
    $cdo -setattribute,${var}@long_name="Canopy Covered Area Percentage" \
        -setattribute,${var}@units="%" \
        -mul -expr,"${var}=100.*(pheno_fract_fpc_box)" \
        ${lnd_2d} fract_box.nc ${var}.nc

    ;;
  cProduct )
    # Carbon mass per unit area that has been removed from the environment through land use change.
    # COMMENT: The variable only contains a crop harvest pool.
if [[ no_var_missing == true ]]; then
    $cdo -setattribute,${var}@long_name="Carbon Mass in Products of Land-Use Change" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${C_mol2kg}*carbon_c_crop_harvest_ta_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
fi
    ;;
  cSoil )
    # Carbon mass in the full depth of the soil model.
    $cdo -setattribute,${var}@long_name="Carbon Mass in Model Soil Pool" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${C_mol2kg}*(carbon_c_sum_humus_ta_box+carbon_c_sum_litter_bg_ta_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  cSoilPools )
    # For models with multiple soil carbon pools, report each pool here. If
    # models also have vertical discretisation these should be aggregated.
    $cdo -setattribute,C_acid_soluble_green@long_name="Carbon mass of acid soluble matter from green plant material in Soil" \
        -setattribute,C_acid_soluble_green@units="kg m-2" \
        -div -expr,"C_acid_soluble_green=${C_mol2kg}*(carbon_c_acid_bg1_ta_box)" \
        ${lnd_2d} fract_land_box.nc cSoilPool01.tmp
    $cdo -setattribute,C_water_soluble_green@long_name="Carbon mass of water soluble matter from green plant material in Soil" \
        -setattribute,C_water_soluble_green@units="kg m-2" \
        -div -expr,"C_water_soluble_green=${C_mol2kg}*(carbon_c_water_bg1_ta_box)" \
        ${lnd_2d} fract_land_box.nc cSoilPool02.tmp
    $cdo -setattribute,C_ethanol_soluble_green@long_name="Carbon Mass of ethanol soluble matter from green plant material in Soil" \
        -setattribute,C_ethanol_soluble_green@units="kg m-2" \
        -div -expr,"C_ethanol_soluble_green=${C_mol2kg}*(carbon_c_ethanol_bg1_ta_box)" \
        ${lnd_2d} fract_land_box.nc cSoilPool03.tmp
    $cdo -setattribute,C_non_soluble_green@long_name="Carbon Mass of non soluble matter from green plant material in Soil" \
        -setattribute,C_non_soluble_green@units="kg m-2" \
        -div -expr,"C_non_soluble_green=${C_mol2kg}*(carbon_c_nonsoluble_bg1_ta_box)" \
        ${lnd_2d} fract_land_box.nc cSoilPool04.tmp
    $cdo -setattribute,C_acid_soluble_woody@long_name="Carbon Mass of acid soluble matter from woody material in Soil" \
        -setattribute,C_acid_soluble_woody@units="kg m-2" \
        -div -expr,"C_acid_soluble_woody=${C_mol2kg}*(carbon_c_acid_bg2_ta_box)" \
        ${lnd_2d} fract_land_box.nc cSoilPool05.tmp
    $cdo -setattribute,C_water_soluble_woody@long_name="Carbon Mass of water soluble matter from woody material in Soil" \
        -setattribute,C_water_soluble_woody@units="kg m-2" \
        -div -expr,"C_water_soluble_woody=${C_mol2kg}*(carbon_c_water_bg2_ta_box)" \
        ${lnd_2d} fract_land_box.nc cSoilPool06.tmp
    $cdo -setattribute,C_ethanol_soluble_woody@long_name="Carbon Mass of ethanol soluble matter from woody material in Soil" \
        -setattribute,C_ethanol_soluble_woody@units="kg m-2" \
        -div -expr,"C_ethanol_soluble_woody=${C_mol2kg}*(carbon_c_ethanol_bg2_ta_box)" \
        ${lnd_2d} fract_land_box.nc cSoilPool07.tmp
    $cdo -setattribute,C_non_soluble_woody@long_name="Carbon Mass of non soluble matter from woody material in Soil" \
        -setattribute,C_non_soluble_woody@units="kg m-2" \
        -div -expr,"C_non_soluble_woody=${C_mol2kg}*(carbon_c_nonsoluble_bg2_ta_box)" \
        ${lnd_2d} fract_land_box.nc cSoilPool08.tmp
    $cdo -setattribute,C_humus_fine@long_name="Carbon Mass of humus from green matter in Soil" \
        -setattribute,C_humus_fine@units="kg m-2" \
        -div -expr,"C_humus_fine=${C_mol2kg}*(carbon_c_humus_1_ta_box)" \
        ${lnd_2d} fract_land_box.nc cSoilPool09.tmp
    $cdo -setattribute,C_humus_coarse@long_name="Carbon Mass of humus from woody matter in Soil" \
        -setattribute,C_humus_coarse@units="kg m-2" \
        -div -expr,"C_humus_coarse=${C_mol2kg}*(carbon_c_humus_2_ta_box)" \
        ${lnd_2d} fract_land_box.nc cSoilPool10.tmp

    $cdo merge cSoilPool??.tmp ${var}.nc
    ;;
  cVeg )
    # Carbon mass per unit area in vegetation.
    $cdo -setattribute,${var}@long_name="Carbon Mass in Vegetation" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${C_mol2kg}*(carbon_c_sum_veg_ta_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  cVegGrass )
    # mean where natural_grasses (mask=grassFrac)
    $cdo -setattribute,${var}@long_name="Carbon Mass in Vegetation on Grass Tiles" \
        -setattribute,${var}@units="kg m-2" \
        -expr,"${var}=${C_mol2kg}*(fract_pft07*carbon_c_sum_veg_ta_pft07+fract_pft08*carbon_c_sum_veg_ta_pft08) \
                                 /(fract_pft07+fract_pft08)" \
        ${lnd_2d} ${var}.nc
    ;;
  cVegShrub )
    # mean where shrubs (mask=shrubFrac)
    $cdo -setattribute,${var}@long_name="Carbon Mass in Vegetation on Grass Tiles" \
        -setattribute,${var}@units="kg m-2" \
        -expr,"${var}=${C_mol2kg}*(fract_pft05*carbon_c_sum_veg_ta_pft05+fract_pft06*carbon_c_sum_veg_ta_pft06) \
                                 /(fract_pft05+fract_pft06)" \
        ${lnd_2d} ${var}.nc
    ;;
  cVegTree )
    # mean where trees (mask=treeFrac)
    $cdo -setattribute,${var}@long_name="Carbon Mass in Vegetation on Tree Tiles" \
        -setattribute,${var}@units="kg m-2" \
        -expr,"${var}=${C_mol2kg}*(fract_pft01*carbon_c_sum_veg_ta_pft01+fract_pft02*carbon_c_sum_veg_ta_pft02 \
                                  +fract_pft03*carbon_c_sum_veg_ta_pft03+fract_pft04*carbon_c_sum_veg_ta_pft04) \
                                 /(fract_pft01+fract_pft02+fract_pft03+fract_pft04)" \
        ${lnd_2d} ${var}.nc
    ;;
  cropFrac )
    # Percentage of entire grid cell that is covered by crop.
    $cdo -setattribute,${var}@long_name="Percentage Crop Cover" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*(fract_pft11+fract_pft12)*pheno_fract_fpc_max_box" ${lnd_frac} ${var}.nc
    ;;
  cropFracC3 )
    # Percentage of entire grid cell that is covered by C3 crops.
    $cdo -setattribute,${var}@long_name="Percentage Crop Cover" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*fract_pft11*pheno_fract_fpc_max_box" ${lnd_frac} ${var}.nc
    ;;
  cropFracC4 )
    # Percentage of entire grid cell that is covered by C4 crops.
    $cdo -setattribute,${var}@long_name="Percentage Crop Cover" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*fract_pft12*pheno_fract_fpc_max_box" ${lnd_frac} ${var}.nc
    ;;
  depthsl )
    # Total (cumulative) thickness of all soil layers. This is the sum of individual thicknesses of all soil layers.
    $cdo -setattribute,${var}@long_name="Total (Cumulative) Thickness of All Soil Layers" \
        -setattribute,${var}@units="m" \
        -setvar,${var} -ifthen with_land.nc -selvar,soil_depth ${bc_soil} ${var}.nc
    ;;
  evspsblsoi )
    # includes sublimation
    $cdo -setattribute,${var}@long_name="Water Evaporation from Soil" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=-1.*hydro_evapotrans_soil_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  fFireAll )
    # From all sources, Including natural, anthropogenic and Land-use change. Only total fire
    # emissions can be compared to observations.
    $cdo -setattribute,${var}@long_name="Carbon Mass Flux into Atmosphere Due to CO2 Emission from Fire Including All Sources" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=${CO2_to_C}*carbon_co2_l2a_fire_ta_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  gpp )
    # The rate of synthesis of biomass from inorganic precursors by autotrophs ("producers")
    # expressed as the mass of carbon which it contains. For example, photosynthesis in plants
    # or phytoplankton. The producers also respire some of this biomass and the difference is
    # referred to as the net primary production.
    # [kgC m-2 s-1]  mean where land
    $cdo -setattribute,${var}@long_name="Carbon Mass Flux out of Atmosphere Due to Gross Primary Production on Land" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=${C_mol2kg}*carbon_gpp_ydaymean_ta_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  gppGrass )
    # [kgC m-2 s-1]  mean where natural_grasses (mask=grassFrac)
    $cdo -setattribute,${var}@long_name="Gross Primary Production on Grass Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=${C_mol2kg}*(fract_pft07*carbon_gpp_ydaymean_ta_pft07+fract_pft08*carbon_gpp_ydaymean_ta_pft08) \
                                 /(fract_pft07+fract_pft08)" \
        ${lnd_2d} ${var}.nc
    ;;
  gppShrub )
    # [kgC m-2 s-1]  mean where shrubs (mask=shrubFrac)
    $cdo -setattribute,${var}@long_name="Gross Primary Production on Shrub Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=${C_mol2kg}*(fract_pft05*carbon_gpp_ydaymean_ta_pft05+fract_pft06*carbon_gpp_ydaymean_ta_pft06) \
                                 /(fract_pft05+fract_pft06)" \
        ${lnd_2d} ${var}.nc
    ;;
  gppTree )
    # [kgC m-2 s-1]  mean where trees (mask=treeFrac)
    $cdo -setattribute,${var}@long_name="Gross Primary Production on Tree Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=${C_mol2kg}*(fract_pft01*carbon_gpp_ydaymean_ta_pft01+fract_pft02*carbon_gpp_ydaymean_ta_pft02  \
                                  +fract_pft03*carbon_gpp_ydaymean_ta_pft03+fract_pft04*carbon_gpp_ydaymean_ta_pft03) \
                                 /(fract_pft01+fract_pft02+fract_pft03+fract_pft04)" \
        ${lnd_2d} ${var}.nc
    ;;
  grassFrac )
    # Percentage of entire grid cell that is covered by natural grass.
    $cdo -setattribute,${var}@long_name="Natural Grass Area Percentage" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*(fract_pft07+fract_pft08)*pheno_fract_fpc_max_box" ${lnd_frac} ${var}.nc
    ;;
  grassFracC3 )
    # Percentage of entire grid cell that is covered by C3 natural grass.
    $cdo -setattribute,${var}@long_name="C3 Natural Grass Area Percentage" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*fract_pft07*pheno_fract_fpc_max_box" ${lnd_frac} ${var}.nc
    ;;
  grassFracC4 )
    # Percentage of entire grid cell that is covered by C4 natural grass.
    $cdo -setattribute,${var}@long_name="C4 Natural Grass Area Percentage" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*fract_pft08*pheno_fract_fpc_max_box" ${lnd_frac} ${var}.nc
    ;;
  lai )
    # A ratio obtained by dividing the total upper leaf surface area of vegetation by the
    # (horizontal) surface area of the land on which it grows.
    $cdo -setattribute,${var}@long_name="Leaf Area Index" \
        -setattribute,${var}@units="1" \
        -div -expr,"${var}=pheno_lai_ta_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  landCoverFrac )
    # The categories may differ from model to model, depending on their PFT definitions. This may
    # include natural PFTs, anthropogenic PFTs, bare soil, lakes, urban areas, etc. Sum of all
    # should equal the fraction of the grid-cell that is land.
    for pft in bare glac lake pft01 pft02 pft03 pft04 pft05 pft06 pft07 pft08 pft09 pft10 pft11 pft12; do
      case ${pft} in
        pft* )
          $cdo -setattribute,${pft}@units="%" \
               -expr,"${pft}=100.*fract_${pft}*pheno_fract_fpc_max_box" ${lnd_frac} ${pft}.tmp
          ;;
        bare )
          $cdo -setattribute,${pft}@units="%" \
               -expr,"${pft}=100.*(1.-pheno_fract_fpc_max_box-fract_glac_box)*fract_land" \
               ${lnd_frac} ${pft}.tmp
          ;;
        glac )
          $cdo -setattribute,${pft}@units="%" \
               -expr,"${pft}=100.*fract_${pft}" ${lnd_frac} ${pft}.tmp
          ;;
        lake )
          # For some reason some glacier grid cells have a lake fraction. This (hopefully) has
          # no physical effect, but we need to remove those lakes here.
          $cdo -setattribute,${pft}@units="%" \
               -expr,"${pft}=100.*fract_${pft}*(1.-fract_glac_box)" ${lnd_frac} ${pft}.tmp
          ;;
      esac
      case ${pft} in
        pft01 ) $cdo setattribute,${pft}@long_name="Tropical evergreen trees"  ${pft}.tmp ${pft}.nc ;;
        pft02 ) $cdo setattribute,${pft}@long_name="Tropical deciduous trees"  ${pft}.tmp ${pft}.nc ;;
        pft03 ) $cdo setattribute,${pft}@long_name="Temperate evergreen trees" ${pft}.tmp ${pft}.nc ;;
        pft04 ) $cdo setattribute,${pft}@long_name="Temperate deciduous trees" ${pft}.tmp ${pft}.nc ;;
        pft05 ) $cdo setattribute,${pft}@long_name="Raingreen shrub"           ${pft}.tmp ${pft}.nc ;;
        pft06 ) $cdo setattribute,${pft}@long_name="Deciduous shrub"           ${pft}.tmp ${pft}.nc ;;
        pft07 ) $cdo setattribute,${pft}@long_name="C3-grass"                  ${pft}.tmp ${pft}.nc ;;
        pft08 ) $cdo setattribute,${pft}@long_name="C4-grass"                  ${pft}.tmp ${pft}.nc ;;
        pft09 ) $cdo setattribute,${pft}@long_name="C3-pasture"                ${pft}.tmp ${pft}.nc ;;
        pft10 ) $cdo setattribute,${pft}@long_name="C4-pasture"                ${pft}.tmp ${pft}.nc ;;
        pft11 ) $cdo setattribute,${pft}@long_name="C3-crop"                   ${pft}.tmp ${pft}.nc ;;
        pft12 ) $cdo setattribute,${pft}@long_name="C4-crop"                   ${pft}.tmp ${pft}.nc ;;
        bare )  $cdo setattribute,${pft}@long_name="Bare soil"                 ${pft}.tmp ${pft}.nc ;;
        glac )  $cdo setattribute,${pft}@long_name="Glacier"                   ${pft}.tmp ${pft}.nc ;;
        lake )  $cdo setattribute,${pft}@long_name="Lakes"                     ${pft}.tmp ${pft}.nc ;;
      esac
      rm ${pft}.tmp
    done
    $cdo -setattribute,long_name="Percentage of Area by Vegetation or Land-Cover Category" \
        -merge bare.nc glac.nc lake.nc \
               pft01.nc pft02.nc pft03.nc pft04.nc pft05.nc pft06.nc \
               pft07.nc pft08.nc pft09.nc pft10.nc pft11.nc pft12.nc \
        ${var}.nc
    ;;
  mrfso )
    # the mass (summed over all all layers) of frozen water.
    $cdo -setattribute,${var}@long_name="Soil Frozen Water Content" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${rhoice}*vertsum(hydro_ice_soil_sl_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  mrro )
    # computed as the total runoff (including "drainage" through the base of the soil model) leaving the
    # land portion of the grid cell divided by the land area in the grid cell.
    # Note: We need hydro_runoff_land, excluding lakes and glaciers, where runoff is calculated as
    #  "runoff(glacier;lake) = rain + snow + evapotrans" just to close the water balance.
    $cdo -setattribute,${var}@long_name="Total Runoff" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=hydro_runoff_land+hydro_drainage_land" \
        ${lnd_2d} ${var}.nc
    ;;
  mrrob )
    # Not in CAP7
    $cdo -setattribute,${var}@long_name="Subsurface Runoff" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=hydro_drainage_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  mrros )
    # the total surface runoff leaving the land portion of the grid cell.
    $cdo -setattribute,${var}@long_name="Surface Runoff" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=hydro_runoff_land" \
        ${lnd_2d} ${var}.nc
    ;;
  mrsfl )
    # in each soil layer, the mass of water in ice phase. Reported as "missing" for grid cells occupied entirely by "sea"
    $cdo -setattribute,${var}@long_name="Frozen Water Content of Soil Layer" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${rhoice}*hydro_ice_soil_sl_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  mrsll )
    # in each soil layer, the mass of water in liquid phase. Reported as "missing" for grid cells occupied entirely by "sea"
    $cdo -setattribute,${var}@long_name="Liquid Water Content of Soil Layer" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${rhoh2o}*hydro_wtr_soil_sl_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  mrso )
    # the mass per unit area (summed over all soil layers) of water in all phases.
    $cdo -setattribute,${var}@long_name="Total Soil Moisture Content" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=vertsum(${rhoh2o}*hydro_wtr_soil_sl_box+${rhoice}*hydro_ice_soil_sl_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  mrsofc )
    # reported "where land": divide the total water holding capacity of all the soil in the grid
    # cell by the land area in the grid cell; reported as "missing" where the land fraction is 0.
    $cdo -setattribute,${var}@long_name="Capacity of Soil to Store Water (Field Capacity)" \
        -setattribute,${var}@units="kg m-2" \
        -ifthen with_land.nc -expr,"${var}=${rhoh2o}*(soil_depth*soil_field_cap)" \
        ${bc_soil} ${var}.nc
    ;;
  mrsol )
    # in each soil layer, the mass of water in all phases, including ice. Reported as "missing"
    # for grid cells occupied entirely by "sea"
    $cdo -setattribute,${var}@long_name="Total Water Content of Soil Layer" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=(${rhoh2o}*hydro_wtr_soil_sl_box+${rhoice}*hydro_ice_soil_sl_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  mrsos )
    # the mass of water in all phases in a thin surface soil layer
    # Note: uppermost 10 cm: layer 1 (6.5cm) + 3.5cm of layer 2
    $cdo -setattribute,${var}@long_name="Moisture in Upper Portion of Soil Column" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${rhoh2o}*(sellevidx(hydro_wtr_soil_sl_box,1)+0.1378*sellevidx(hydro_wtr_soil_sl_box,2))  \
                          +${rhoice}*(sellevidx(hydro_ice_soil_sl_box,1)+0.1378*sellevidx(hydro_ice_soil_sl_box,2))" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  nbp )
    # This is the net mass flux of carbon between land and atmosphere calculated as photosynthesis
    # MINUS the sum of plant and soil respiration, carbon fluxes from fire, harvest, grazing and
    # land use change. Positive flux is into the land.
    $cdo -setattribute,${var}@long_name="Carbon Mass Flux out of Atmosphere Due to Net Biospheric Production on Land" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=-${CO2_to_C}*(carbon_co2_l2a_npp_ta_box+carbon_co2_l2a_resp_ta_box \
                                       +carbon_co2_l2a_herb_ta_box+carbon_co2_l2a_fire_ta_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  nep )
    # Net Ecosystem Exchange
    # Standard name:
    # surface_net_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_all_land_processes_excluding_anthropogenic_land_use_change
    # ATTENTION: carbon_co2_l2a_resp_ta_box includes a crop harvest flux. It cannot be removed here,
    #   as the corresponding jsbach variable "Cflx_crop_harvest_2_atmos" is only defined locally and
    #   not available in the output. Nep thus removed from list of available CMOR variables.
    $cdo -setattribute,${var}@long_name="Net Carbon Mass Flux out of Atmosphere Due to Net Ecosystem Productivity on Land" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=-${CO2_to_C}*(carbon_co2_l2a_npp_ta_box+carbon_co2_l2a_resp_ta_box \
                                        +carbon_co2_l2a_herb_ta_box+carbon_co2_l2a_fire_ta_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  netAtmosLandCO2Flux )
    # Flux of carbon as carbon dioxide into the land. This flux should be reproducible by
    # differencing the sum of all carbon pools (cVeg, cLitter, cSoil, and cProducts or
    # equivalently cLand) from one time step to the next, except in the case of lateral transfer
    # of carbon due to harvest, riverine transport of dissolved organic and/or inorganic carbon, or
    # any other process (in which case the lateral_carbon_transfer_over_land term, see below, will
    # be zero data).
    # Same as nbp.
    $cdo -setattribute,${var}@long_name="Net Flux of CO2 Between Atmosphere and Land (Positive into Land) as a Result of All Processes [kgC m-2 s-1]" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=-${CO2_to_C}*(carbon_co2_l2a_npp_ta_box+carbon_co2_l2a_resp_ta_box \
                                        +carbon_co2_l2a_herb_ta_box+carbon_co2_l2a_fire_ta_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  npp )
    # "Production of carbon" means the production of biomass expressed as the mass of carbon which
    # it contains. Net primary production is the excess of gross primary production (rate of
    # synthesis of biomass from inorganic precursors) by autotrophs ("producers"), for example,
    # photosynthesis in plants or phytoplankton, over the rate at which the autotrophs themselves
    # respire some of this biomass. "Productivity" means production per unit area.
    $cdo -setattribute,${var}@long_name="Net Primary Production on Land as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=-${CO2_to_C}*(carbon_co2_l2a_npp_ta_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  nppGrass )
    # Total NPP of grass in the grid cell
    $cdo -setattribute,${var}@long_name="Net Primary Production on Grass Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=-${CO2_to_C}*(fract_pft07*carbon_co2_l2a_npp_ta_pft07+fract_pft08*carbon_co2_l2a_npp_ta_pft08) \
                                  /(fract_pft07+fract_pft08)" \
        ${lnd_2d} ${var}.nc
    ;;
  nppShrub )
    # Total NPP of shrub in the grid cell
    $cdo -setattribute,${var}@long_name="Net Primary Production on Shrub Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=-${CO2_to_C}*(fract_pft05*carbon_co2_l2a_npp_ta_pft05+fract_pft06*carbon_co2_l2a_npp_ta_pft06) \
                                  /(fract_pft05+fract_pft06)" \
        ${lnd_2d} ${var}.nc
    ;;
  nppTree )
    # Total NPP of tree in the grid cell
    $cdo -setattribute,${var}@long_name="Net Primary Production on Tree Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=-${CO2_to_C}*(fract_pft01*carbon_co2_l2a_npp_ta_pft01+fract_pft02*carbon_co2_l2a_npp_ta_pft02  \
                                   +fract_pft03*carbon_co2_l2a_npp_ta_pft03+fract_pft04*carbon_co2_l2a_npp_ta_pft04) \
                                  /(fract_pft01+fract_pft02+fract_pft03+fract_pft04)" \
        ${lnd_2d} ${var}.nc
    ;;
  orog )
    # height above the geoid; as defined here, "the geoid" is a surface of constant geopotential
    # that, if the ocean were at rest, would coincide with mean sea level. Under this definition,
    # the geoid changes as the mean volume of the ocean changes (e.g., due to glacial melt, or
    # global warming of the ocean). Reported here is the height above the present-day geoid
    # (0.0 over ocean).
    $cdo -setattribute,${var}@long_name="Surface Altitude" \
        -setattribute,${var}@units="m" \
        -expr,"${var}=elevation" ${bc_sso} ${var}.nc
    ;;
  pastureFrac )
    # fraction of entire grid cell that is covered by anthropogenic pasture.
    $cdo -setattribute,${var}@long_name="Percentage of Land Which Is Anthropogenic Pasture" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*(fract_pft09+fract_pft10)*pheno_fract_fpc_max_box" ${lnd_frac} ${var}.nc
    ;;
  pastureFracC3 )
    # Percentage of entire grid cell covered by C3 pasture
    $cdo -setattribute,${var}@long_name="C3 Pasture Area Percentage" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*fract_pft09*pheno_fract_fpc_max_box" ${lnd_frac} ${var}.nc
    ;;
  pastureFracC4 )
    # Percentage of entire grid cell covered by C4 pasture
    $cdo -setattribute,${var}@long_name="C4 Pasture Area Percentage" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*fract_pft10*pheno_fract_fpc_max_box" ${lnd_frac} ${var}.nc
    ;;
  ra )
    # Carbon mass flux per unit area into atmosphere due to autotrophic respiration on land
    # (respiration by producers) [see rh for heterotrophic production]
    $cdo -setattribute,${var}@long_name="Carbon Mass Flux into Atmosphere Due to Autotrophic (Plant) Respiration on Land" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=${C_mol2kg}*assimi_dark_respiration_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  raGrass )
    # Total RA of grass in the grid cell
    $cdo -setattribute,${var}@long_name="Autotrophic Respiration on Grass Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=${C_mol2kg}*(fract_pft07*assimi_dark_respiration_pft07+fract_pft08*assimi_dark_respiration_pft08) \
                                 /(pheno_fract_fpc_max_box*(fract_pft07+fract_pft08))" \
        ${lnd_2d} ${var}.nc
    ;;
  raShrub )
    # Total RA of shrubs in the grid cell
    $cdo -setattribute,${var}@long_name="Autotrophic Respiration on Shrub Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=${C_mol2kg}*(fract_pft05*assimi_dark_respiration_pft05+fract_pft06*assimi_dark_respiration_pft06) \
                                 /(pheno_fract_fpc_max_box*(fract_pft05+fract_pft06))" \
        ${lnd_2d} ${var}.nc
    ;;
  raTree )
    # Carbon mass flux per unit area into atmosphere due to autotrophic respiration on land
    # (respiration by producers) [see rh for heterotrophic production]
    $cdo -setattribute,${var}@long_name="Carbon Mass Flux into Atmosphere Due to Autotrophic (Plant) Respiration on Land" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=${C_mol2kg}*(fract_pft01*assimi_dark_respiration_pft01+fract_pft02*assimi_dark_respiration_pft02  \
                                  +fract_pft03*assimi_dark_respiration_pft03+fract_pft04*assimi_dark_respiration_pft04) \
                                 /(pheno_fract_fpc_max_box*(fract_pft01+fract_pft02+fract_pft03+fract_pft04))" \
        ${lnd_2d} ${var}.nc
    ;;
  residualFrac )
    # fraction of entire grid cell that is land and is covered by "non-vegetation" and "non-bare-soil" (e.g., urban, ice, lakes, etc.)
    # ATTENTION: The examples are not inline with our definitions!! We exclude lakes from the land!!
    $cdo -setattribute,${var}@long_name="Percentage of Grid Cell That Is Land but neither Vegetation Covered nor Bare Soil" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*fract_glac" \
        ${lnd_frac} ${var}.nc
    ;;
  rh )
    # Carbon mass flux per unit area into atmosphere due to heterotrophic respiration on land
    # (respiration by consumers)
    $cdo -setattribute,${var}@long_name="Total Heterotrophic Respiration on Land as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=${CO2_to_C}*(carbon_co2_l2a_resp_ta_box + carbon_co2_l2a_herb_ta_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  rhGrass )
    # Total RH of grass in the grid cell
if [[ no_var_missing == true ]]; then
    $cdo -setattribute,${var}@long_name="Heterotrophic Respiration on Grass Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=${CO2_to_C}*( fract_pft07*(carbon_co2_l2a_resp_ta_pft07+carbon_co2_l2a_herb_ta_pft07) \
                                   +fract_pft08*(carbon_co2_l2a_resp_ta_pft08+carbon_co2_l2a_herb_ta_pft08)) \
                                 /(pheno_fract_fpc_max_box*(fract_pft07+fract_pft08))" \
        ${lnd_2d} ${var}.nc
fi
    ;;
  rhShrub )
    # Total RH of Shrubs in the grid cell
if [[ no_var_missing == true ]]; then
    $cdo -setattribute,${var}@long_name="Heterotrophic Respiration on Shrub Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=${CO2_to_C}*( fract_pft05*(carbon_co2_l2a_resp_ta_pft05+carbon_co2_l2a_herb_ta_pft05) \
                                   +fract_pft06*(carbon_co2_l2a_resp_ta_pft06+carbon_co2_l2a_herb_ta_pft06)) \
                                 /(pheno_fract_fpc_max_box*(fract_pft05+fract_pft06))" \
        ${lnd_2d} ${var}.nc
fi
    ;;
  rhTree )
    # Total RH of trees in the grid cell
if [[ no_var_missing == true ]]; then
    $cdo -setattribute,${var}@long_name="Heterotrophic Respiration on Tree Tiles as Carbon Mass Flux" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -expr,"${var}=${CO2_to_C}*( fract_pft01*(carbon_co2_l2a_resp_ta_pft01+carbon_co2_l2a_herb_ta_pft01) \
                                   +fract_pft02*(carbon_co2_l2a_resp_ta_pft02+carbon_co2_l2a_herb_ta_pft02)) \
                                   +fract_pft03*(carbon_co2_l2a_resp_ta_pft03+carbon_co2_l2a_herb_ta_pft03)) \
                                   +fract_pft04*(carbon_co2_l2a_resp_ta_pft04+carbon_co2_l2a_herb_ta_pft04)) \
                                 /(pheno_fract_fpc_max_box*(fract_pft01+fract_pft02+fract_pft03+fract_pft04))" \
        ${lnd_2d} ${var}.nc
fi
    ;;
  rootd )
    # report the maximum soil depth reachable by plant roots (if defined in model), i.e., the
    # maximum soil depth from which they can extract moisture; report as "missing" where the land
    # fraction is 0.
    $cdo -setattribute,${var}@long_name="Maximum Root Depth" \
        -setattribute,${var}@units="m" \
        -setvar,${var} -ifthen with_land.nc -selvar,root_depth ${bc_soil} ${var}.nc
    ;;
  sftgif )
    # fraction of grid cell occupied by "permanent" ice (i.e., glaciers).
    $cdo -setattribute,${var}@long_name="Land Ice Area Percentage" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*fract_glac" ${lnd_frac} ${var}.nc
    ;;
  sftlaf )
    # Fraction of horizontal land grid cell area occupied by lake.
    # Note: due to an inconsistency of current initial files we need to exclude lakes from glacier cells.
    $cdo -setattribute,${var}@long_name="Fraction of the Grid Cell Occupied by Lake" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*fract_lake*(1.-fract_glac_box)" ${lnd_frac} ${var}.nc
    ;;
  shrubFrac )
    # Percentage of entire grid cell that is covered by shrub.
    $cdo -setattribute,${var}@long_name="Percentage Cover by Shrub" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*(fract_pft05+fract_pft06)*pheno_fract_fpc_max_box" ${lnd_2d} ${var}.nc
    ;;
  slthick )
    # Thickness of Soil Layers
    $cdo -setattribute,${var}@long_name="Thickness of Soil Layers" \
        -setattribute,${var}@units="m" \
        -expr,"${var}=hydro_soil_depth_sl_land" ${lnd_2d} ${var}.nc
    ;;
  snc )
    # Fraction of each grid cell that is occupied by snow that rests on land portion of cell.
    $cdo -setattribute,${var}@long_name="Snow Area Percentage" \
        -setattribute,${var}@units="%" \
        -mul -expr,"${var}=100.*hydro_fract_snow_soil_box" \
        ${lnd_2d} fract_box.nc ${var}.nc
    ;;
  snd )
    # where land over land, this is computed as the mean thickness of snow in the land portion of
    # the grid cell (averaging over the entire land portion, including the snow-free fraction).
    # Reported as missing where the land fraction is 0.
    $cdo -setattribute,${var}@long_name="Snow Depth" \
        -setattribute,${var}@units="m" \
        -div -expr,"${var}=vertsum(sse_snow_depth_sl_box)" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  snw )
    # Computed as the mass of surface snow on the land portion of the grid cell divided by the land
    # area in the grid cell; reported as missing where the land fraction is 0; excluded is snow on
    # vegetation canopy or on sea ice.
    $cdo -setattribute,${var}@long_name="Surface Snow Amount" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${rhoh2o}*hydro_weq_snow_soil_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  snwc )
    # canopy_snow_amount
    $cdo -setattribute,${var}@long_name="Snow Water Equivalent Intercepted by the Vegetation" \
        -setattribute,${var}@units="kg m-2" \
        -div -expr,"${var}=${rhoh2o}*hydro_weq_snow_soil_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  tran )
    # Transpiration (may include dew formation as a negative flux).
    $cdo -setattribute,${var}@long_name="Transpiration" \
        -setattribute,${var}@units="kg m-2 s-1" \
        -div -expr,"${var}=hydro_transpiration_box" \
        ${lnd_2d} fract_land_box.nc ${var}.nc
    ;;
  treeFrac )
    # fraction of entire grid cell that is covered by trees.
    $cdo -setattribute,${var}@long_name="Tree Cover Percentage" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*(fract_pft01+fract_pft02+fract_pft03+fract_pft04)*pheno_fract_fpc_max_box" \
        ${lnd_2d} ${var}.nc
    ;;
  tsl )
    # Temperature of each soil layer. Reported as "missing" for grid cells occupied entirely by "sea".
    $cdo -setattribute,${var}@long_name="Temperature of Soil" \
        -setattribute,${var}@units="K" \
        -ifthen with_land.nc -expr,"${var}=sse_t_soil_sl_land" \
        ${lnd_2d} ${var}.nc
    ;;
  vegFrac )
    # Percentage of grid cell that is covered by vegetation.This SHOULD be the sum of tree, grass
    # (natural and pasture), crop and shrub fractions.
    $cdo -setattribute,${var}@long_name="Total Vegetated Percentage Cover" \
        -setattribute,${var}@units="%" \
        -expr,"${var}=100.*pheno_fract_fpc_max_box*fract_land" ${lnd_2d} ${var}.nc
    ;;
  * )
    echo ERROR: ${var} not supported.
    exit 1
    ;;
esac
done

echo ""
echo "$0: Variable processing done."
echo ""

# Move to CMOR output directory
[[ -d ${cmordir} ]] || mkdir -p ${cmordir}
for var in ${varlist_all}; do
  if [[ -f ${var}.nc ]]; then
    mv ${var}.nc ${cmordir}
  else
    echo "WARNING: ${var}.nc missing"
  fi
done

#------------------------------------------------------------------------------
# Consistency tests
#    All test arrays should be close to zero.
#------------------------------------------------------------------------------
if [[ ${do_tests} == yes ]]; then

  echo ""
  echo "$0: Consistency tests ..."
  echo ""

  for test in ${test_list}; do
    echo "Running ${test}"

    case $test in
      test01 ) # cropFrac = cropFracC3 + cropFracC4
        if [[ -f ${cmordir}/cropFrac.nc && ${cmordir}/cropFracC3.nc && ${cmordir}/cropFracC4.nc ]]; then
          $cdo add ${cmordir}/cropFracC3.nc \
                   ${cmordir}/cropFracC4.nc  cropFrac.test
          $cdo sub ${cmordir}/cropFrac.nc    cropFrac.test $test.nc
        fi
        ;;
      test02 ) # cSoil =  SUM(csoilPools(:))
        if [[ -f ${cmordir}/cSoil.nc && -f ${cmordir}/cSoilPools.nc ]]; then
          $cdo -expr,"cSoil=C_acid_soluble_green+C_water_soluble_green+C_ethanol_soluble_green+C_non_soluble_green \
                            +C_acid_soluble_woody+C_water_soluble_woody+C_ethanol_soluble_woody+C_non_soluble_woody \
                            +C_humus_fine+C_humus_coarse" ${cmordir}/cSoilPools.nc cSoil.test
          $cdo sub ${cmordir}/cSoil.nc cSoil.test $test.nc
        fi
        ;;
      test03 ) # cLand = cSoil + cLitter + cVeg + cProduct
        # This test does not work, as a crop harvest pool is included in "carbon_c_sum_natural_ta_box" (i.e. in cLand)
        if [[ -f ${cmordir}/cLand.nc && -f ${cmordir}/cSoil.nc && -f ${cmordir}/cLitter.nc \
           && -f ${cmordir}/cVeg.nc && -f ${cmordir}/cProduct.nc ]]; then
          $cdo -O enssum ${cmordir}/cSoil.nc    \
                         ${cmordir}/cLitter.nc  \
                         ${cmordir}/cVeg.nc     \
                         ${cmordir}/cProduct.nc  cLand.test
          $cdo sub ${cmordir}/cLand.nc cLand.test $test.nc
        fi
        ;;
      test04 ) # gpp = npp + ra
        # Leads to small mainly positive values as we are only using the ydaymean variable for GPP,
        # while ra and NPP are actuel fluxes, with the latter being corrected with ydaymean NPP_pot.
        if [[ -f ${cmordir}/npp.nc && -f ${cmordir}/gpp.nc && -f ${cmordir}/ra.nc ]]; then
          $cdo add ${cmordir}/npp.nc \
                   ${cmordir}/ra.nc   gpp.test
          $cdo sub ${cmordir}/gpp.nc  gpp.test $test.nc
        fi
        ;;
      test05 ) # grassFrac = grassFracC3 + grassFracC4
        if [[ -f ${cmordir}/grassFrac.nc && -f ${cmordir}/grassFracC3.nc && -f ${cmordir}/grassFracC4.nc ]]; then
          $cdo add ${cmordir}/grassFracC3.nc \
                   ${cmordir}/grassFracC4.nc  grassFrac.test
          $cdo sub ${cmordir}/grassFrac.nc    grassFrac.test $test.nc
        fi
        ;;
      test06 ) # pastureFrac = pastureFracC3 + pastureFracC4
        if [[ -f ${cmordir}/pastureFrac.nc && -f ${cmordir}/pastureFracC3.nc && -f ${cmordir}/pastureFracC4.nc ]]; then
          $cdo add ${cmordir}/pastureFracC3.nc \
                   ${cmordir}/pastureFracC4.nc  pastureFrac.test
          $cdo sub ${cmordir}/pastureFrac.nc    pastureFrac.test $test.nc
        fi
        ;;
      test07 ) # landCoverFrac: sum of all cover types == 100 %
        if [[ -f ${cmordir}/landCoverFrac.nc ]]; then
           $cdo -expr,"landCoverFrac= \
                       pft01+pft02+pft03+pft04+pft05+pft06+pft07+pft08+pft09+pft10+pft11+pft12 \
                      +bare+lake+glac" ${cmordir}/landCoverFrac.nc landCoverFrac.test
           $cdo -expr,"landCoverFrac=100.*fract_box" ${lnd_2d} fract_box.test
           $cdo sub landCoverFrac.test fract_box.test $test.nc
        fi
        ;;
      test08 ) # fract_land = vegFrac + baresoilFrac + residualFrac
        if [[ -f ${cmordir}/vegFrac.nc && -f ${cmordir}/baresoilFrac.nc  && -f ${cmordir}/residualFrac.nc ]]; then
          # Note: In current initial files a few glacier cells have a lake fraction! As there are no fractional glaciers,
          # these fractions need to be counted as land.
          $cdo expr,"fract_land=100.*fract_land+100.*fract_land*(fract_glac_box*fract_lake_box)" ${lnd_frac} fract_land.nc
          $cdo -add ${cmordir}/residualFrac.nc \
               -add ${cmordir}/vegFrac.nc      \
                    ${cmordir}/baresoilFrac.nc \
               fract_land.test
          $cdo sub fract_land.nc fract_land.test $test.nc
        fi
        ;;
      test09 ) # vegFrac = treeFrac + shrubFrac + grassFrac + cropFrac + pastureFrac
        if [[ -f ${cmordir}/vegFrac.nc && -f ${cmordir}/treeFrac.nc && -f ${cmordir}/shrubFrac.nc && -f ${cmordir}/grassFrac.nc \
             && -f ${cmordir}/cropFrac.nc && -f ${cmordir}/pastureFrac.nc ]]; then
          $cdo -O enssum ${cmordir}/treeFrac.nc  \
                         ${cmordir}/shrubFrac.nc  \
                         ${cmordir}/grassFrac.nc  \
                         ${cmordir}/cropFrac.nc   \
                         ${cmordir}/pastureFrac.nc  vegFrac.test
          $cdo sub ${cmordir}/vegFrac.nc            vegFrac.test $test.nc
        fi
        ;;
      test10 ) # vegfrac = c3PftFrac + c4PftFrac
        if [[ -f ${cmordir}/vegFrac.nc && -f ${cmordir}/c3PftFrac.nc && -f ${cmordir}/c4PftFrac.nc ]]; then
          $cdo add ${cmordir}/c3PftFrac.nc \
                   ${cmordir}/c4PftFrac.nc  vegFrac.test
          $cdo sub ${cmordir}/vegFrac.nc    vegFrac.test $test.nc
        fi
        ;;
      test11 ) # landCoverFrac:pft07 = grassFracC3
        if [[ -f ${cmordir}/landCoverFrac.nc && ${cmordir}/grassFracC3.nc ]]; then
          $cdo sub ${cmordir}/grassFracC3.nc \
                   -selvar,pft07 ${cmordir}/landCoverFrac.nc $test.nc
        fi
        ;;
      test12 ) # landCoverFrac:pft08 = grassFracC4
        if [[ -f ${cmordir}/landCoverFrac.nc && ${cmordir}/grassFracC4.nc ]]; then
          $cdo sub ${cmordir}/grassFracC4.nc \
                   -selvar,pft08 ${cmordir}/landCoverFrac.nc $test.nc
        fi
        ;;
      test13 ) # landCoverFrac:pft09 = pastureFracC3
        if [[ -f ${cmordir}/landCoverFrac.nc && ${cmordir}/pastureFracC3.nc ]]; then
          $cdo sub ${cmordir}/pastureFracC3.nc \
                   -selvar,pft09 ${cmordir}/landCoverFrac.nc $test.nc
        fi
        ;;
      test14 ) # landCoverFrac:pft10 = pasturecropFracC4
        if [[ -f ${cmordir}/landCoverFrac.nc && ${cmordir}/pastureFracC4.nc ]]; then
          $cdo sub ${cmordir}/pastureFracC4.nc \
                   -selvar,pft10 ${cmordir}/landCoverFrac.nc $test.nc
        fi
        ;;
      test15 ) # landCoverFrac:pft11 = cropFracC3
        if [[ -f ${cmordir}/landCoverFrac.nc && ${cmordir}/cropFracC3.nc ]]; then
          $cdo sub ${cmordir}/cropFracC3.nc \
                   -selvar,pft11 ${cmordir}/landCoverFrac.nc $test.nc
        fi
        ;;
      test16 ) # landCoverFrac:pft12 = cropFracC4
        if [[ -f ${cmordir}/landCoverFrac.nc && ${cmordir}/cropFracC4.nc ]]; then
          $cdo sub ${cmordir}/cropFracC4.nc \
                   -selvar,pft12 ${cmordir}/landCoverFrac.nc $test.nc
        fi
        ;;
      test17 ) # landCoverFrac:bare = baresoilFrac
        if [[ -f ${cmordir}/landCoverFrac.nc && ${cmordir}/baresoilFrac.nc ]]; then
          $cdo sub ${cmordir}/baresoilFrac.nc \
                   -selvar,bare ${cmordir}/landCoverFrac.nc $test.nc
        fi
        ;;
      test18 ) # landCoverFrac:glac = sftgif
        if [[ -f ${cmordir}/landCoverFrac.nc && ${cmordir}/sftgif.nc ]]; then
          $cdo sub ${cmordir}/sftgif.nc \
                   -selvar,glac ${cmordir}/landCoverFrac.nc $test.nc
        fi
        ;;
      test19 ) # landCoverFrac:lake = sftlaf
        if [[ -f ${cmordir}/landCoverFrac.nc && ${cmordir}/sftlaf.nc ]]; then
          $cdo sub ${cmordir}/sftlaf.nc \
                   -selvar,lake ${cmordir}/landCoverFrac.nc $test.nc
        fi
        ;;
      test20 ) # mrro = mrros + mrrob
        if [[ -f ${cmordir}/mrro.nc && -f ${cmordir}/mrros.nc && -f ${cmordir}/mrrob.nc ]]; then
          $cdo add ${cmordir}/mrros.nc ${cmordir}/mrrob.nc mrro.test
          $cdo sub ${cmordir}/mrro.nc mrro.test $test.nc
        fi
        ;;
      test21 ) # mrso = mrsll + mrsfl
        if [[ -f ${cmordir}/mrso.nc && -f ${cmordir}/mrsll.nc && -f ${cmordir}/mrsfl.nc ]]; then
          $cdo -vertsum -add ${cmordir}/mrsll.nc ${cmordir}/mrsfl.nc mrso.test
          $cdo sub ${cmordir}/mrso.nc mrso.test $test.nc
        fi
        ;;
      test22 ) # mrso = mrsol
        if [[ -f ${cmordir}/mrso.nc && -f ${cmordir}/mrsol.nc ]]; then
          $cdo -vertsum ${cmordir}/mrsol.nc mrso.test
          $cdo sub ${cmordir}/mrso.nc mrso.test $test.nc
        fi
        ;;
      test23 ) # mrfso = mrsfl
        if [[ -f ${cmordir}/mrfso.nc && -f ${cmordir}/mrsfl.nc ]]; then
          $cdo -vertsum ${cmordir}/mrsfl.nc mrfso.test
          $cdo sub ${cmordir}/mrfso.nc mrfso.test $test.nc
        fi
        ;;
      test24 ) # nbp = gpp - ra - rh - fFireAll (- harvest - grazing - landuse)
        # Note: grazing and crop harvest fluxes are included in rh, no flux due to land use change.
        if [[ -f ${cmordir}/nbp.nc && -f ${cmordir}/gpp.nc && -f ${cmordir}/ra.nc \
           && -f ${cmordir}/rh.nc && -f ${cmordir}/fFireAll.nc ]]; then
          $cdo -O enssum        ${cmordir}/gpp.nc      \
                      -mulc,-1. ${cmordir}/ra.nc       \
                      -mulc,-1. ${cmordir}/rh.nc       \
                      -mulc,-1. ${cmordir}/fFireAll.nc  nbp.test
          $cdo sub ${cmordir}/nbp.nc nbp.test $test.nc
        fi
        ;;
      test25 ) # nep = gpp - ra - (rh - crop_harvest) - fFireAll - grazing
        # Note: grazing and crop harvest fluxes are included i/cVegn rh. As nep excludes fluxes from
        # anthropogenic land use we would need to add crop harvest here, the variable is however
        # not available as model output.
        if [[ -f ${cmordir}/nep.nc && -f ${cmordir}/gpp.nc && -f ${cmordir}/ra.nc \
           && -f ${cmordir}/rh.nc && -f ${cmordir}/fFireAll.nc ]]; then
          $cdo -O enssum        ${cmordir}/gpp.nc      \
                      -mulc,-1. ${cmordir}/ra.nc       \
                      -mulc,-1. ${cmordir}/rh.nc       \
                      -mulc,-1. ${cmordir}/fFireAll.nc  nep.test
           $cdo sub ${cmordir}/nep.nc nep.test $test.nc
        fi
        ;;
      test26 ) # d/dt(cLand) = npp - fFireAll - rh - grazing - fCLandToOcean - fAnthDisturb - fProductDecomp
        # Note: grazing and crop harvest fluxes are included in rh, fCLandTocean=0., fAnthDisturb=0.
        # Note: Differences due to 15 day time shift between monthly mean fluxes and monthly mean state
        if [[ -f ${cmordir}/cLand.nc && -f ${cmordir}/npp.nc && -f ${cmordir}/fFireAll.nc && -f ${cmordir}/rh.nc ]]; then
          $cdo -O enssum        ${cmordir}/npp.nc      \
                      -mulc,-1. ${cmordir}/fFireAll.nc \
                      -mulc,-1. ${cmordir}/rh.nc        fluxes_cLand.test
          $cdo sub -deltat ${cmordir}/cLand.nc  -mulc,86400. -muldpm fluxes_cLand.test  $test.nc
        fi
        ;;
      test27 ) # depthsl = vertsum(slthick)
        if [[ -f ${cmordir}/depthsl.nc && -f ${cmordir}/slthick.nc ]]; then
          $cdo sub ${cmordir}/depthsl.nc -vertsum ${cmordir}/slthick.nc $test.nc
        fi
        ;;
    esac
  done
  echo ""
  echo "$0: Test array generation done."
  echo ""

  # Move test arrays to tests directory
  [[ -d ${testdir} ]] || mkdir -p ${testdir}
  for test in ${test_list}; do
    if [[ -f ${test}.nc ]]; then
      mv ${test}.nc ${testdir}
    else
      echo "WARNING: ${test}.nc missing"
    fi
  done

  echo ""
  echo "$0: Test output moved to ${testdir}"
fi

# Clean up
cd ..
rm -rf ${workdir}

echo ""
echo "$0: CMOR-like output moved to ${cmordir}"
echo "$0:   - done -"
exit 0

