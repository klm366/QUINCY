#!/bin/bash 


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


#-----------------------------------------------------------------------------
# The script reads 6-hourly CRU-JRA v3.5 forcing over 1901- 2024 (longwave radiation,
# shortwave radiation, tmax, tmin, precipitation, wind speed, specific humidity),
# computes daily statistics, inserts Feb 29 for leap years by duplicating Feb 28,
# reorders longitude from [-180,180] to [0,360], merges all variables into a
# single file, and remap conservatively onto the ICON grid.
# Optionally patches the land fraction file to mask ocean cells that have no!!!!!!!!!!!!!!!!!!!!!!
# valid climate data after remapping.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#----------------------------------------------------------------------------
#
###############################################################################
### Batch Queuing System is SLURM
#SBATCH --job-name=crujra
#SBATCH --output=crujra.o%j
#SBATCH --error=crujra.o%j
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=5120 
#SBATCH --ntasks=1
#SBATCH --account mj0143
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --time=07:30:00



module load cdo nco

DataPath="/work/mj0143/icon_land_data_pool/preprocessing/climate_forcing/crujra/rev01"
grid_def_r2b4=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc
grid_def_r2b5=/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc

outdir="./crujra_climate"
tmpdir="./tmp_crujra_climate"

mkdir -p "$outdir" "$tmpdir"

# Process CRUJRA data 1901-2024
for year in $(seq 1901 2024); do
  echo "Processing year $year"

  if (( year % 4 == 0 )); then leap=1; else leap=0; fi

  # 1. LONGWAVE
  cdo daymean \
    "${DataPath}/dlwrf/crujra.v3.5d.dlwrf.${year}.365d.noc.nc" \
    "${tmpdir}/longwave_365_${year}.nc"
  echo "done longwave"

  # 2. SHORTWAVE
  cdo daymean \
    "${DataPath}/Radiation/tswrf_gcb2025_${year}.nc" \
    "${tmpdir}/shortwave_365_${year}.nc"
  echo "done shortwave"

  # 3. TMAX
  cdo subc,273.15 \
    -daymax "${DataPath}/tmax/crujra.v3.5d.tmax.${year}.365d.noc.nc" \
    "${tmpdir}/tmax_365_${year}.nc"
  echo "done tmax"

  # 4. TMIN
  cdo subc,273.15 \
    -daymin "${DataPath}/tmin/crujra.v3.5d.tmin.${year}.365d.noc.nc" \
    "${tmpdir}/tmin_365_${year}.nc"
  echo "done tmin"

  # 5. PRECIP
  cdo divc,86400.0 \
    -daysum "${DataPath}/pre/crujra.v3.02.5d.pre.${year}.365d.noc.nc" \
    "${tmpdir}/precip_365_${year}.nc"
  echo "done precip"

  # 6. WIND SPEED
  cdo daymean \
    -expr,'wsspeed=sqrt(ugrd*ugrd+vgrd*vgrd)' \
    -merge \
      "${DataPath}/ugrd/crujra.v3.5d.ugrd.${year}.365d.noc.nc" \
      "${DataPath}/vgrd/crujra.v3.5d.vgrd.${year}.365d.noc.nc" \
    "${tmpdir}/wsspeed_365_${year}.nc"
  echo "done wsspeed"

  # 7. QAIR
  cdo daymean \
    "${DataPath}/spfh/crujra.v3.5d.spfh.${year}.365d.noc.nc" \
    "${tmpdir}/qair_365_${year}.nc"
  echo "done qair"

  # 8. Add Feb 29 for leap years (duplicate Feb 28 = timestep 59)
  vars=(longwave shortwave tmax tmin precip wsspeed qair)

  for var in "${vars[@]}"; do
    src="${tmpdir}/${var}_365_${year}.nc"
    dst="${tmpdir}/${var}_final_${year}.nc"
    if (( leap == 1 )); then
      cdo cat \
        -seltimestep,1/59  "$src" \
        -seltimestep,59    "$src" \
        -seltimestep,60/365 "$src" \
        "$dst"
      cdo settaxis,${year}-01-01,00:00:00,1day "$dst" "${dst%.nc}_fixed.nc"
      mv "${dst%.nc}_fixed.nc" "$dst"
    else
      cp "$src" "$dst"
    fi
  done

  # 9. Rename variables + set metadata
  cdo chname,dlwrf,longwave  "${tmpdir}/longwave_final_${year}.nc" "${tmpdir}/lw_r_${year}.nc"
  cdo chname,tswrf,shortwave "${tmpdir}/shortwave_final_${year}.nc" "${tmpdir}/sw_r_${year}.nc"
  cp "${tmpdir}/tmax_final_${year}.nc"    "${tmpdir}/tmax_r_${year}.nc"
  cp "${tmpdir}/tmin_final_${year}.nc"    "${tmpdir}/tmin_r_${year}.nc"
  cdo chname,pre,precip      "${tmpdir}/precip_final_${year}.nc"   "${tmpdir}/precip_r_${year}.nc"
  cp "${tmpdir}/wsspeed_final_${year}.nc" "${tmpdir}/wsspeed_r_${year}.nc"
  cdo chname,spfh,qair       "${tmpdir}/qair_final_${year}.nc"     "${tmpdir}/qair_r_${year}.nc"

  ncatted -O \
    -a standard_name,longwave,o,c,"surface_downwelling_longwave_flux_in_air" \
    -a long_name,longwave,o,c,"Surface Downwelling Longwave Radiation" \
    -a units,longwave,o,c,"W m-2" \
    "${tmpdir}/lw_r_${year}.nc"
  ncatted -O \
    -a standard_name,shortwave,o,c,"surface_downwelling_shortwave_flux_in_air" \
    -a long_name,shortwave,o,c,"Surface Downwelling Shortwave Radiation" \
    -a units,shortwave,o,c,"W m-2" \
    "${tmpdir}/sw_r_${year}.nc"
  ncatted -O \
    -a standard_name,tmax,o,c,"air_temperature" \
    -a long_name,tmax,o,c,"Daily Maximum Near-Surface Air Temperature" \
    -a units,tmax,o,c,"degC" \
    "${tmpdir}/tmax_r_${year}.nc"
  ncatted -O \
    -a standard_name,tmin,o,c,"air_temperature" \
    -a long_name,tmin,o,c,"Daily Minimum Near-Surface Air Temperature" \
    -a units,tmin,o,c,"degC" \
    "${tmpdir}/tmin_r_${year}.nc"
  ncatted -O \
    -a standard_name,precip,o,c,"precipitation_flux" \
    -a long_name,precip,o,c,"Precipitation" \
    -a units,precip,o,c,"kg m-2 s-1" \
    "${tmpdir}/precip_r_${year}.nc"
  ncatted -O \
    -a standard_name,wsspeed,o,c,"wind_speed" \
    -a long_name,wsspeed,o,c,"Near-Surface Wind Speed" \
    -a units,wsspeed,o,c,"m s-1" \
    "${tmpdir}/wsspeed_r_${year}.nc"
  ncatted -O \
    -a standard_name,qair,o,c,"specific_humidity" \
    -a long_name,qair,o,c,"Near-Surface Specific Humidity" \
    -a units,qair,o,c,"kg kg-1" \
    "${tmpdir}/qair_r_${year}.nc"

  # 10. Reorder longitude [-180,180] -> [0,360] and merge
  for var in lw sw tmax tmin precip wsspeed qair; do
    cdo sellonlatbox,0,360,-90,90 \
      "${tmpdir}/${var}_r_${year}.nc" \
      "${tmpdir}/${var}_lon360_${year}.nc"
  done

  cdo merge \
    "${tmpdir}/lw_lon360_${year}.nc" \
    "${tmpdir}/sw_lon360_${year}.nc" \
    "${tmpdir}/tmax_lon360_${year}.nc" \
    "${tmpdir}/tmin_lon360_${year}.nc" \
    "${tmpdir}/precip_lon360_${year}.nc" \
    "${tmpdir}/wsspeed_lon360_${year}.nc" \
    "${tmpdir}/qair_lon360_${year}.nc" \
    "${tmpdir}/climate_${year}_grid.nc"

  # 11. Remap to ICON grid
  for grid_tag in r2b4 r2b5; do
    if [[ "$grid_tag" == "r2b4" ]]; then
      grid_def="$grid_def_r2b4"
      grid_outdir="${outdir}/r2b4"
      institute_tag="mpim"
    else
      grid_def="$grid_def_r2b5"
      grid_outdir="${outdir}/r2b5"
      institute_tag="edzw"
    fi
    mkdir -p "$grid_outdir"

    cdo remapcon,"$grid_def" \
      "${tmpdir}/climate_${year}_grid.nc" \
      "${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_${year}.nc"

    echo "Year $year remapped -> ${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_${year}.nc"
  done

  # Cleanup
  rm -f "${tmpdir}"/*_${year}.nc

done 

# Fill 1850-1900 by recycling from 1901-1932
echo "Filling years 1850-1900 by recycling..."

for grid_tag in r2b4 r2b5; do
  if [[ "$grid_tag" == "r2b4" ]]; then
    grid_outdir="${outdir}/r2b4"
    institute_tag="mpim"
  else
    grid_outdir="${outdir}/r2b5"
    institute_tag="edzw"
  fi

  # 1900 is not a leap year but 1932 is -> use 1931 instead
  echo "Copying 1931 -> 1900"
  cdo setyear,1900 \
    "${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_1931.nc" \
    "${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_1900.nc"

  # Loop 1899 down to 1850
  for y in $(seq 1899 -1 1850); do
    a=$(( 1932 - y ))
    b=$(( a % 32 ))
    c=$(( 1932 - b ))
    echo "Copying $c -> $y"
    cdo setyear,${y} \
      "${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_${c}.nc" \
      "${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_${y}.nc"
  done

done  # end grid loop

echo "All years processed."
